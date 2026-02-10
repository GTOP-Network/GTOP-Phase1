#!/bin/bash


main(){

alignment
TR_calling
SNV_calling
SV_calling

}

#Function to align reads to human hg38 reference fasta
function alignment(){
    #Generate index file for reference (hg38)
    pbmm2 index reference.fasta reference.mmi --preset CCS
    #Align reads
    pbmm2 align --preset CCS --sort -j 32 -J 32 --log-level INFO reference.mmi ${input_bam} ${align_bam} 

}

# function to genotype the tandem repeat
function TR_calling(){

    # 1. Genotype the tandem repeat
    trgt genotype --genome reference.fasta --repeats ${repeat_bed} --reads ${align_bam} --output-prefix ${output_dir}/${sampleID}
    bcftools sort -Ob -o ${output_dir}/${sampleID}.sorted.vcf.gz ${output_dir}/${sampleID}.vcf.gz 
    bcftools index ${output_dir}/${sampleID}.sorted.vcf.gz --threads 8
    # 2. The trgt merge command can be used to merge TRGT VCFs into a single multi-sample VCF.
    trgt merge --vcf ${output_dir}/*.sorted.vcf.gz  \
		--genome reference.fasta \
		--output ${output_dir}/GTOP_LRS_TR.raw.vcf.gz \
		--output-type z --force-single 

    # 3. TR Variant Filtering 
    # 3.1 Filtering by missingness rate (max-missing 0.80)
    vcftools --gzvcf ${output_dir}/GTOP_LRS_TR.raw.vcf.gz --max-missing 0.80 --recode --recode-INFO-all --out ${output_dir}/GTOP_LRS_TR.raw.miss80
    # 3.2 removed TRs without major allele frequency <0.95 and allele count > 1 
    echo "Step 3.2a: Calculating AF and filtering by MAF <0.95 and AC > 1..."
    bcftools +fill-tags ${output_dir}/GTOP_LRS_TR.raw.miss80.recode.vcf -Oz -o ${output_dir}/GTOP_LRS_TR.raw.miss80.AF.vcf.gz --threads 32 -- -t AF 
    echo "Step 3.2b: Filtering by max-maf 0.95 and min-alleles 2..."
    vcftools --gzvcf ${output_dir}/GTOP_LRS_TR.raw.miss80.AF.vcf.gz --max-maf 0.95 --min-alleles 2 --recode --recode-INFO-all --stdout | bgzip -c > ${output_dir}/GTOP_LRS_TR.raw.miss80.AF95.AC1.vcf.gz
	bcftools view -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY \
    ${output_dir}/GTOP_LRS_TR.raw.miss80.AF95.AC1.vcf.gz -Oz -o ${output_dir}/GTOP_LRS_TR.raw.miss80.AF95.AC1.autosomeXY.vcf.gz


}

# function to genotype the structure variants
function SV_calling(){

    # Tool 1: sniffles2 (INS/DEL/DUP/TRA/BND discovery)
    # 1.1 Genotype the SV 
    sniffles --minsupport 4 --minsvlen 50 --threads 32 \
            --input ${align_bam} \
            --tandem-repeats human_GRCh38_no_alt_analysis_set.trf.bed \ 
            --reference reference.fasta \
            --output-rnames --vcf ${output_dir}/sniffles2/${sampleID}.vcf \
            --snf ${output_dir}/sniffles2/snf/${sampleID}.snf
    bcftools sort -Ob -o $sampleID.sorted.vcf.gz $sampleID.vcf.gz
    bcftools index $sampleID.sorted.vcf.gz
    # 1.2 merge multi-sample VCF 
    ls ${output_dir}/sniffles2/snf/*snf | awk -F'\' '{path = $0;  basename = $NF; split(basename, arr, ".");  print path "\t" arr[1];}' > GTOP_SV_snf.tsv
    sniffles  --input GTOP_SV_snf.tsv --reference reference.fasta --threads 32 --allow-overwrite --vcf ${output_dir}/GTOP_LRS_SV.sniffles.raw.vcf

    # Tool 2: pbsv (INS/DEL/DUP/TRA/BND discovery)
    # 2.1 Genotype the SV		
    pbsv discover --tandem-repeats ${trbed} ${align_bam}  ${output_dir}/SV_Calling/pbsv/${sampleID}.svsig.gz
    # 2.2 merge multi-sample VCF	
    pbsv call -j 32 --ccs -t '${TYPE}' -m 50  log reference.fasta ${output_dir}/SV_Calling/pbsv/${sampleID}.*.svsig.gz ${output_dir}/SV_Calling/pbsv/GTOP_LRS_SV.pbsv.vcf


    # Tool 3: cuteSV (INS/DEL/DUP/TRA/BND discovery)
    # 3.1 Genotype the SV		
    /lustre/home/rfding/anaconda3/envs/work/bin/cuteSV --genotype -l 50 -s 5 \
            --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 \
            --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5  \
            --min_size 50 --threads 32 \
            ${align_bam} \
            reference.fasta \
            -S ${sampleID} ${output_dir}/cuteSV/${sampleID}.vcf ${output_dir}/cuteSV/${sampleID}/temp
    ## 3.2 reGenotype and merge multi-sample VCF
    ls ${output_dir}/cuteSV/*.vcf > ${output_dir}/cuteSV/cuteSV.vcf_list
    /media/london_C/alps2/sunjs/Software/SURVIVOR/Debug/SURVIVOR merge ${output_dir}/cuteSV/cuteSV.vcf_list 500 1 1 -1 -1 -1 ${output_dir}/cuteSV/GTOP_LRS_SV.cuteSV.raw.vcf
	cuteFC ${align_bam} reference.fasta -S ${sampleID}  ${output_dir}/cuteSV/${sampleID}/temp \
        -Ivcf ${output_dir}/cuteSV/GTOP_LRS_SV.cuteSV.raw.vcf \
        -l 50 -s 5 \
        --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 \
        --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5 \
        --threads 32 
    ls *.reGenotype.sorted.vcf.gz > cuteSV_161INDs.merge_vcf_list.txt
    bcftools merge -m none --force-samples -l cuteSV_161INDs.merge_vcf_list.txt -Oz -o GTOP_LRS_SV.pbsv.merged.reGenotyped.vcf.gz

	# Tool 4: SVision ( Complex structure variants discovery)
    /lustre/home/rfding/anaconda3/envs/py3.6/bin/SVision -o ${output_dir}/${sampleID} -b ${align_bam} -m svision-cnn-model.ckpt -g reference.fasta -n ${sampleID} -t 32 -s 5 --graph --qname



}

# function to genotype the SNPs and small Indels(< 50bp)
function SNV_calling(){

    # 1. Genotype the SNV
    #docker pull google/deepvariant
    #make sure there is a index file for your hg38.fasta file, under your $INPUT_DIR directory. 
    INPUT_DIR=/media/iceland_B/share/Datasets/Asia_gtex/LRS_WGS_alignment
    OUTPUT_DIR=/media/london_B/lixing/2024-08-29-TR-AsianGTEX/2024-12-11-SNV-INDEL-LongReads/deepvariant
    INTERMEDIATE_DIRECTORY=/media/london_B/lixing/2024-08-29-TR-AsianGTEX/2024-12-11-SNV-INDEL-LongReads/deepvariant/intermediate_results_dir/${sampleID}
    BIN_VERSION="1.8.0"
    docker run -v "${INPUT_DIR}:${INPUT_DIR}" \
                -v "${OUTPUT_DIR}:${OUTPUT_DIR}" \
                google/deepvariant:"${BIN_VERSION}" \
                /opt/deepvariant/bin/run_deepvariant \
                --model_type PACBIO \
                --ref  reference.fasta \
                --reads ${sampleID} \
                --output_vcf "${OUTPUT_DIR}/${sampleID}.vcf.gz" \
                --output_gvcf "${OUTPUT_DIR}/${sampleID}.g.vcf.gz" \
                --num_shards 32 \
                --intermediate_results_dir "${INTERMEDIATE_DIRECTORY}"
    
    # 2. Merge multi-sample VCF
    /media/london_B/lixing/software/glnexus_cli --config DeepVariantWGS ${OUTPUT_DIR}/*.g.vcf.gz > ${OUTPUT_DIR}/GTOP_LRS_SNV.raw.bcf
    bcftools view ${OUTPUT_DIR}/GTOP_LRS_SNV.raw.bcf | bgzip -@ 24 -c > ${OUTPUT_DIR}/GTOP_LRS_SNV.raw.vcf.gz
    tabix -p vcf ${OUTPUT_DIR}/GTOP_LRS_SNV.raw.vcf.gz
    bcftools norm -m -both -d none -freference.fasta ${OUTPUT_DIR}/GTOP_LRS_SNV.raw.vcf.gz | bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%ALT' -Oz -o ${OUTPUT_DIR}/GTOP_LRS_SNV.raw.bialleic.vcf.gz
    # 3. Filtering
    # 3.1 QUAL>=20 && F_MISSING<=0.15 && abs(strlen(REF)-strlen(ALT))<=50
    bcftools +fill-tags ${OUTPUT_DIR}/GTOP_LRS_SNV.raw.bialleic.vcf.gz  -Oz -o ${OUTPUT_DIR}/GTOP_LRS_SNV.raw.bialleic.annotation.vcf.gz -- -t all 
    tabix -p vcf ${OUTPUT_DIR}/GTOP_LRS_SNV.raw.bialleic.annotation.vcf.gz
    bcftools view -i 'QUAL>=20 && F_MISSING<=0.15 && abs(strlen(REF)-strlen(ALT))<=50' ${OUTPUT_DIR}/GTOP_LRS_SNV.raw.bialleic.annotation.vcf.gz -Oz -o ${OUTPUT_DIR}/GTOP_LRS_SNV.raw.bialleic.Q20.missing15.len50.vcf.gz
    tabix -p vcf ${OUTPUT_DIR}/GTOP_LRS_SNV.raw.bialleic.Q20.missing15.len50.vcf.gz
    
    # 3.2 hwe 1e-6
    vcftools --gzvcf ${OUTPUT_DIR}/GTOP_LRS_SNV.raw.bialleic.Q20.missing15.len50.vcf.gz --hwe 1e-6 --recode --recode-INFO-all --stdout | bgzip -c > ${OUTPUT_DIR}/GTOP_LRS_SNV.raw.bialleic.Q20.missing15.len50.hwe.vcf.gz
    tabix -p vcf  ${OUTPUT_DIR}/GTOP_LRS_SNV.raw.bialleic.Q20.missing15.len50.hwe.vcf.gz

    # 3.3 RSID annotation
    bcftools annotate -a GCF_000001405.40.vcf.gz -c ID  ${OUTPUT_DIR}/GTOP_LRS_SNV.raw.bialleic.Q20.missing15.len50.hwe.vcf.gz -Oz -o   ${OUTPUT_DIR}/GTOP_LRS_SNV.filtered.vcf.gz --threads 32 && tabix -p vcf ${OUTPUT_DIR}/GTOP_LRS_SNV.filtered.vcf.gz

}

main