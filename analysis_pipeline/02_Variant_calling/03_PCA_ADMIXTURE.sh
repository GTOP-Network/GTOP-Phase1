#!/bin/bash


CUR_DIR="path/to/work/dir"
GTOP_VCF='path/to/GTOP/vcf/dir' # The path to the GTOP small variant VCF file
DIR_VCF_1KG='path/to/1KGP/vcf/dir' # The directory containing 1000 Genomes VCF files (one per chromosome) from http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage
DIR_VCF_GTEx='path/to/GTEX/vcf/dir' # The path to the GTEx VCF file


mkdir -p ${CUR_DIR}/input/101_vcf_chrs
for chr in {1..22} X Y
    do
         bcftools view -r chr${chr} ${GTOP_VCF} -Oz -o ${CUR_DIR}/input/101_vcf_chrs/GTOP_LRS_SmallVariant.chr${chr}.vcf.gz && tabix -p vcf  ${CUR_DIR}/input/101_vcf_chrs/GTOP_LRS_SmallVariant.chr${chr}.vcf.gz &
done   


###
mkdir -p ${CUR_DIR}/input/102_vcf_chrs_1KG

for chr in {1..22}
    do
        bcftools merge ${CUR_DIR}/input/101_vcf_chrs/GTOP_LRS_SmallVariant.chr${chr}.vcf.gz ${DIR_VCF_1KG}/1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz ${DIR_VCF_GTEx}/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.chr${chr}.vcf.gz \\
        -Oz -o ${CUR_DIR}/input/102_vcf_chrs_1KG/GTOP_1KG_GTEx.SNP_INDEL.chr${chr}.vcf.gz --threads 4 && tabix -p vcf  ${CUR_DIR}/input/102_vcf_chrs_1KG/GTOP_1KG_GTEx.SNP_INDEL.chr${chr}.vcf.gz &
done

###
mkdir -p ${CUR_DIR}/input/103_intertect_snps
for chr in {1..22}
    do
        bcftools isec -n=3 -c all ${CUR_DIR}/input/101_vcf_chrs/GTOP_LRS_SmallVariant.chr${chr}.vcf.gz   \
        ${DIR_VCF_1KG}/1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz  \
        ${DIR_VCF_GTEx}/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.chr${chr}.vcf.gz   \
        -p  ${CUR_DIR}/input/103_intertect_snps/GTOP_1KG_GTEx.intersect.chr${chr} --threads 4 && echo "done for chr${chr}" &
done


###
mkdir -p ${CUR_DIR}/input/104_vcf_chrs_1KG
for chr in {1..22}
    do
        bcftools view -R ${CUR_DIR}/input/103_intertect_snps/GTOP_1KG_GTEx.intersect.chr${chr}/sites.txt ${CUR_DIR}/input/102_vcf_chrs_1KG/GTOP_1KG_GTEx.SNP_INDEL.chr${chr}.vcf.gz | bcftools norm -m +any - | bcftools view -m 2 -M 2 -v snps - | bcftools annotate -x INFO,^FORMAT/GT,QUAL,FILTER - -Oz -o ${CUR_DIR}/input/104_vcf_chrs_1KG/GTOP_1KG_GTEx.SNP.chr${chr}.vcf.gz &&   tabix -p vcf  ${CUR_DIR}/input/104_vcf_chrs_1KG/GTOP_1KG_GTEx.SNP.chr${chr}.vcf.gz &
done


###
mkdir -p ${CUR_DIR}/input/105_vcf
line=''
for chr in {1..22}
    do
        line=${line}' '${CUR_DIR}/input/104_vcf_chrs_1KG/GTOP_1KG_GTEx.SNP.chr${chr}.vcf.gz
done
 bcftools concat $line -Oz -o ${CUR_DIR}/input/105_vcf/GTOP_1KG_GTEx.SNP.chrs.vcf.gz && tabix -p vcf ${CUR_DIR}/input/105_vcf/GTOP_1KG_GTEx.SNP.chrs.vcf.gz


###
bcftools view -S ${CUR_DIR}/input/GTOP_1KG_unrelated_GTEx.population.list.txt   ${CUR_DIR}/input/105_vcf/GTOP_1KG_GTEx.SNP.chrs.vcf.gz  -Oz -o ${CUR_DIR}/input/106_1KG_unrel/GTOP_1KG_GTEx.SNP.chrs.vcf.gz  && tabix -p vcf ${CUR_DIR}/input/106_1KG_unrel/GTOP_1KG_GTEx.SNP.chrs.vcf.gz   
plink --vcf  ${CUR_DIR}/input/106_1KG_unrel/GTOP_1KG_GTEx.SNP.chrs.vcf.gz --maf 0.05 --bp-space 100000 --keep-allele-order --recode vcf-iid bgz --out  ${CUR_DIR}/input/106_1KG_unrel/GTOP_1KG_GTEx.SNP.para_100K
mkdir -p input/107_PCA
QTLtools pca --vcf ${CUR_DIR}/input/106_1KG_unrel/GTOP_1KG_GTEx.SNP.para_100K.vcf.gz --center --scale --maf 0.05 --distance 100000 --out  ${CUR_DIR}/input/107_PCA/GTOP_1KG_GTEx.SNP.para_100K.SNP

###
mkdir -p ${CUR_DIR}/input/109_PCA_EastAsia_1KG_GTEx
bcftools view -S ${CUR_DIR}/input/GTOP_1KG_GTEx_EAS.population.list.txt  ${CUR_DIR}/input/106_1KG_unrel/GTOP_1KG_GTEx.SNP.para_100K.vcf.gz -Oz -o  ${CUR_DIR}/input/109_PCA_EastAsia_1KG_GTEx/GTOP_1KG_GTEx_EastAsian.SNP.para_100K.vcf.gz  && tabix -p vcf  ${CUR_DIR}/input/109_PCA_EastAsia_1KG_GTEx/GTOP_1KG_GTEx_EastAsian.SNP.para_100K.vcf.gz
QTLtools pca --vcf  ${CUR_DIR}/input/109_PCA_EastAsia_1KG_GTEx/GTOP_1KG_GTEx_EastAsian.SNP.para_100K.vcf.gz --center --scale --maf 0.05 --distance 100000 --out   ${CUR_DIR}/input/109_PCA_EastAsia_1KG_GTEx/GTOP_1KG_GTEx_EastAsian.SNP.para_100K.qtltools


###
mkdir -p ${CUR_DIR}/input/108_admixture
plink --vcf  ${CUR_DIR}/input/106_1KG_unrel/cohort_1KG_GTEx.SNP.para_100K.vcf.gz --make-bed --recode --out  ${CUR_DIR}/input/108_admixture/cohort_1KG_GTEx.SNP.para_100K
plink --bfile  ${CUR_DIR}/input/108_admixture/cohort_1KG_GTEx.SNP.para_100K --pca --out  ${CUR_DIR}/input/108_admixture/cohort_1KG_GTEx.SNP.para_100K
mkdir -p ${CUR_DIR}/input/108_admixture/repeat1
for K in {2..10}
    do
        admixture --cv ${CUR_DIR}/input/108_admixture/cohort_1KG_GTEx.SNP.para_100K.bed $K -s time  | tee   ${CUR_DIR}/input/108_admixture/repeat1/log${K}.out
done
    




