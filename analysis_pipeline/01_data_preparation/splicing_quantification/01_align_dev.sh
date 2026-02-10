#!/bin/bash

#==========#
# Get args #
#==========#

workDir=$1     # Working directory where outputs will be generated
VCFdir=$2   # Directory containing per doner VCF files
starIndex=$3   # Path to the STAR index directory
fastqDir=$4    # Path to directory containing FASTQ files, organized into subdirectories based on sequencing batch
outDir=$5    # Path to directory containing outputfiles
mkdir -p $workDir/output/wasp_mapping

allsample=(` ls $fastqDir/ |cut -d_ -f1|sort|uniq`)

for sample in ${allsample[@]};
do
INDS=`echo $sample|cut -d- -f2`


## align RNA-seq data to reference with STAR

STAR \
        --runMode alignReads \
        --runThreadN 32 \
        --genomeDir $starIndex \
        --readFilesIn $fastqDir/${sample}_1.fq.gz $fastqDir/${sample}_2.fq.gz \
        --readFilesCommand zcat \
        --outFileNamePrefix $outDir/$sample. \
        --twopassMode Basic \
        --outFilterMultimapNmax 20 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverLmax 0.1 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --outFilterType BySJout \
        --outFilterScoreMinOverLread 0.33 \
        --outFilterMatchNmin 0 \
        --outFilterMatchNminOverLread 0.33 \
        --limitSjdbInsertNsj 1200000 \
        --limitOutSJcollapsed 5000000 \
        --outSAMstrandField intronMotif \
        --outFilterIntronMotifs None \
        --alignSoftClipAtReferenceEnds Yes \
        --outSAMtype BAM Unsorted \
        --outSAMunmapped Within \
        --outSAMattrRGline ID:rg1 SM:sm1 \
        --outSAMattributes NH HI AS nM NM ch vW \
        --waspOutputMode SAMtag \
        --varVCFfile $VCFdir/$INDS.vcf.gz \
        --winAnchorMultimapNmax 50 \
        --chimSegmentMin 15 \
        --chimJunctionOverhangMin 15 \
        --chimOutType Junctions WithinBAM SoftClip \
        --chimMainSegmentMultNmax 1 \
        --chimOutJunctionFormat 0


#======================#
# Keep WASP pass reads #
#======================#

samtools view -h -q 255 ${outDir}/$sample.Aligned.sortedByCoord.out.bam | grep -v "vW:i:[2-7]" | samtools view -h -b > ${outDir}/$sample.Aligned.sortedByCoord.WASPpass.out.bam
samtools index ${outDir}/$sample.Aligned.sortedByCoord.WASPpass.out.bam

done


