#!/bin/bash

#==============================================#
# RNA_mapping_quantification #
#==============================================#

# >>> define functions <<<

function run_QC(){
	rawdir=path/to/sample
	baseDir=/path/to/basedir
	outDir=$baseDir/QC_RNA
	if [ ! -d "$outDir" ]
	then
		mkdir -p $outDir
	fi
	i=1
	for SAMPLE in `cat $baseDir/sample_list`
	do
		echo -e "${i}:\t$SAMPLE"
		fastqc $rawdir/${SAMPLE}/${SAMPLE}_1.fq.gz -o $outDir
		fastqc $rawdir/${SAMPLE}/${SAMPLE}_2.fq.gz -o $outDir
		i=$(( i + 1 ))
	done
}

function run_mapping(){
	dir=path/to/dir
	start_index=$dir/path/to/align_index
	ref_gtf=/path/to/gencode.v47.annotation.gtf
	BAM_dir=/path/to/Mapping_out
	module load bioinformatics/STAR-2.7.3a

# STAR-2.7.3a

	for SAMP in `cat $dir/task_use.txt`
	do
		echo $SAMP
		STAR --runThreadN 8 --genomeDir $start_index --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 20 --sjdbGTFfile $ref_gtf --limitOutSJcollapsed 1000000 --readFilesCommand zcat --outFileNamePrefix $BAM_dir/${SAMP}. --readFilesIn $dir/raw_fastq/${SAMP}_1.fq.gz $dir/raw_fastq/${SAMP}_2.fq.gz
	done

}

function run_rnaseqc(){
	dir=/path/to/dir
	bamDir=/path/to/output
	outDir=${dir}/RNA_seqc_out
	curDir=`pwd`
	i=1
	for SAMP in `cat $curDir/task_use.txt`
	do
		echo -e "$i\t$SAMP"
		$dir/bin/rnaseqc --sample=$SAMP --gene-length=100 --coverage --stranded="RF" $dir/path/to/gencode.v47.genes.gtf $bamDir/${SAMP}.Aligned.sortedByCoord.out.bam $outDir
		i=$(( i+1 ))
	done

}

