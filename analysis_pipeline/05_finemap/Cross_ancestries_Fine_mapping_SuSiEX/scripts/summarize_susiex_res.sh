#!/bin/bash

dir=`pwd`
target_dir=$1
Tissue=$2
outdir=$dir/output


cp $outdir/header.summary $outdir/${Tissue}.susiex.summary.txt
# summarize summary file
for f in `ls $target_dir/*.summary`
do
	fname=`basename $f .summary`
	gene=${fname#susiex_*}
	echo $gene
	N=`cat $f|grep -v "^#"| grep -v "NULL" | grep -v "CS_ID"|wc -l`
	if [ $N -gt 0 ]
	then
		cat $f|grep -v "^#"| grep -v "NULL" | grep -v "CS_ID"| awk -v GENE=$gene '{print $0"\t"GENE}' >> $outdir/${Tissue}.susiex.summary.txt
	fi
done

cp $outdir/header.cs $outdir/${Tissue}.susiex.cs.txt
# summarize cs file
for f in `ls $target_dir/*.cs`
do
	fname=`basename $f .cs`
	gene=${fname#susiex_*}
	echo $gene
	N=`cat $f|grep -v "NULL"|grep -v "CS_ID"|wc -l`
	if [ $N -gt 0 ]
	then
		cat $f|grep -v "CS_ID"|awk -v GENE=$gene '{print $0"\t"GENE}' >> $outdir/${Tissue}.susiex.cs.txt
	fi
done


# extract analyzed genes
for f in `ls $target_dir/*.summary`
do
	fname=`basename $f .summary`
	gene=${fname#susiex_*}
	echo $gene >> $outdir/${Tissue}.susiex_tested_genes.txt
done
