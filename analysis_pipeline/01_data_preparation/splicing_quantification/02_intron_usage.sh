#!/bin/bash

########################### BACKGROUND ###########################
# Intron clustering performed using the protocol specified in    #
# leafcutter2                                                    #
# https://github.com/cfbuenabadn/leafcutter2.                    #
##################################################################

STAR_DIR=$1 #directory containing WASP pass bam
Meta_data=$2 Meta_data=$2  # The second positional argument passed to the script, which specifies a 
                           # metadata file containing two columns: 
                           # - Column 1: "tissue"
                           # - Column 2: "tissuecode"
                           # The file is typically tab-separated (TSV) without header
juncDir=$3 #directory containing junction files
intronDir=$4 #directory containing intron files
GTFfile=$5 #gtf file
reference=$6 #reference genome 


mkdir -p $juncDir $intronDir



#Estimate intron usage per sample (BAM) using regtools
for BAM in STAR_DIR
do
    #Variable prep
    library=`echo "${BAM}" | awk -F"/" '{print $(NF)}'`
    bamfile=$STAR_DIR/${BAM}


    #Generate regtools to quantify intron usage
    echo "Generating ${library}.junc"
    regtools junctions extract -s XS -a 8 -m 50 -M 500000 ${bamfile} -o $juncDir/$library.junc

done

#Estimate intron usage per tissues
cat $Meta_data|while read $tis $tissuecode;
do

mkdir -p $intronDir/$tis
#Write junc to txt file for clustering
ls $juncDir/*.junc|grep $tissuecode > $juncDir/$tis.junc_files.txt


#Cluster intron usage estimates using leafcutter2
echo "Generating intron usage file for ${tis} "
python ~/software/leafcutter2/scripts/leafcutter2.py \
-j $tis.junc_files.txt \
-r $intronDir/$tis \
-A ${GTFfile} \
-G ${reference} \
--checkchrom --keepleafcutter1 \
--minclureads 30 --mincluratio 0.001 --maxintronlen 500000

done