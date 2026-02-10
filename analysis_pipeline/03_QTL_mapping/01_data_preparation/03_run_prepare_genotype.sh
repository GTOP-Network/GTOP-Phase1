#!/bin/bash

main(){
	s1_make_plink_format
}

function s1_make_plink_format(){
	curr_dir=`pwd`
	VCF=/path/to/gtop_vcf_file.maf005.vcf.gz
	outdir=$curr_dir/output/genotype
	plink2 --output-chr chrM \
        	--vcf $VCF \
        	--vcf-half-call m \
        	--out $dir/GTOP.GT
}



# main
main
