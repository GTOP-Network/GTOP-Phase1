#!/bin/bash

vcfFile=$1 # Separate genotype VCF files were generated for each variant type, including small variants, SVs, and TRs. 
outDir=$2 # Directory to write output files to

vep --input_file ${vcfFile} \
  --output_file ${outDir}/GTOP.${variant_type}.vep_anno.vcf \
  --format vcf \
  --species homo_sapiens \
  --cache \
  --offline \
  --dir_cache ensembl-vep/GRCH38 \
  --assembly GRCh38 \
  --plugin LoF,loftee_path:loftee,human_ancestor_fa:human_ancestor.fa.gz,conservation_file:loftee.sql \
  --dir_pluginsloftee \
  --tab 
  
vep \
  --input_file ${vcfFile}  \
  --output_file ${outDir}/GTOP.${variant_type}.vep_Regulatory_anno.vcf \
  --regulatory \
  --custom homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20240230.chr.bed.gz,Regulatory_Build,BED,overlap,0 \
  --custom Homo_sapiens.GRCh38.motif_features.v113.bed.gz,TF,BED,overlap,0 


