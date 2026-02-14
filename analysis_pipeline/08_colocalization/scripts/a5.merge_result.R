#!/usr/bin/env Rscript

# Load required libraries
library(data.table)

# ------------------------------------------------------------
# Function: validate and parse command-line arguments
# Expected input order:
#   1. INPUTDIR (e.g., "slurm/coloc/task1/output")
#   2. OUTPUTDIR (e.g., "output")
# ------------------------------------------------------------
argvs <- commandArgs(trailingOnly = TRUE)

INPUTDIR <- argvs[1]
OUTPUTDIR <- argvs[2]

coloc_files <- list.files(INPUTDIR)
coloc_res <- rbindlist(lapply(coloc_files, function(x){
  f_info <-  fread(file.path(INPUTDIR, x))
  f_info
}))

for(i in unique(coloc_res$xQTL_type)){
  sub_data <- coloc_res[xQTL_type==i]
  sub_data <- sub_data %>% arrange(GWAS_name)
  fwrite(sub_data, file = paste0(OUTPUTDIR, "/", i, "_coloc.txt"), sep = "\t")
}
