
setwd("/media/bora_A/zhangt/2025-05-07-EAS_specific_xQTL-Project/")

library(data.table)
library(dplyr)
library(tidyr)

# "Gallbladder"
tissue_change <- c("Adipose"="Adipose_Visceral_Omentum",
                   "Adrenal_Gland"="Adrenal_Gland",
                   "Liver"="Liver",
                   "Muscle"="Muscle_Skeletal",
                   "Pancreas_Body"="Pancreas",
                   "Pancreas_Head"="Pancreas",
                   "Pancreas_Tail"="Pancreas",
                   "Skin"="Skin_Not_Sun_Exposed_Suprapubic",
                   "Spleen"="Spleen",
                   "Whole_Blood"="Whole_Blood")

xQTL_type <- "eQTL"

GTOP_variants <- vroom::vroom("2025-06-11-specific_xQTL/input/GTOP/SNV/GTOP_chrpos_rsid.txt", col_names = FALSE) %>% as.data.frame()
colnames(GTOP_variants) <- c("variant_id", "rsid")


GTEx_gene_info <- fread("2025-06-11-specific_xQTL/data/tensorqtl/GTExv8_eQTL/gencode.v26.gene.loci")
GTOP_gene_info <- fread("2025-06-11-specific_xQTL/input/GTOP/annotation/GTOP_gencode_v47.gene_info.txt")
GTEx_genes <- unique(gsub("\\.[0-9]+", "", GTEx_gene_info$V1))
GTOP_genes <- unique(gsub("\\.[0-9]+", "", GTOP_gene_info$V5))
overlapped_genes <- intersect(GTEx_genes, GTOP_genes)


######################
## permutation -------
######################

## GTEx ---------
GTEx_permutation_dir <- sprintf("2025-06-11-specific_xQTL/data/tensorqtl/GTExv8_%s/permutation", xQTL_type)

pbmcapply::pbmclapply(unique(tissue_change), function(tmp_tissue){
  cat(tmp_tissue, "\n")
  tmp_data <- vroom::vroom(paste0(GTEx_permutation_dir, "/", tmp_tissue, ".txt.gz")) %>% 
                             as.data.frame()
  tmp_data$variant_id <- gsub("_b38", "", tmp_data$variant_id)
  if(xQTL_type=="eQTL"){
    tmp_data$phenotype_id <- gsub("\\.\\d+", "", tmp_data$gene_id)
    tmp_data <- tmp_data[tmp_data$phenotype_id%in%overlapped_genes, ]
  }else{
    tmp_data$phenotype_id <- gsub("\\.\\d+", "", tmp_data$phenotype_id)
    tmp_data$phenotype_id <- gsub(":clu_[0-9]+", "", tmp_data$phenotype_id)
  }
  tmp_data$rsid <- tmp_data$rs_id_dbSNP151_GRCh38p7
  tmp_data <- tmp_data[, c("variant_id", "phenotype_id", "pval_nominal", "slope", "slope_se")]
  fwrite(tmp_data, file = sprintf("2025-09-30-mash/output/SNV_%s/permutation/GTEx_%s.txt", xQTL_type, tmp_tissue), sep="\t")
}, mc.cores = 2, mc.preschedule = FALSE)


## GTOP
GTOP_permutation_dir <- sprintf("2025-06-11-specific_xQTL/data/tensorqtl/SNV_%s/permutation/", xQTL_type)

pbmcapply::pbmclapply(names(tissue_change), function(tmp_tissue){
  cat(tmp_tissue, "\n")
  tmp_data <- vroom::vroom(paste0(GTOP_permutation_dir, "/", tmp_tissue, ".txt.gz")) %>% 
    as.data.frame()
  tmp_data$rsid <- tmp_data$variant_id
  if(xQTL_type == "eQTL"){
    tmp_data$phenotype_id <- gsub("\\.\\d+", "", tmp_data$phenotype_id)
    tmp_data <- tmp_data[tmp_data$phenotype_id%in%overlapped_genes, ]
  }else{
    tmp_data$phenotype_id <- gsub("\\.\\d+", "", tmp_data$phenotype_id)
    tmp_data$phenotype_id <- gsub(":clu_[0-9]+", "", tmp_data$phenotype_id)
  }

  tmp_data$variant_id <- NULL
  tmp_data <- inner_join(tmp_data, GTOP_variants, by="rsid")
  
  tmp_data <- tmp_data[, c("variant_id", "phenotype_id", "pval_nominal", "slope", "slope_se")]
  fwrite(tmp_data, file = sprintf("2025-09-30-mash/output/SNV_%s/permutation/GTOP_%s.txt", xQTL_type, tmp_tissue), sep="\t")
}, mc.cores = 5, mc.preschedule = FALSE)


##################################
## Nominal -----------------------
##################################

## GTEx ---------
GTEx_nominal_dir <- sprintf("2025-06-11-specific_xQTL/data/tensorqtl/GTExv8_%s/nominal", xQTL_type)

pbmcapply::pbmclapply(unique(tissue_change)[c(8)], function(tmp_tissue){
  cat(tmp_tissue, "\n")
  tmp_data <- vroom::vroom(paste0(GTEx_nominal_dir, "/", tmp_tissue, ".txt.gz")) %>% as.data.frame()
  
  tmp_data$variant_id <- gsub("_b38", "", tmp_data$variant_id)
  if(xQTL_type=="eQTL"){
    tmp_data$phenotype_id <- gsub("\\.\\d+", "", tmp_data$gene_id)
    tmp_data <- tmp_data[tmp_data$phenotype_id%in%overlapped_genes, ]
  }else{
    tmp_data$phenotype_id <- gsub("\\.\\d+", "", tmp_data$phenotype_id)
    tmp_data$phenotype_id <- gsub(":clu_[0-9]+", "", tmp_data$phenotype_id)
  }

  tmp_data <- tmp_data[, c("variant_id", "phenotype_id", "pval_nominal", "slope", "slope_se")]
  fwrite(tmp_data, file = sprintf("2025-09-30-mash/output/SNV_%s/nominal/GTEx_%s.txt.gz", 
                                  xQTL_type, tmp_tissue), sep="\t")
}, mc.cores = 1, mc.preschedule = FALSE)


## GTOP ---------
GTOP_nominal_dir <- sprintf("2025-06-11-specific_xQTL/data/tensorqtl/SNV_%s/nominal/", xQTL_type)

pbmcapply::pbmclapply(names(tissue_change), function(tmp_tissue){

  cat(tmp_tissue, "\n")
  tmp_data <- fread(paste0(GTOP_nominal_dir, "/", tmp_tissue, ".txt.gz")) %>% as.data.frame()
  
  tmp_data$rsid <- tmp_data$variant_id
  if(xQTL_type == "eQTL"){
    tmp_data$phenotype_id <- gsub("\\.\\d+", "", tmp_data$phenotype_id)
    tmp_data <- tmp_data[tmp_data$phenotype_id%in%overlapped_genes, ]
  }else{
    tmp_data$phenotype_id <- gsub("\\.\\d+", "", tmp_data$phenotype_id)
    tmp_data$phenotype_id <- gsub(":clu_[0-9]+", "", tmp_data$phenotype_id)
  }
  
  tmp_data$variant_id <- NULL
  tmp_data <- inner_join(tmp_data, GTOP_variants, by="rsid")
  
  tmp_data <- tmp_data[, c("variant_id", "phenotype_id", "pval_nominal", "slope", "slope_se")]
  
  fwrite(tmp_data, file = sprintf("2025-09-30-mash/output/SNV_%s/nominal/GTOP_%s.txt.gz", xQTL_type, tmp_tissue), sep="\t")
}, mc.cores = 3, mc.preschedule = FALSE)
