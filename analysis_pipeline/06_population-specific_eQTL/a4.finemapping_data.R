
library(data.table)
library(tidyverse)

GTEx_annotation <- fread("2025-06-11-specific_xQTL/data/tensorqtl/GTExv8_eQTL/gencode.v26.gene.loci")
GTOP_annotation <- fread("2025-06-11-specific_xQTL/input/GTOP/annotation/GTOP_gencode_v47.gene_info.txt")
GTEx_annotation$gene_tmpid <- gsub("\\..+", "", GTEx_annotation$V1)
GTOP_annotation$gene_tmpid <- gsub("\\..+", "", GTOP_annotation$V5)
overlapped_genes <- intersect(GTOP_annotation$gene_tmpid, GTEx_annotation$gene_tmpid)


## xQTL type
xQTL_type <- "eQTL"

eQTL_fm_freq <- fread("2025-06-11-specific_xQTL/output/data/eQTL/finemapping_add_GTEx_freq_all.txt.gz")
eQTL_fm_id <- unique(eQTL_fm_freq[, .(variant_id=chr_pos_ref_alt,
                                      phenotype_id=gene_name)])
eQTL_fm_id$phenotype_id <- case_when(
  grepl("^ENSG", eQTL_fm_id$phenotype_id) ~ gsub("\\.[0-9]+", "", eQTL_fm_id$phenotype_id),
  TRUE ~ gsub("\\.[0-9]+", "", eQTL_fm_id$phenotype_id)
)

eQTL_fm_ol_id <- eQTL_fm_id[eQTL_fm_id$phenotype_id%in%overlapped_genes, ]


## mash files
mash_nominal_dir <- "2025-09-30-mash/output/SNV_eQTL/nominal/" 
tissue_list <- list.files(mash_nominal_dir)

pbmcapply::pbmclapply(tissue_list, function(tmp_file){
  
  cat(tmp_file, "\n")
  
  tmp_data <- fread(paste0(mash_nominal_dir, "/", tmp_file))
  
  tmp_subdata <- merge(tmp_data, eQTL_fm_ol_id, by=c("variant_id", "phenotype_id"))
  
  fwrite(tmp_subdata, file = sprintf("2025-09-30-mash/output/SNV_%s/finemapping/%s", xQTL_type, tmp_file), sep="\t")
}, mc.cores = 4, mc.preschedule = FALSE)

