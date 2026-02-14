
library(data.table)

algwas_names <-  fread("/media/bora_A/zhangt/2025-05-07-EAS_specific_xQTL-Project/2025-09-28-coloc/input/GWAS_qtype.txt")

GWAS_loci_info <- fread("2025-09-28-coloc/input/EAS_allsentinel.info.txt")
GWAS_loci_info <- unique(GWAS_loci_info)
colnames(GWAS_loci_info)[1] <- "GWAS_name"

GWAS_loci_info <- GWAS_loci_info[GWAS_name%in%algwas_names$GWAS_name]

table(GWAS_loci_info$GWAS_name%in% c("Height-Sakaue-2021", "Systemic@lupus@erythematosus-Wang-2021"))

GWAS_loci_info <- GWAS_loci_info[!GWAS_name%in% c("Height-Sakaue-2021","Systemic@lupus@erythematosus-Wang-2021")]
GWAS_loci_info$loci_start <- GWAS_loci_info$position - 1e6
GWAS_loci_info$loci_start[GWAS_loci_info$loci_start < 1] <- 1
GWAS_loci_info$loci_end   <- GWAS_loci_info$position + 1e6

GWAS_loci_info$overlap_MHC <- "FALSE"
GWAS_loci_info$overlap_MHC[GWAS_loci_info$chr == "chr6" & 
                             GWAS_loci_info$loci_start <= 33480577 & 
                             GWAS_loci_info$loci_end >= 28510120] <- "TRUE"

table(GWAS_loci_info$overlap_MHC)

GWAS_loci_info <- GWAS_loci_info[overlap_MHC=="FALSE", ]
GWAS_loci_info$overlap_MHC <- NULL

GWAS_loci_info <- GWAS_loci_info[chr%in%paste0("chr", 1:22)]
# 4505

fwrite(GWAS_loci_info, file = "2025-09-28-coloc/input/EAS_GWAS_noMHC_all_sentinal.txt", sep = "\t")
