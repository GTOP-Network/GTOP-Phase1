#!/usr/bin/env Rscript

# Load required libraries
library(ieugwasr)
library(data.table)

# ------------------------------------------------------------
# Function: validate and parse command-line arguments
# Expected input order:
#   1. QTLTYPE (e.g., "SNV_eQTL")
#   2. TISSUENAME (e.g., "Pancreas_Body")
#   3. GENENAME (numeric, e.g., ENSG00000260196)
#   4. POPNAME (numeric, e.g., EAS)
# ------------------------------------------------------------
argvs <- commandArgs(trailingOnly = TRUE)

# Check number of arguments
if (length(argvs) != 3) {
  stop("Error: Exactly 3 command-line arguments are required: GWASNAME CHRNAME SENTINELPOS SENTINELSNP")
}

QTLTYPE <- argvs[1]
TISSUENAME <- argvs[2]
GENENAME <- argvs[3]

if(QTLTYPE %in% c("SNV_eQTL", "SNV_juQTL", "SNV_tuQTL", "SNV_pQTL")){
  POPNAME <- "EAS"
}else if(QTLTYPE %in% c("GTExv8_eQTL", "GTExv8_juQTL")){
  POPNAME <- "EUR"
}else if(QTLTYPE %in% c("JCTF_eQTL", "JCTF_juQTL")){
  POPNAME <- "EAS"
}else if(QTLTYPE %in% c("MAGE_eQTL", "MAGE_juQTL")){
  POPNAME <- "MAGE"
}else{
  stop("Error: input")
}

# ------------------------------------------------------------
# Read xQTL summary statistics
# ------------------------------------------------------------
QTL_data <- fread(sprintf("/flashfs1/scratch.global/tzhang/2025-06-11-specific_xQTL/data/tensorqtl/%s/nominal_slim/split/%s/%s.txt.gz",
                              QTLTYPE, TISSUENAME, GENENAME))
QTL_data$varbeta <- QTL_data$se * QTL_data$se
QTL_data$variant_id <- gsub("_b38", "", QTL_data$variant_id)
QTL_data$variant_id[QTL_data$variant_id==""] <- QTL_data$chr_pos_ref_alt[QTL_data$variant_id==""]
CHRNAME <- strsplit(QTL_data$chr_pos_ref_alt[1], "_")[[1]][1]
setDF(QTL_data)
rownames(QTL_data) <- QTL_data$variant_id

# ------------------------------------------------------------
# Compute LD matrix using 1000G EAS reference panel
# ------------------------------------------------------------
# bfile_path <- sprintf("/media/bora_A/zhangt/src/data/1000G/five_ancestry_groups/EAS/splitbychr/1000G.EAS.maf01.%s", CHRNAME)
# plink_bin  <- "/media/bora_A/zhangt/src/bin/plink"

bfile_path <- sprintf("/lustre/home/tzhang/src/1000G_%s_hg38/splitbychr/1000G.%s.maf01.%s", POPNAME, POPNAME, CHRNAME)
plink_bin  <- "/lustre/home/tzhang/src/plink"

# Verify that PLINK binary exists
if (!file.exists(plink_bin)) {
  stop("Error: PLINK binary not found at specified path.")
}

# Verify that PLINK bfile prefix exists (check .bed/.bim/.fam)
if (!all(file.exists(paste0(bfile_path, c(".bed", ".bim", ".fam"))))) {
  stop(sprintf("Error: PLINK reference files missing for chromosome: %s", bfile_path))
}

# Compute LD matrix
QTL_EAS_ld <- tryCatch({
  ld_matrix(
    variants = QTL_data$variant_id,
    bfile = bfile_path,
    plink_bin = plink_bin,
    with_alleles = FALSE
  )
}, error = function(e) {
  stop("Error during LD matrix computation: ", conditionMessage(e))
})

# Ensure LD matrix matches available SNPs
common_snps <- intersect(rownames(QTL_EAS_ld), rownames(QTL_data))
if (length(common_snps) == 0) {
  stop("Error: No overlapping SNPs between GWAS data and LD reference panel.")
}

# Subset GWAS data and LD matrix to common SNPs only
QTL_data_sub <- QTL_data[common_snps, , drop = FALSE]
QTL_EAS_ld_sub   <- QTL_EAS_ld[common_snps, common_snps, drop = FALSE]

# ------------------------------------------------------------
# Prepare dataset for coloc fine-mapping (SuSiE)
# ------------------------------------------------------------
d2 <- list(
  beta    = QTL_data_sub$beta,
  varbeta = QTL_data_sub$varbeta,
  pval    = QTL_data_sub$p_value,
  snp     = common_snps,
  LD      = as.matrix(QTL_EAS_ld_sub),
  N       = QTL_data_sub$N[1],
  type    = "quant",
  sdY     = 1
)

# ------------------------------------------------------------
# Run SuSiE fine-mapping with comprehensive error handling
# ------------------------------------------------------------
output_prefix <- sprintf(
  "input/finemapping_QTL/%s/%s/%s",
  QTLTYPE, TISSUENAME, GENENAME
)

# Ensure output directory exists
dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)

s2 <- tryCatch({
  runsusie(d2, maxit = 10000, repeat_until_convergence = FALSE)
}, error = function(e) {
  message("SuSiE failed to converge: ", conditionMessage(e))
  NULL
})

# Save result based on outcome
if (is.null(s2)) {
  # Save input data when SuSiE fails entirely
  save(d2, file = paste0(output_prefix, ".noConverged.RData"))
} else {
  susie_summary <- tryCatch(summary(s2), error = function(e) NULL)
  if (is.null(susie_summary) || is.null(susie_summary$cs)) {
    # Save input data when no credible sets are returned
    save(d2, file = paste0(output_prefix, ".noCS.RData"))
  } else {
    # Save successful SuSiE result
    save(s2, d2, file = paste0(output_prefix, ".finemapping.RData"))
  }
}

message("Processing completed for: ", TISSUENAME, " at ", GENENAME)
