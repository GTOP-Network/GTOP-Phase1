#!/usr/bin/env Rscript

# Load required libraries
library(ieugwasr)
library(data.table)

# ------------------------------------------------------------
# Function: validate and parse command-line arguments
# Expected input order:
#   1. GWASNAME (e.g., "250.2")
#   2. CHRNAME (e.g., "chr11")
#   3. SENTINELPOS (numeric, e.g., 17387083)
#   4. SENTINELSNP (e.g., "rs5215" or "chr11_17387083_A_T")
# ------------------------------------------------------------
argvs <- commandArgs(trailingOnly = TRUE)

# Check number of arguments
if (length(argvs) != 4) {
  stop("Error: Exactly 4 command-line arguments are required: GWASNAME CHRNAME SENTINELPOS SENTINELSNP")
}

GWASNAME      <- as.character(argvs[1])
CHRNAME       <- as.character(argvs[2])
SENTINELPOS   <- as.numeric(argvs[3])
SENTINELSNP   <- as.character(argvs[4])

# Validate numeric conversion
if (is.na(SENTINELPOS)) {
  stop("Error: SENTINELPOS must be a valid numeric position.")
}

# Define locus window: ±1 Mb around sentinel SNP
LOCISTART <- max(1, SENTINELPOS - 1e6)
LOCIEND   <- SENTINELPOS + 1e6

# ------------------------------------------------------------
# Load GWAS trait type from metadata file
# ------------------------------------------------------------
gwas_qtype_file <- "input/GWAS_qtype.txt"
if (!file.exists(gwas_qtype_file)) {
  stop("Error: GWAS_qtype.txt file not found at expected path.")
}

GWAS_qtype_df <- fread(gwas_qtype_file, stringsAsFactors = FALSE)
if (!GWASNAME %in% GWAS_qtype_df$GWAS_name) {
  stop(sprintf("Error: GWAS name '%s' not found in GWAS_qtype.txt.", GWASNAME))
}
GWAS_qtype <- GWAS_qtype_df[GWAS_name == GWASNAME, Type]

# ------------------------------------------------------------
# Read GWAS summary statistics for the specified chromosome
# ------------------------------------------------------------
gwas_input_path <- sprintf("input/susie_coloc_split/%s/%s_coloc.%s.txt", 
                           GWASNAME, GWASNAME, CHRNAME)
if (!file.exists(gwas_input_path)) {
  stop(sprintf("Error: GWAS input file not found: %s", gwas_input_path))
}

GWAS_data <- fread(gwas_input_path, stringsAsFactors = FALSE)

# Subset to locus window and remove rows with missing standard errors
GWASloci_data <- GWAS_data[position >= LOCISTART & position <= LOCIEND]
GWASloci_data <- GWASloci_data[!is.na(se)]

# Ensure required columns exist
required_cols <- c("position", "se", "beta", "maf", "P", "rsid", "N")
missing_cols <- setdiff(required_cols, names(GWASloci_data))
if (length(missing_cols) > 0) {
  stop(sprintf("Error: Missing required columns in GWAS data: %s", paste(missing_cols, collapse = ", ")))
}

# Compute variance of beta (needed for coloc)
GWASloci_data[, varbeta := se^2]

# Standardize rsid format: keep rsIDs as-is; convert others to chr_pos format
GWASloci_data[, rsid_clean := ifelse(
  grepl("^rs\\d+", rsid), 
  rsid, 
  paste0("chr", gsub(":", "_", gsub("^chr", "", rsid, ignore.case = TRUE)))
)]

# Convert to data.frame for compatibility with coloc/ieugwasr
setDF(GWASloci_data)
GWASloci_data <- GWASloci_data[!duplicated(GWASloci_data$rsid_clean),]
rownames(GWASloci_data) <- GWASloci_data$rsid_clean

# ------------------------------------------------------------
# Compute LD matrix using 1000G EAS reference panel
# ------------------------------------------------------------
# bfile_path <- sprintf("/media/bora_A/zhangt/src/data/1000G/five_ancestry_groups/EAS/splitbychr/1000G.EAS.maf01.%s", CHRNAME)
# plink_bin  <- "/media/bora_A/zhangt/src/bin/plink"

bfile_path <- sprintf("/lustre/home/tzhang/src/1000G_EAS_hg38/splitbychr/1000G.EAS.maf01.%s", CHRNAME)
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
GWAS_EAS_ld <- tryCatch({
  ld_matrix(
    variants = GWASloci_data$rsid_clean,
    bfile = bfile_path,
    plink_bin = plink_bin,
    with_alleles = FALSE
  )
}, error = function(e) {
  stop("Error during LD matrix computation: ", conditionMessage(e))
})

# Ensure LD matrix matches available SNPs
common_snps <- intersect(rownames(GWAS_EAS_ld), rownames(GWASloci_data))
if (length(common_snps) == 0) {
  stop("Error: No overlapping SNPs between GWAS data and LD reference panel.")
}

# Subset GWAS data and LD matrix to common SNPs only
GWASloci_data_sub <- GWASloci_data[common_snps, , drop = FALSE]
GWAS_EAS_ld_sub   <- GWAS_EAS_ld[common_snps, common_snps, drop = FALSE]

# ------------------------------------------------------------
# Prepare dataset for coloc fine-mapping (SuSiE)
# ------------------------------------------------------------
d1 <- list(
  beta    = GWASloci_data_sub$beta,
  varbeta = GWASloci_data_sub$varbeta,
  MAF     = GWASloci_data_sub$maf,
  pval    = GWASloci_data_sub$P,
  snp     = rownames(GWASloci_data_sub),
  LD      = as.matrix(GWAS_EAS_ld_sub),
  N       = GWASloci_data_sub$N[1],
  type    = GWAS_qtype
)

# ------------------------------------------------------------
# Run SuSiE fine-mapping with comprehensive error handling
# ------------------------------------------------------------
output_prefix <- sprintf(
  "input/finemapping_GWAS/%s_%s",
  GWASNAME, SENTINELSNP
)

# Ensure output directory exists
dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)

s1 <- tryCatch({
  runsusie(d1, maxit = 10000, repeat_until_convergence = FALSE)
}, error = function(e) {
  message("SuSiE failed to converge: ", conditionMessage(e))
  NULL
})

# Save result based on outcome
if (is.null(s1)) {
  # Save input data when SuSiE fails entirely
  save(d1, file = paste0(output_prefix, ".noConverged.RData"))
} else {
  susie_summary <- tryCatch(summary(s1), error = function(e) NULL)
  
  if (is.null(susie_summary) || is.null(susie_summary$cs)) {
    # Save input data when no credible sets are returned
    save(d1, file = paste0(output_prefix, ".noCS.RData"))
  } else {
    # Save successful SuSiE result
    save(s1, d1, file = paste0(output_prefix, ".finemapping.RData"))
  }
}

message("Processing completed for: ", GWASNAME, " at ", SENTINELSNP)
