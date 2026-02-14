#!/usr/bin/env Rscript

setwd("/media/bora_A/zhangt/2025-05-07-EAS_specific_xQTL-Project")

# Load required libraries
library(data.table)
library(tidyverse)

## load data
SNV_freq_1KG_res <- fread("2025-06-11-specific_xQTL/output/data/processing/GTOP_1KG_SNV_freq.txt.gz")
table(abs(SNV_freq_1KG_res$af_GTOP-SNV_freq_1KG_res$af_EAS)<0.2)


GTOP_freq <- fread("2025-06-11-specific_xQTL/output/data/processing/GTOP_SNV_information.txt.gz")
GTEx_freq <- fread("2025-06-11-specific_xQTL/input/GTEx/GTEx_EUR.afreq")
GTEx_freq$chr_pos_ref_alt <- gsub("_b38", "", GTEx_freq$ID)

GTOP_GTEx_freq1 <- merge(GTOP_freq, GTEx_freq[, .(chr_pos_ref_alt, af_GTEx_raw = ALT_FREQS)], by="chr_pos_ref_alt", all.x=TRUE)
GTOP_GTEx_freq <- merge(GTOP_GTEx_freq1, SNV_freq_1KG_res[, .(chr_pos_ref_alt, af_EUR)], by="chr_pos_ref_alt", all.x=TRUE)

# NFE_freq <- fread("../src/data/1000G/five_ancestry_groups/NFE/NFE.afreq")


GTOP_GTEx_freq$af_EUR <- as.numeric(GTOP_GTEx_freq$af_EUR)
GTOP_GTEx_freq$af_GTEx_raw[GTOP_GTEx_freq$af_GTEx_raw == 0] <- NA

GTOP_GTEx_freq$af_GTEx <- case_when(
  is.na(GTOP_GTEx_freq$af_GTEx_raw) ~ GTOP_GTEx_freq$af_EUR,
  TRUE ~ GTOP_GTEx_freq$af_GTEx_raw
)

table(is.na(GTOP_GTEx_freq$af_GTEx_raw))
table(is.na(GTOP_GTEx_freq$af_GTEx))

GTOP_GTEx_freq <- GTOP_GTEx_freq[!is.na(GTOP_GTEx_freq$af_GTEx), ]

table(abs(na.omit(GTOP_GTEx_freq$af_GTEx_raw-GTOP_GTEx_freq$af_EUR))<0.2)

# GTOP_GTEx_freq <- merge(SNV_freq_1KG_res, GTEx_freq[, .(chr_pos_ref_alt, af_GTEx=ALT_FREQS)], 
#                            by="chr_pos_ref_alt", all.x = TRUE)
# GTOP_GTEx_freq$af_GTEx[is.na(GTOP_GTEx_freq$af_GTEx)] <- 0
# table(GTOP_GTEx_freq$af_GTEx==0 & GTOP_GTEx_freq$af_EUR != 0)
# GTOP_GTEx_freq <- GTOP_GTEx_freq[abs(GTOP_GTEx_freq$af_GTEx-GTOP_GTEx_freq$af_EUR)<0.2, ]

# GTOP_GTEx_freq$af_GTEx[GTOP_GTEx_freq$af_GTEx==0] <- GTOP_GTEx_freq$af_EUR[GTOP_GTEx_freq$af_GTEx==0]
# table(GTOP_GTEx_freq$af_GTEx==0)
nrow(GTOP_GTEx_freq)

GTOP_GTEx_freq$Type_GTEx <- "C"
index1 <- GTOP_GTEx_freq$af_GTEx < 0.01 | GTOP_GTEx_freq$af_GTEx > 0.99
index2 <- GTOP_GTEx_freq$af_GTEx == 0 | GTOP_GTEx_freq$af_GTEx == 1
GTOP_GTEx_freq$Type_GTEx[index1] <- "R"
GTOP_GTEx_freq$Type_GTEx[index2] <- "U"
table(GTOP_GTEx_freq$Type_GTEx)

## 6035957 6664773 6647083

fwrite(GTOP_GTEx_freq, "2025-06-11-specific_xQTL/output/data/processing/GTOP_1KG_GTEx_SNV_freq.txt.gz", sep = "\t")

# GTOP_GTEx_freq <- fread("2025-06-11-specific_xQTL/output/data/processing/GTOP_1KG_GTEx_SNV_freq.txt.gz")

## GTEx allele frequency in fina-mapping result
onehot_reshape <- function(idata, icol_name) {
  # 
  col_sym <- enquo(icol_name)
  col_name_str <- rlang::as_name(col_sym)
  new_col_name <- paste0("new_", col_name_str)
  # 1. one-hot
  onehot_data <- idata %>%
    dplyr::select(!!col_sym) %>%
    model.matrix(~ . - 1, data = .)
  
  oh_cols <- colnames(onehot_data)
  if (length(oh_cols) == 0) {
    stop("No one-hot encoded columns generated.")
  }
  
  # 2. weight
  pip_weights <- idata $ pip
  onehot_weighted <- sweep(onehot_data, 1, pip_weights, `*`)
  
  # 3. merge meta data
  combined <- idata %>%
    dplyr::select(variant_proxy, !!col_sym) %>%
    as.data.frame() %>%
    cbind(onehot_weighted)
  
  # 4. variant_proxy
  result <- combined %>%
    group_by(variant_proxy) %>%
    summarise(
      across(all_of(oh_cols), sum, .names = "{.col}"),
      .groups = "drop"
    )
  
  # 5. max weighted value
  value_matrix <- as.matrix(result[oh_cols])
  max_idx <- max.col(value_matrix, ties.method = "first")
  max_names <- oh_cols[max_idx]
  
  # 6. level name
  new_labels <- sub(paste0("^", col_name_str), "", max_names)
  
  # 7. add information
  result <- result %>%
    mutate(!!new_col_name := new_labels)
  
  result_final <- idata %>%
    # dplyr::select(-all_of(c(col_name_str))) %>% 
    inner_join(result %>% dplyr::select(variant_proxy, !!new_col_name), 
               by = "variant_proxy")
  
  return(result_final)
}


eQTL_fm_res <- fread("2025-06-11-specific_xQTL/output/data/eQTL/finemapping_merge_variants.txt")
eQTL_fm_freq <- merge(eQTL_fm_res, GTOP_GTEx_freq[, .(variant_id, chr_pos_ref_alt, af_GTOP, 
                                                          af_GTEx, Type_GTEx)], 
                      by="variant_id")
eQTL_fm_freq <- onehot_reshape(idata = eQTL_fm_freq, icol_name = "Type_GTEx")
length(unique(eQTL_fm_freq$variant_proxy))

fwrite(eQTL_fm_freq, "2025-06-11-specific_xQTL/output/data/eQTL/finemapping_add_GTEx_freq_all.txt.gz", sep = "\t")


sQTL_fm_res <- fread("2025-06-11-specific_xQTL/output/data/sQTL/finemapping_merge_variants.txt")
sQTL_fm_freq <- merge(sQTL_fm_res, GTOP_GTEx_freq[, .(variant_id, chr_pos_ref_alt, af_GTOP, 
                                                         af_GTEx, Type_GTEx)], 
                      by="variant_id")
sQTL_fm_freq <- onehot_reshape(idata = sQTL_fm_freq, icol_name = "Type_GTEx")
length(unique(sQTL_fm_freq$variant_proxy))

fwrite(sQTL_fm_freq, "2025-06-11-specific_xQTL/output/data/sQTL/finemapping_add_GTEx_freq_all.txt.gz", sep = "\t")
