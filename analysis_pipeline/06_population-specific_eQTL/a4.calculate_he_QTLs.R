
library(data.table)
library(tidyverse)

mashr_res_dir <- "2025-09-30-mash/output"
each_xQTLtype <- "SNV_eQTL"

## beta se values
strong_beta_se <- vroom::vroom(sprintf("%s/%s/Strong2/strong_beta_se.MashR_input.txt.gz", 
                                       mashr_res_dir, each_xQTLtype)) %>% 
  as.data.frame() %>% 
  column_to_rownames("pair_id")

raw_z.stat <- strong_beta_se[, seq(3, ncol(strong_beta_se),3)]
colnames(raw_z.stat) <- gsub("_zval", "", colnames(raw_z.stat))

raw_beta <- strong_beta_se[, seq(1, ncol(strong_beta_se),3)]
colnames(raw_beta) <- gsub("_slope", "", colnames(raw_beta))

## lfsr
posterior_list <- readRDS(sprintf("%s/%s/top_pairs2/m.s_zval.RDS", 
                                  mashr_res_dir, each_xQTLtype))
lfsr_mash <- posterior_list$result$lfsr

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
new_tissue_change <- paste0("GTEx_", tissue_change)
names(new_tissue_change) <- paste0("GTOP_", names(tissue_change))
z_score_cutoff <- -qnorm(1e-3/2)
-qnorm(1e-5/2)

eQTL_fm_freq <- fread("2025-06-11-specific_xQTL/output/data/eQTL/finemapping_add_GTEx_freq_all.txt.gz")

sig_res_list <- lapply(names(new_tissue_change), function(x){
  ## tissue information
  col1 <- x
  col2 <- new_tissue_change[x]
  
  ##
  tissue_fm <- eQTL_fm_freq[eQTL_fm_freq$tissue==gsub("GTOP_","",x),]
  tissue_fm$phenotype_id <- case_when(
    grepl("ENSG", tissue_fm$phenotype_id) ~ gsub("\\.[0-9]+", "", tissue_fm$phenotype_id),
    TRUE ~ tissue_fm$phenotype_id
  )
  
  ##
  tmp_lfsr <- as.data.frame(lfsr_mash[, c(col1, col2)])
  
  tmp_zscore <- raw_z.stat[rownames(tmp_lfsr), colnames(tmp_lfsr)]
  colnames(tmp_zscore) <- paste0("Z_", colnames(tmp_zscore))
  
  tmp_beta <- raw_beta[rownames(tmp_lfsr), colnames(tmp_lfsr)]
  colnames(tmp_beta) <- paste0("beta_", colnames(tmp_beta))
  
  tmp_lfsr <- cbind(tmp_lfsr, tmp_zscore, tmp_beta)
  tmp_lfsr$chr_pos_ref_alt <- sapply(strsplit(rownames(tmp_lfsr), split = ","), function(x){x[1]})
  tmp_lfsr$phenotype_id <- sapply(strsplit(rownames(tmp_lfsr), split = ","), function(x){x[2]})
  tmp_lfsr$ID <- rownames(tmp_lfsr)
  
  tmp_lfsr1 <- inner_join(tmp_lfsr, as.data.frame(tissue_fm[,.(chr_pos_ref_alt, variant_proxy, phenotype_id)]),
                    by=c("chr_pos_ref_alt", "phenotype_id"))
    
  # tmp_lfsr3 <- tmp_lfsr1 %>%
  #   group_by(phenotype_id, variant_proxy) %>%
  #   summarise(
  #     TF1 = all(.data[[col1]] < 0.05, na.rm = TRUE),
  #     TF2 = all(.data[[col2]] > 0.05, na.rm = TRUE),
  #     # TF3 = all(!is.na(.data[[paste0("Z_", col1)]])),
  #       # all(.data[[paste0("Z_", col1)]] > z_score_cutoff, na.rm = TRUE),
  #     # TF4 = all(!is.na(.data[[paste0("Z_", col2)]]), na.rm = TRUE),
  #     TF3 = nrow(.data)==1 & all(!is.na(.data[[paste0("Z_", col2)]])),
  #     TF4 = nrow(.data)>1 & sum(!is.na(.data[[paste0("Z_", col2)]]))>1,
  #     .groups = "drop"
  #   )
  # 
  tmp_lfsr2 <- tmp_lfsr1 %>%
    group_by(phenotype_id, variant_proxy) %>%
    summarise(
      TF1 = all(.data[[col1]] < 0.05, na.rm = TRUE),
      TF2 = all(.data[[col2]] > 0.05, na.rm = TRUE),
      TF3 = (n() <= 1 & all(!is.na(.data[[paste0("Z_", col2)]][1]))) | (n() >= 2 & sum(!is.na(.data[[paste0("Z_", col2)]])) >= 2),
      .groups = "drop"
    )
  
  tmp_lfsr2 <- tmp_lfsr2[rowSums(tmp_lfsr2[, 3:5])==3, ]
  # tmp_lfsr2 <- tmp_lfsr2[rowSums(tmp_lfsr2[, 3:6])==4, ]
  tmp_lfsr3 <- inner_join(tmp_lfsr1, tmp_lfsr2[, c("phenotype_id", "variant_proxy")], 
                          by=c("phenotype_id", "variant_proxy"))
  
  return(data.frame("tissue"=x, "xQTL"=tmp_lfsr3$ID,
                    "gene"=tmp_lfsr3$phenotype_id,
                    "variant_id"=tmp_lfsr3$chr_pos_ref_alt,
                    "variant_proxy"=tmp_lfsr3$variant_proxy,
                    "GTOP_lfsr"=tmp_lfsr3[,1], 
                    "GTEx_lfsr"=tmp_lfsr3[,2],
                    "GTOP_Z"=tmp_lfsr3[,3],
                    "GTEx_Z"=tmp_lfsr3[,4],
                    "GTOP_beta"=tmp_lfsr3[,5],
                    "GTEx_beta"=tmp_lfsr3[,6]))
})


sig_res_df <- do.call('rbind', sig_res_list)
sig_res_df %>% group_by(tissue) %>% summarise(count=length(unique(gene)))

gene_annotation  <- fread("2025-06-11-specific_xQTL/input/GTOP/annotation/GTOP_gencode_v47.gene_info.txt") %>% as.data.frame()
rownames(gene_annotation) <- case_when(
  grepl("^ENSG", gene_annotation$V5) ~ gsub("\\.[0-9]+", "", gene_annotation$V5),
  TRUE ~ gene_annotation$V5
)

sig_res_df$gene_symbol <- gene_annotation[sig_res_df$gene, "V6"]
sig_res_df$gene_type <- gene_annotation[sig_res_df$gene, "V7"]
sig_res_df$tissue <- gsub("GTOP_", "", sig_res_df$tissue)

fwrite(sig_res_df, "2025-06-11-specific_xQTL/output/data/he_QTL_res.txt", sep = "\t")


GTOP_coloc <- fread("2025-09-28-coloc/output/result/GTOP_EASGWAS_coloc_result.txt")
GTOP_SMR <- fread("2025-09-28-coloc/output/result/GTOP_EASGWAS_SMR_result.txt")

sig_res_subdf <- merge(sig_res_df, GTOP_SMR[, .(Disease_trait, topSNP, 
                                                  tissue=paste0("GTOP_", tissue), gene_symbol)], 
                       by=c("tissue", "gene_symbol"))
