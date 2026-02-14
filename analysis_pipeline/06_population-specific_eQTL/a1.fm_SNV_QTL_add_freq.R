#!/usr/bin/env Rscript

setwd("/media/bora_A/zhangt/2025-05-07-EAS_specific_xQTL-Project")

# Load required libraries
library(data.table)
library(tidyverse)

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


clean_anno <- fread("2025-06-11-specific_xQTL/input/GTOP/annotation/GTOP_gencode_v47.gene_info.txt") %>% as.data.frame()
rownames(clean_anno) <- clean_anno$V5

CHRPOS_RSID <- fread("2025-06-11-specific_xQTL/input/GTOP/SNV/GTOP_chrpos_rsid.txt", 
                     header = FALSE, col.names = c("chr_pos_ref_alt", "variant_id"))


SNV_freq_1KG_res <- fread("2025-06-11-specific_xQTL/output/data/processing/GTOP_1KG_SNV_freq.txt.gz")
dim(SNV_freq_1KG_res)
table(SNV_freq_1KG_res$Type)
table(abs(SNV_freq_1KG_res$af_GTOP-SNV_freq_1KG_res$af_EAS)<0.2)

SNV_freq_1KG_res$Type_AFR <- case_when(
  SNV_freq_1KG_res$af_AFR == 0 | SNV_freq_1KG_res$af_AFR == 1 ~ "U",
  SNV_freq_1KG_res$af_AFR < 0.01 | SNV_freq_1KG_res$af_AFR > 0.99 ~ "R",
  TRUE ~ "C"
)
SNV_freq_1KG_res$Type_EUR <- case_when(
  SNV_freq_1KG_res$af_EUR == 0 | SNV_freq_1KG_res$af_EUR == 1 ~ "U",
  SNV_freq_1KG_res$af_EUR < 0.01 | SNV_freq_1KG_res$af_EUR > 0.99 ~ "R",
  TRUE ~ "C"
)
SNV_freq_1KG_res$Type_SAS <- case_when(
  SNV_freq_1KG_res$af_SAS == 0 | SNV_freq_1KG_res$af_SAS == 1 ~ "U",
  SNV_freq_1KG_res$af_SAS < 0.01 | SNV_freq_1KG_res$af_SAS > 0.99 ~ "R",
  TRUE ~ "C"
)
SNV_freq_1KG_res$Type_AMR <- case_when(
  SNV_freq_1KG_res$af_AMR == 0 | SNV_freq_1KG_res$af_AMR == 1 ~ "U",
  SNV_freq_1KG_res$af_AMR < 0.01 | SNV_freq_1KG_res$af_AMR > 0.99 ~ "R",
  TRUE ~ "C"
)

## finemapping result -----------------------------------------------
susie_res <- fread("2025-06-11-specific_xQTL/output/data/SNV_xQTL_finemapping.txt")

length(unique(paste0(susie_res$tissue, "_", susie_res$gene_name, "_", susie_res$cs_id)[susie_res$xQTL_type=="eQTL"]))
length(unique(paste0(susie_res$tissue, "_", susie_res$gene_name, "_", susie_res$cs_id)[susie_res$xQTL_type!="eQTL"]))

length(unique(susie_res$gene_name[susie_res$xQTL_type=="eQTL"]))
length(unique(susie_res$gene_name[susie_res$xQTL_type!="eQTL"]))

susie_res$xQTL_type <- case_when(
  susie_res$xQTL_type %in% c("juQTL", "tuQTL") ~ "sQTL",
  TRUE ~ "eQTL"
)


for(qtl_type in c("eQTL","sQTL")){
  
  output_dir <- sprintf("2025-06-11-specific_xQTL/output/data/%s", qtl_type)
  if(!file.exists(output_dir)){dir.create(output_dir, recursive = T)}
  
  ## eQTL/sQTL
  susie_subres <- susie_res[xQTL_type==qtl_type]
  susie_subres <- susie_subres %>% 
    arrange(tissue, phenotype_id, cs_size, desc(pip), variant_id)
  table(unique(susie_subres$variant_id)%in%SNV_freq_1KG_res$variant_id)
  
  ## selected variant with highest pip value
  susie_subres <- susie_subres  %>% 
    group_by(phenotype_id, gene_name, gene_symbol, gene_type, 
             cs_id, cs_size, 
             tissue, 
             xQTL_type) %>%
    mutate(variant_proxy=variant_id[pip==max(pip)][1])
  setDT(susie_subres)
  
  fm_index <- paste0(susie_subres$tissue, "_", 
                              susie_subres$cs_id, "_",
                              susie_subres$gene_name)

  lead_fm_info <- unique(susie_subres[variant_proxy==variant_id, .(phenotype_id, gene_name, cs_id, tissue, variant_proxy, pip)])
  
  ## split into genes
  lead_each_gene <- split(lead_fm_info, lead_fm_info$gene_name)
  
  merge_cs_res <- pbmcapply::pbmclapply(1:length(lead_each_gene), function(row_index){
    ## 
    top_gene <- names(lead_each_gene)[row_index]
    current_subset <- lead_each_gene[[row_index]]
    gene_variant <- current_subset[gene_name == top_gene]
    
    variant_list <- list()
    remain_variants <- unique(gene_variant$variant_proxy)
    
    while(length(remain_variants) > 0) {
      current_qtl <- gene_variant[variant_proxy %in% remain_variants]
      top_variant <- current_qtl$variant_proxy[which.max(current_qtl$pip)][1]
      
      top_variant_group <- susie_subres[gene_name==top_gene & variant_id==top_variant]
      
      flanking_index <- paste0(top_variant_group$tissue, "_", 
                               top_variant_group$cs_id, "_",
                               top_variant_group$gene_name)
      
      flanking_variant1 <- susie_subres[fm_index%in%flanking_index]
      high_ld_snps <- setdiff(intersect(remain_variants, flanking_variant1$variant_id), top_variant)
      
      if(length(high_ld_snps) >= 1){
        variant_list[[top_variant]] <- paste0(c(top_variant, high_ld_snps), collapse = ",")
        remain_variants <- setdiff(remain_variants, c(high_ld_snps, top_variant))
      }else{
        variant_list[[top_variant]] <- top_variant
        remain_variants <- setdiff(remain_variants, c(top_variant))
      }
    }
    variant_list_info <- lapply(variant_list, function(x){strsplit(x, ",")[[1]]})
    variant_list_df <- data.table("gene_name" = top_gene,
                                  "variant_proxy" = unlist(variant_list_info),
                                  "variant_proxy_c" = rep(names(variant_list_info), sapply(variant_list_info, length)))
    variant_list_df
  }, mc.cores = 40, mc.preschedule = FALSE)
  
  merge_cs_df <- do.call('rbind', merge_cs_res) %>% as.data.table()
  
  fwrite(merge_cs_df, file = sprintf("2025-06-11-specific_xQTL/output/data/%s/finemapping_variant_change.txt", 
                                     qtl_type), sep = "\t")
  
  susie_subres_c <- merge(susie_subres, merge_cs_df, by=c("variant_proxy", "gene_name"))
  
  length(unique(susie_subres_c$variant_proxy))
  length(unique(susie_subres_c$variant_proxy_c))
  

  fm_res <- susie_subres_c[, .(variant_proxy=variant_proxy_c, variant_id, 
                             phenotype_id, pip, cs_id, cs_size, tissue,
                             xQTL_type, gene_name, gene_symbol, gene_type)]
  
  fwrite(fm_res, file = sprintf("2025-06-11-specific_xQTL/output/data/%s/finemapping_merge_variants.txt", 
                                               qtl_type), sep = "\t")
  
  ## SNV -------
  # fm_res <- fread(sprintf("2025-06-11-specific_xQTL/output/data/%s/finemapping_merge_variants.txt", qtl_type))
  
  fm_freq_res <- merge(fm_res, SNV_freq_1KG_res[, .(variant_id, Type_AFR, Type_EUR, Type_SAS, Type_AMR, Type)], 
                  by="variant_id")
  
  fm_freq_res <- onehot_reshape(idata = fm_freq_res, icol_name = "Type_AFR")
  fm_freq_res <- onehot_reshape(idata = fm_freq_res, icol_name = "Type_EUR")
  fm_freq_res <- onehot_reshape(idata = fm_freq_res, icol_name = "Type_SAS")
  fm_freq_res <- onehot_reshape(idata = fm_freq_res, icol_name = "Type_AMR")
  
  fm_freq_res$new_Type <- case_when(
    fm_freq_res$new_Type_AFR %in% c("R", "U") &
      fm_freq_res$new_Type_EUR %in% c("R", "U") &
      fm_freq_res$new_Type_SAS %in% c("R", "U") &
      fm_freq_res$new_Type_AMR %in% c("R", "U") ~ "R2",
    fm_freq_res$new_Type_AFR %in% c("R", "U") |
      fm_freq_res$new_Type_EUR %in% c("R", "U") |
      fm_freq_res$new_Type_SAS %in% c("R", "U") |
      fm_freq_res$new_Type_AMR %in% c("R", "U") ~ "R1",
    TRUE ~ "C"
  )
  

  fwrite(fm_freq_res, file = sprintf("2025-06-11-specific_xQTL/output/data/%s/finemapping_SNV_add_four_freq.txt", qtl_type), sep = "\t")
}


## fm_res
eQTL_fm_res <- fread(sprintf("2025-06-11-specific_xQTL/output/data/%s/finemapping_merge_variants.txt", "eQTL"))

eQTL_fm_proxy <- eQTL_fm_res %>% group_by(gene_name) %>% summarise(cs_num=length(unique(variant_proxy)))
eQTL_fm_proxy$cs_num[eQTL_fm_proxy$cs_num>=7] <- 7
eQTL_fm_proxy <- eQTL_fm_proxy %>% group_by(cs_num) %>% summarise(phenotype_count=length(unique(gene_name)))

sQTL_fm_res <- fread(sprintf("2025-06-11-specific_xQTL/output/data/%s/finemapping_merge_variants.txt", "sQTL"))

sQTL_fm_proxy <- sQTL_fm_res %>% group_by(gene_name) %>% summarise(cs_num=length(unique(variant_proxy)))
sQTL_fm_proxy$cs_num[sQTL_fm_proxy$cs_num>=7] <- 7
sQTL_fm_proxy <- sQTL_fm_proxy %>% group_by(cs_num) %>% summarise(phenotype_count=length(unique(gene_name)))

fm_proxy_count <- rbind(eQTL_fm_proxy %>% mutate("xQTL_type"="eQTL"),
                        sQTL_fm_proxy %>% mutate("xQTL_type"="sQTL"))

ggplot(fm_proxy_count, aes(x=cs_num, y=phenotype_count, fill=xQTL_type)) +
  geom_col(position = "dodge") +
  scale_x_continuous(breaks = 1:7, labels = c(1:6,">=7")) +
  theme_classic() +
  labs(x="Number of independent QTLs", y="Gene number") +
  scale_fill_manual(values = c("eQTL"="#9ac294", "sQTL"="#7b86a7"))

ggsave("2025-06-11-specific_xQTL/output/result/Gene_number_of_independent_QTL.pdf", width = 4, height = 3)

fwrite(fm_proxy_count, "/media/london_A/mengxin/GTOP_code/extend/extend_6/input/Extended_fig6a_data.txt", sep = "\t")


## eQTL for EAS-specific SNVs
SNV_eQTL_fm_freq <- fread("2025-06-11-specific_xQTL/output/data/eQTL/finemapping_SNV_add_four_freq.txt") 
length(unique(SNV_eQTL_fm_freq$variant_proxy)) # 11,121

SNV_eQTL_fm_freq_lead <- unique(SNV_eQTL_fm_freq[Type_AFR==new_Type_AFR & Type_EUR==new_Type_EUR &
                                                   Type_SAS==new_Type_SAS & Type_AMR==new_Type_AMR])
length(unique(SNV_eQTL_fm_freq_lead$variant_proxy)) # 11,105

SNV_eQTL_fm_freq_lead <- SNV_eQTL_fm_freq_lead %>%
  group_by(variant_proxy, gene_name) %>%
  mutate(
    is_same = variant_id == variant_proxy
  ) %>%
  dplyr::slice(which.max(ifelse(is_same, 1, -pip))) %>%
  dplyr::select(-is_same)
setDT(SNV_eQTL_fm_freq_lead)

SNV_eQTL_fm_freq_lead <- unique(SNV_eQTL_fm_freq_lead[, 
                          .(SNP=variant_id, AFR=new_Type_AFR, EUR=new_Type_EUR, 
                            SAS=new_Type_SAS, AMR=new_Type_AMR)])
length(unique(SNV_eQTL_fm_freq_lead$SNP))
nrow(SNV_eQTL_fm_freq_lead) # 11,057

SNV_eQTL_fm_freq_lead[, (2:5) := lapply(.SD, function(x) ifelse(x == "U", 0, ifelse(x == "R", 0.01, ifelse(x == "C", 0.1, x)))), .SDcols = 2:5]
