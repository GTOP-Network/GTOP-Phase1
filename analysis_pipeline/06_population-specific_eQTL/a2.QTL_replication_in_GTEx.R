
library(data.table)
library(tidyverse)

setwd("/media/bora_A/zhangt/2025-05-07-EAS_specific_xQTL-Project/")

## compare with GTEx
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

GTEx_support_res <- pbmcapply::pbmclapply(seq_len(length(tissue_change)), function(f_index){
  GTEx_V8_egenes <- fread(sprintf("2025-06-11-specific_xQTL/data/tensorqtl/GTExv8_eQTL/permutation/%s.txt.gz", 
                                  tissue_change[f_index]))
  GTEx_V8_egenes <- GTEx_V8_egenes[qval < 0.05]
  
  GTEx_V8_sig <- fread(sprintf("2025-06-11-specific_xQTL/data/tensorqtl/GTExv8_eQTL/nominal_slim/nom_thresh/%s.nom_thresh.txt.gz", 
                               tissue_change[f_index]))
  GTEx_V8_sig <- GTEx_V8_sig[, c("phenotype_id", "chr_pos_ref_alt")]
  
  GTEx_V8_sig <- inner_join(GTEx_V8_sig, GTEx_V8_egenes[, .(phenotype_id=gene_id, gene_symbol=gene_name)], by="phenotype_id")
  GTEx_V8_sig$tissue <- names(tissue_change)[f_index]
  return(GTEx_V8_sig)
}, mc.cores = 10, mc.preschedule = FALSE)
GTEx_support_df <- rbindlist(GTEx_support_res)
GTEx_support_df$gene_tmpid <- gsub("\\..+", "", GTEx_support_df$phenotype_id)
GTEx_support_df$chr_pos_ref_alt <- gsub("_b38", "", GTEx_support_df$chr_pos_ref_alt)


GTEx_expGene_res <- pbmcapply::pbmclapply(seq_len(length(tissue_change)), function(f_index){
  GTEx_V8_egenes <- fread(sprintf("2025-06-11-specific_xQTL/data/tensorqtl/GTExv8_eQTL/permutation/%s.txt.gz", 
                                  tissue_change[f_index]))
  GTEx_V8_egenes <- GTEx_V8_egenes
  
  GTEx_V8_egenes$tissue <- names(tissue_change)[f_index]
  return(GTEx_V8_egenes)
}, mc.cores = 10, mc.preschedule = FALSE)
GTEx_expGene_df <- rbindlist(GTEx_expGene_res)
GTEx_expGene_df$gene_tmpid <- gsub("\\..+", "", GTEx_expGene_df$gene_id)
GTEx_expGene_id <- paste0(GTEx_expGene_df$tissue, "_", GTEx_expGene_df$gene_tmpid)


GTEx_eGene_res <- pbmcapply::pbmclapply(seq_len(length(tissue_change)), function(f_index){
  GTEx_V8_egenes <- fread(sprintf("2025-06-11-specific_xQTL/data/tensorqtl/GTExv8_eQTL/permutation/%s.txt.gz", 
                                  tissue_change[f_index]))
  GTEx_V8_egenes <- GTEx_V8_egenes[qval < 0.05]
  
  GTEx_V8_egenes$tissue <- names(tissue_change)[f_index]
  return(GTEx_V8_egenes)
}, mc.cores = 10, mc.preschedule = FALSE)
GTEx_eGene_df <- rbindlist(GTEx_eGene_res)
GTEx_eGene_df$gene_tmpid <- gsub("\\..+", "", GTEx_eGene_df$gene_id)


save(GTEx_expGene_df, GTEx_eGene_df, GTEx_support_df, file = "2025-06-11-specific_xQTL/output/data/QC/GTEx_exp_eGene_info.RData")


## Gene information
GTEx_annotation <- fread("2025-06-11-specific_xQTL/data/tensorqtl/GTExv8_eQTL/gencode.v26.gene.loci")
GTOP_annotation <- fread("2025-06-11-specific_xQTL/input/GTOP/annotation/GTOP_gencode_v47.gene_info.txt")

GTEx_annotation$gene_tmpid <- gsub("\\..+", "", GTEx_annotation$V1)
GTOP_annotation$gene_tmpid <- gsub("\\..+", "", GTOP_annotation$V5)

ggvenn::ggvenn(list("GTEx"=GTEx_annotation$gene_tmpid, "GTOP"=GTOP_annotation$gene_tmpid))

overlapped_genes <- intersect(GTOP_annotation$gene_tmpid, GTEx_annotation$gene_tmpid)


GTOP_eGene_files <- list.files("/media/bora_A/zhangt/2025-05-07-EAS_specific_xQTL-Project/2025-06-11-specific_xQTL/data/tensorqtl/SNV_eQTL/clean_xGene")
GTOP_eGene_info <- rbindlist(lapply(GTOP_eGene_files, function(tmp_file){
   f_info <- fread(paste0("/media/bora_A/zhangt/2025-05-07-EAS_specific_xQTL-Project/2025-06-11-specific_xQTL/data/tensorqtl/SNV_eQTL/clean_xGene/", tmp_file))
   f_info[, .(chr_pos_ref_alt, phenotype_id, tissue=gsub(".txt", "", tmp_file))]
}))

GTOP_eGene_info$gene_tmpid <- gsub("\\..+", "", GTOP_eGene_info$phenotype_id)
table(GTOP_eGene_info$gene_tmpid %in% overlapped_genes )

GTOP_eGene_info <- GTOP_eGene_info[gene_tmpid%in%overlapped_genes]
# fwrite(GTOP_eGene_info, file = "2025-06-11-specific_xQTL/output/data_plot/data/GTOP_eGene_info.txt", sep = "\t")

GTOP_sigQTL_files <- list.files("/media/bora_A/zhangt/2025-05-07-EAS_specific_xQTL-Project/2025-06-11-specific_xQTL/data/tensorqtl/SNV_eQTL/nominal_slim/nom_thresh")
GTOP_sigQTL_info <- rbindlist(lapply(GTOP_sigQTL_files, function(tmp_file){
  f_info <- fread(paste0("/media/bora_A/zhangt/2025-05-07-EAS_specific_xQTL-Project/2025-06-11-specific_xQTL/data/tensorqtl/SNV_eQTL/nominal_slim/nom_thresh/", tmp_file))
  f_info[, .(chr_pos_ref_alt, phenotype_id, tissue=gsub(".nom_thresh.txt.gz", "", tmp_file))]
}))
GTOP_sigQTL_info$gene_tmpid <- gsub("\\..+", "", GTOP_sigQTL_info$phenotype_id)

save(GTOP_eGene_info, GTOP_sigQTL_info, file = "2025-06-11-specific_xQTL/output/data/QC/GTOP_eGene_info.RData")


## Overlap
SNV_sigQTL_subinfo <- merge(GTOP_sigQTL_info, GTEx_support_df, by=c("gene_tmpid", "chr_pos_ref_alt", "tissue"))

GTOP_eGene_info$GTEx_support <- "-"
tmp_id <- paste0(GTOP_eGene_info$tissue, "_", GTOP_eGene_info$gene_tmpid)
table(tmp_id%in%paste0(SNV_sigQTL_subinfo$tissue, "_", SNV_sigQTL_subinfo$gene_tmpid))
GTOP_eGene_info$GTEx_support[tmp_id%in%paste0(SNV_sigQTL_subinfo$tissue, "_", SNV_sigQTL_subinfo$gene_tmpid)] <- "Yes"

GTOP_eGene_info$GTEx_expGene <- "-"
tmp_id <- paste0(GTOP_eGene_info$tissue, "_", GTOP_eGene_info$gene_tmpid)
GTEx_expGene_id <- paste0(GTEx_expGene_df$tissue, "_", GTEx_expGene_df$gene_tmpid)
table(tmp_id%in%GTEx_expGene_id)
GTOP_eGene_info$GTEx_expGene[tmp_id%in%GTEx_expGene_id] <- "Yes"

GTOP_eGene_info$GTEx_eGene <- "-"
tmp_id <- paste0(GTOP_eGene_info$tissue, "_", GTOP_eGene_info$gene_tmpid)
GTEx_eGene_id <- paste0(GTEx_eGene_df$tissue, "_", GTEx_eGene_df$gene_tmpid)
table(tmp_id%in%GTEx_eGene_id)
GTOP_eGene_info$GTEx_eGene[tmp_id%in%GTEx_eGene_id] <- "Yes"

GTE_comparison <- GTOP_eGene_info %>% 
  group_by(tissue, GTEx_expGene, GTEx_eGene, GTEx_support) %>% 
  summarise(count=length(unique(phenotype_id))) %>% 
  as.data.table()

GTE_comparison$Type <- paste0(GTE_comparison$GTEx_expGene,
                              "_", GTE_comparison$GTEx_eGene, "_", GTE_comparison$GTEx_support)
table(GTE_comparison$Type)
GTE_comparison$Type <- factor(GTE_comparison$Type, levels = c(
  "-_-_-", "Yes_-_-", "Yes_Yes_-", "Yes_Yes_Yes"
))

GTE_comparison <- GTE_comparison %>% group_by(tissue) %>% mutate(allcount=sum(count))
GTE_comparison <- GTE_comparison %>% arrange(desc(allcount))
GTE_comparison <- GTE_comparison[GTE_comparison$tissue!="Gallbladder", ]
GTE_comparison$ntissue <- as.character(paste0(GTE_comparison$tissue, "\n(n = ", GTE_comparison$allcount, ")"))
GTE_comparison$ntissue <- factor(GTE_comparison$ntissue, levels = rev(unique(GTE_comparison$ntissue)))

GTE_comparison$ratio <- GTE_comparison$count/GTE_comparison$allcount

mean(GTE_comparison$ratio[GTE_comparison$Type=="-_-_-"])
mean(GTE_comparison$ratio[GTE_comparison$Type=="Yes_Yes_-"])
mean(GTE_comparison$ratio[GTE_comparison$Type=="Yes_-_-"] + GTE_comparison$ratio[GTE_comparison$Type=="Yes_Yes_-"])

color_df <- fread("/media/pacific/share/Datasets/Asian_GTEx/Metainfo/input/GMTiP_tissue_code_and_colors.csv")
color_vec <- paste0("#", color_df$Tissue_Color_Code)
names(color_vec) <- color_df$Tissue

tissue_df <- unique(GTE_comparison[, c("tissue", "ntissue")]) %>% as.data.frame()
rownames(tissue_df) <- tissue_df$ntissue

ggplot(GTE_comparison, aes(x = count, y=ntissue, fill = Type)) +
  geom_col(position = "fill") +
  # facet_wrap(tissue~., nrow = 1) +
  labs(x = "Proportion of eGenes", y = "") +
  theme_classic() +
  theme(
    axis.line = element_line(color = "black"),
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(color = "black", size = 12)
  ) +
  scale_fill_manual(values = c("#9d3929", "#fc9272", "#ab889a", "#3578ac","1")) +
  scale_color_manual(values = color_vec)
