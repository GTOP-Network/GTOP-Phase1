#==================================#
# long & short RNA #
# Figure-2 #
#==================================#
library(data.table)
library(ggplot2)
setwd("/path/to/GTOP_code/fig-2/input/")
# Fig.2a LRS_RNA_MDS ------------------------------------------------------



THEME_SIZE <- 12
TITLE_SIZE <- 14
POINT_SIZE <- 1.8
LINE_WIDTH <- 0.6
ALPHA_NORMAL <- 0.7
ALPHA_OUTLIER <- 0.9

unified_theme <- theme_bw(base_size = THEME_SIZE) +
  theme(
    plot.title = element_text(hjust = 0.5, size = TITLE_SIZE, face = "bold"),
    axis.title = element_text(size = THEME_SIZE, face = "bold"),
    axis.text = element_text(size = THEME_SIZE - 1),
    legend.title = element_text(size = THEME_SIZE, face = "bold"),
    legend.text = element_text(size = THEME_SIZE - 1),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  )

load("sample_annot_full_1586.RData")

final_annot <- sample_annot_full[, c("sample_id", "Subject", "Tissue","Batch", "Tissue_Color_Code")]
colnames(final_annot) <- c("sample", "individual", "tissue", "batch", "color")
rownames(final_annot) <- final_annot$sample

meta <- final_annot[, c("sample", "individual", "tissue", "batch", "color")]
rownames(meta) <- meta$sample

# Color setup -------------------------------------------------------------

all_tissue <- sort(unique(meta$tissue))
tissue_colors <- unique(meta[, c("tissue", "color")])
tissue_colors <- setNames(paste0("#", tissue_colors$color), tissue_colors$tissue)

mds_df_corrected<-fread("LR_Sample_MDS_coor.txt")
mds_coor<-fread("LR_Sample_MDS_var.txt")

prop1_corr<-mds_coor[1,1]
prop2_corr<-mds_coor[1,2]
# MDS plot  ---------------------------------------------------------------

gg_mds_corr <- ggplot(mds_df_corrected, aes(x = MDS1, y = MDS2, color = Tissue)) +
  geom_point(size = 2, alpha = 1) +
  scale_color_manual(values = tissue_colors) +
  scale_y_continuous(
    breaks = seq(-75, 75, by = 25)
  ) +
  unified_theme +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = paste0("Coordinate 1 (", prop1_corr, ")"),
    y = paste0("Coordinate 2 (", prop2_corr, ")"),
    color = "Tissue"
  )
print(gg_mds_corr)


# Fig.2b LR_transcript_novel_stat_Transcript --------------------------------------------
library(data.table)
library(ggplot2)

dat<-fread("fig2b.LR_transcript_novel_stat_Transcript.txt")
dat_long <- melt(
  dat,
  id.vars = "index",
  variable.name = "Type",
  value.name = "Count"
)
dat_long$Type <- factor(
  dat_long$Type,
  levels = c( "Novel","Annotated", "GENCODE v47")
)
dat_long$index <- factor(
  dat_long$index,
  levels = c("GTOP", "GENCODE", "GTOP + GENCODE")
)
type_colors <- c(
  "GENCODE v47"        = "#9faac1",
  "Annotated"  = "#c3968f",
  "Novel"          = "#9d3929"
)
dat_long$Count<-dat_long$Count/1000
ggplot(dat_long, aes(x = index, y = Count, fill = Type)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = type_colors) +
  labs(
    x = NULL,
    y = "Number of transcripts (x10³)",
    fill = "Category"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x = element_text(
      angle = 50,
      hjust = 1,
      vjust = 1
    ))

# Fig.2c LR_SQANTI3_annot_coding_status -----------------------------------


dat<-fread("fig2c.LR_SQANTI3_annot_coding_status.txt")

dat_long <- melt(
  dat,
  id.vars = "index",
  variable.name = "Type",
  value.name = "Count"
)
dat_long$Type <- factor(
  dat_long$Type,
  levels = rev(c( "Protein-coding","Noncoding", "NMD-sensitive"))
)
dat_long$index <- factor(
  dat_long$index,
  levels = dat$index
)
type_colors <- c(
  "Protein-coding" = "#ac4630",
  "Noncoding"= "#e2a234",
  "NMD-sensitive"= "#559e87"
)
dat_long$Count<-dat_long$Count/1000

ggplot(dat_long, aes(x = index, y = Count, fill = Type)) +
  geom_col(width = 0.7) +
  scale_y_continuous(
    breaks = seq(0, 100, by = 20)
  ) +
  scale_fill_manual(values = type_colors) +
  labs(
    x = NULL,
    y = "Number of transcripts (x10³)",
    fill = "Coding status"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x = element_text(
      angle = 50,
      hjust = 1,
      vjust = 1
    ))


# Fig2.d LR_suppa_splicing_events_stat ------------------------------------

dat<-fread("fig2d.LR_suppa_splicing_events_stat.txt")
dat_long <- melt(
  dat,
  id.vars = "index",
  variable.name = "Type",
  value.name = "Count"
)
dat_long$Type <- factor(
  dat_long$Type,
  levels = rev(c( "GENCODE v47", "Novel"))
)
dat_long$index <- factor(
  dat_long$index,
  levels = dat$index
)
type_colors <- c(
  "GENCODE v47"        = "#9faac1",
  "Novel"          = "#9d3929"
)
dat_long$Count<-dat_long$Count/1000

ggplot(dat_long, aes(x = index, y = Count, fill = Type)) +
  geom_col(width = 0.7) +
  scale_y_continuous(
    breaks = seq(0, 400, by = 100)
  ) +
  scale_fill_manual(values = type_colors) +
  labs(
    x = NULL,
    y = "Number of splicing events (x10³)",
    fill = "Coding status"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x = element_text(
      angle = 50,
      hjust = 1,
      vjust = 1
    ))



# Fig2.e LR_transcript_tissue_specificity ---------------------------------
library(tidyr)
dat<-fread("fig2e.LR_transcript_tissue_specificity.txt")

dat_long <- pivot_longer(
  dat,
  cols = c("Annotated", "Novel"),
  names_to = "Type",
  values_to = "Proportion"
)

ggplot(dat_long, aes(x = index, y = Proportion, color = Type)) +
  geom_point(size = 2) +           
  geom_line(size = 1) +  
  scale_y_continuous(
    breaks = seq(0, 0.30, by = 0.05)
  ) +
  scale_color_manual(values = c("Annotated" = "#c3968f", "Novel" = "#9d3929")) +  
  labs(x = "Tissue number", y = "Proportion of transcripts", color = "Type") +
  theme_classic(base_size = 13)


# Fig2.f MS_all_peptide_validate ------------------------------------------
library(dplyr)
dat<-fread("fig2f.MS_all_peptide_validate.txt")
dat_long <- melt(
  dat,
  id.vars = "index",
  measure.vars = c("Unique peptide", "Shared peptide"),
  variable.name = "Type",
  value.name = "Count"
)
dat_long$Type <- factor(
  dat_long$Type,
  levels = rev(c( "Unique peptide", "Shared peptide"))
)
dat_long$index <- factor(
  dat_long$index,
  levels = dat$index
)
type_colors <- c(
  "Unique peptide"        = "#6b8ec4",
  "Shared peptide"          = "#abc0df"
)

tissue_colors <- setNames(dat[[4]], dat[[1]])


dat_long<-dat_long %>%
  left_join(dat[,c(1,4)],by="index")

point_data <- unique(dat_long[, .(index)]) 

ggplot(dat_long, aes(x = factor(index), y = Count, fill = Type)) +
  geom_col(width = 0.7) +   
  geom_point(
    data = point_data, 
    mapping = aes(x = factor(index), y = -0.05, color = index), 
    inherit.aes = FALSE, 
    size = 4
  ) +
  scale_fill_manual(values = type_colors) +
  scale_color_manual(values = tissue_colors) +
  scale_y_continuous(
    breaks = seq(0.0, 1.0, by = 0.2),
    expand = expansion(mult = c(0.05, 0.1))
  ) +
  labs(
    x = NULL,
    y = "Proportion of isoform with peptides",
    fill = "Type",
    color = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "top")

# Fig2.g WGCNA ------------------------------------------------------------

library(tidyverse)  
library(limma)  
library(RColorBrewer) 
library(vroom)
library(stringr)

load("heart_go_use.RData")
load("merged_expr.RData")
load("sample_annot_full_1586.RData")

heart<-sample_annot_full[str_detect(sample_annot_full$Tissue,"Heart"),]

tissue <- heart$Tissue
names(tissue) <- heart$sample_id


merged_expr_norm <- t(scale(t(merged_expr), center = TRUE, scale = TRUE))
merged_expr_norm <- as.data.frame(merged_expr_norm)

my_colors <- c(colorRampPalette(c("#8089a9", "white"))(50),
               colorRampPalette(c("white", "#bb4633"))(50))
my_breaks <- c(seq(-2, 0, length.out = 51), seq(0.01, 2, length.out = 50))

library(pheatmap)
p<-pheatmap(merged_expr_norm,
            color = my_colors,
            breaks = my_breaks,
            cluster_rows = T,
            cluster_cols = T,
            border_color = NA,
            show_rownames = T,
            show_colnames = T
)

p

row_order <- p$tree_row$order
row_names_sorted <- rownames(merged_expr_norm)[row_order]

#  go result --------------------------------------------------------------


go_result_bar$log10p<- -log10(go_result_bar$p.adjust)
library(ggplot2)
go_result_bar<-go_result_bar[match(row_names_sorted,go_result_bar$Module),]
go_result_bar$Description <- factor(go_result_bar$Description, levels = rev(go_result_bar$Description))

ggplot(go_result_bar, aes(x = Description, y = log10p)) +
  geom_bar(stat = "identity", width = 0.7, fill = "#56B4E9") +  
  coord_flip() +  
  labs(x = "GO Term", y = "-log10(p.adjust)", title = "GO Enrichment") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),          
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )


# proportion novel&annotated ----------------------------------------------


library(data.table)
temp_df<-fread("WGCNA_Heart_separate_geneInfo.tsv")

temp_df$type<-ifelse(str_detect(temp_df$gene_id,"GTOP"),"novel","annotated")
temp_df<-temp_df[temp_df$module_MEnumber %in% go_result_bar$Module,]
table(temp_df$type)
p_novel<-data.frame(table(temp_df$module_MEnumber,temp_df$type))

p_novel <- p_novel %>%
  group_by(Var1) %>%
  mutate(Proportion = Freq / sum(Freq)) %>%
  ungroup()
my_colors <- c("annotated" = "white", "novel" = "#CCB1CC")
p_novel$Var1<-factor(p_novel$Var1,levels = rev(row_names_sorted))

ggplot(p_novel, aes(x = Var1, y = Proportion, fill = Var2)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = my_colors) +  
  labs(
    x = "Var1",
    y = "Proportion",
    fill = "Type"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title.x = element_blank()
  )+
  coord_flip() 


# Fig2.h ASE/ASTS gene number  ---------------------------------------------------------
library(ggpubr)
df_plot <- fread("Fig 2h.txt")
ggplot(df_plot, aes(x=number, y=reorder(class, number)))+
  geom_bar(stat = "identity", fill="#5c86af")+
  #geom_text(aes(x=number+2, label=number))+
  theme_pubr()+
  ylab("")

# Fig.2i the number of significant ase/asts gene for per tissue -----------


df_count.ase_tissue <- fread("Fig 2i.txt")

ggplot(df_count.ase_tissue,aes(x=reorder(Tissue, -Freq),y=Freq)) + geom_bar(stat = "identity",fill="#a7b2c5") + theme_pubr() + 
  xlab("")+
  ylab("# ASE(gene)")+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = .5))

df_plot <- fread("Fig 2i.asts.txt")
df_plotsum <- df_plot %>% 
  group_by(Tissue) %>% 
  dplyr::summarise(total=sum(value)) %>% 
  arrange(desc(total))


df_plot$Tissue <- factor(df_plot$Tissue,levels = df_plotsum$Tissue)

ggplot(df_plot,aes(x=Tissue,y=value,fill=variable)) + geom_bar(stat = "identity") + theme_pubr() + 
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = .5)) + 
  scale_fill_manual(breaks = c("A","B","C"),values = c("#fee5d9","#fc9272","#de2d26"),
                    labels = c("no novel transcript", "with novel transcript", "with novel transcript aFC>1.5")) + 
  ylab("# ASTS(gene)")




