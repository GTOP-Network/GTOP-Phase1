#==============================================#
# GWAS-fdQTL #
# Supp-Figure-28#
#==============================================#


library(ggplot2)
library(tidyverse)
library(ggpubr)

setwd("/media/london_A/mengxin/GTOP_code/supp/supp_fig28")


# Supp.Fig.28a Example Disease associated-fd-QTL --------------------------------------------

## Figure a
load("./input/supp_fig28a_data.RData")

loci_start <- min(GWAS_data$pos)
loci_end <- max(GWAS_data$pos) #- 500000
chr_name <- GWAS_data$chrom[1]
SNP_name <- "rs3782886" #GWAS_data$rsid[which.min(GWAS_data$p)]
GWAS_name <- "Weight"

## GWAS locuszoom
GWAS_locus <- locus(data = as.data.frame(GWAS_data),
                    xrange = c(loci_start, loci_end),
                    seqname = chr_name, index_snp = SNP_name,
                    ens_db = "EnsDb.Hsapiens.v86")
GWAS_locus$data$ld <- d1$LD[,SNP_name][GWAS_locus$data$rsid]^2
GWAS_plot <- gg_scatter(GWAS_locus, pcutoff = FALSE, yzero=T,  labels = SNP_name, legend_pos = "right",
                        LD_scheme = c("#e5e5e5", "#e5e5e5", "#3e70b4", "#3f7d1d",
                                      "orange", "red", "red")) + 
  annotate("text",x=loci_end/10^6, y=max(-log10(GWAS_data$p)),label=GWAS_name, hjust="right") +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank())

## eQTL locuszoom
SNV_eQTL_locus <- locus(data = as.data.frame(SNV_eQTL_data),
                        xrange = c(loci_start, loci_end),
                        seqname = chr_name, index_snp = SNP_name,
                        ens_db = "EnsDb.Hsapiens.v86")
SNV_eQTL_locus$data$ld <- d2$LD[,SNP_name][SNV_eQTL_locus$data$rsid]^2

SNV_eQTL_plot <- gg_scatter(SNV_eQTL_locus, pcutoff = FALSE, yzero=T,  labels = SNP_name, legend_pos = "right",
                            LD_scheme = c("#e5e5e5", "#e5e5e5", "#3e70b4", "#3f7d1d",
                                          "orange", "red", "red")) + 
  annotate("text",x=loci_end/10^6, y=max(-log10(SNV_eQTL_data$p)),label="eQTL_Whole Blood", hjust="right") +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank())

## gene structure
gene_plot <- gg_genetracks(SNV_eQTL_locus, highlight = "BRAP", 
                           filter_gene_name = c("BRAP"), 
                           filter_gene_biotype = c("protein_coding"))

## plot
wrap_plots(list(GWAS_plot, SNV_eQTL_plot, gene_plot), ncol = 1, heights = c(3,3,3,3,1))



# Supp.Fig.28b ------------------------------------------------------------


tmp_df <- fread("./input/supp_fig28b_data.txt")

ggplot(tmp_df[tmp_df$p<0.01,], aes(af_other*100, af_GTOP*100)) +
  ggpointdensity::geom_pointdensity( alpha=1, aes(size=-log10(p)), adjust=1, shape=21 )+
  scale_color_viridis_b()+
  ggdensity::geom_hdr_lines( linetype="dashed", linewidth=0.5)+
  labs(x="Alternative allele frequency (Other populations)", 
       y="Alternative allele frequency (GTOP)")+
  facet_grid(.~type) +
  theme_classic()+
  theme(
    axis.line = element_line(color="black"),
    axis.ticks = element_line(color="black"),
    axis.text = element_text(color="black"),
    legend.position = "right",
  ) + ylim(c(0,100)) +
  ggh4x::coord_axes_inside(labels_inside = F)
