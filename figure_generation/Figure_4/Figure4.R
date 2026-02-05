#==================================#
# Population-specific QTLs #
# Figure-4 #
#==================================#

setwd(".")

library(data.table)
library(tidyverse)
library(reticulate)
library(locuszoomr)
library(EnsDb.Hsapiens.v86)
library(patchwork)


## Fig. 4a, Geographic frequencies of eQTL -------------------------------------

if (!py_module_available("matplotlib")) {
  py_install("matplotlib")
}

if (!py_module_available("geovar")) {
  py_install("git+https://github.com/aabiddanda/geovar")
}

py_run_string('
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# import pkg_resources
from geovar import *

plt.rcParams[\'pdf.fonttype\'] = 42

geovar_test = GeoVar()
geovar_test.add_freq_mat("./input/Fig4a.txt")
geovar_test.geovar_binning()

geovar_plot = GeoVarPlot()
geovar_plot.add_data_geovar(geovar_test)
geovar_plot.filter_data()
geovar_plot.add_cmap()

fig, ax = plt.subplots(1,1,figsize=(3,6))
geovar_plot.plot_geovar(ax)
ax.set_xticklabels(geovar_plot.poplist)
plt.savefig("./Fig4a.pdf", dpi=300, bbox_inches=\'tight\')
')


## Fig 4b, Allele frequency difference -----------------------------------------

freq_plot_data <- fread("./input/Fig4b.txt")

ggplot(freq_plot_data, aes(x = af_GTEx, y = af_GTOP)) +
  geom_point(aes(color = proxy_delta_group), alpha = 0.6) +
  scale_color_manual(
    name = "Frequency Difference",
    values =  c("#f8e4da", "#f1cf9c", "#e6986d", "#cf5e47", "#a43331", "#5c95d6") # , "#475c92"
  ) +
  labs(
    x = "Allele frequency (1KGP European)",
    y = "Allele frequency (GTOP)"
  ) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  theme_classic() +
  theme(
    text = element_text(family = "Arial"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black")
  ) +
  geom_hline(yintercept=c(0.05, 0.95), linetype='dashed', color='grey', size=0.5) +
  geom_abline(linetype='dashed', color='darkgrey', size=0.5)

length(freq_plot_data$variant_proxy[freq_plot_data$proxy_delta_group%in%
                                      c("|ΔAF| > 0.2", "|ΔAF| > 0.3", 
                                        "|ΔAF| > 0.4", "|ΔAF| > 0.5")])/
  length(freq_plot_data$variant_proxy)


## Fig 4c, Distribution of fd-QTL variants -------------------------------------

freq_for_type <- fread("./input/Fig4c.txt")

ggplot(data = freq_for_type, aes(x=proxy_delta_group, y=value, fill=variable)) +
  geom_col(position = 'dodge') +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black", angle = 30, vjust = 1, hjust = 1),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  scale_fill_manual(values = c("#0f3d7b", "#7d8cad")) +
  labs(x="", y="Proportion of fine-mapped variants")


## Fig 4d, Distribution of fd-QTL genes ----------------------------------------

ff_gene_summary <- fread("./input/Fig4d.txt")
color_vec <- readRDS("./input/tissue_color.RDS")

ggplot(ff_gene_summary, aes(x = Gene_count, y = new_two_group)) +
  geom_boxplot(outlier.colour = NA, aes(fill = new_two_group), width=0.8) +
  ggbeeswarm::geom_quasirandom(aes(group=tissue, col=tissue), size=2, width = 0.4) +
  scale_color_manual(values = color_vec) +
  scale_fill_manual(values = c("#c0c3d1", "#959db7")) +
  theme_classic() +
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(x="", y="Gene number")


## Fig 4e, fd-eQTL example -----------------------------------------------------
load("./input/Fig4e.RData")

SNP_name <- GWAS_data$rsid[which.min(GWAS_data$p)]
SV_name <- GWAS_data$rsid[which.min(GWAS_data$p)]
loci_start <- min(GWAS_data$pos) + 200000
loci_end <- max(GWAS_data$pos) - 400000
GWAS_name <- "Urticaria"
chr_name <- GWAS_data$chrom[1]

GWAS_locus <- locus(data = as.data.frame(GWAS_data),
                    xrange = c(loci_start, loci_end),
                    seqname = chr_name, index_snp = SNP_name,
                    ens_db = "EnsDb.Hsapiens.v86")
GWAS_locus <- link_LD(GWAS_locus, token = "b6336e5da5d3", pop = "EAS")
GWAS_plot <- gg_scatter(GWAS_locus, pcutoff = FALSE, yzero=T,  labels = SNP_name, legend_pos = "right",
                        LD_scheme = c("#e5e5e5", "#e5e5e5", "#3e70b4", "#3f7d1d",
                                      "orange", "red", "red")) + 
  annotate("text",x=loci_end/10^6, y=max(-log10(GWAS_data$p)),label=GWAS_name, hjust="right") +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank()) +
  ylim(c(0, max(-log10(GWAS_data$p)+1))) +
  geom_hline(yintercept = -log10(5e-8), linetype = "dashed", size = 1) 

SNV_eQTL_locus <- locus(data = as.data.frame(SNV_eQTL_data),
                        xrange = c(loci_start, loci_end),
                        seqname = chr_name, index_snp = SNP_name,
                        ens_db = "EnsDb.Hsapiens.v86")
SNV_eQTL_locus <- link_LD(SNV_eQTL_locus, token = "b6336e5da5d3", pop = "EAS")
SNV_eQTL_plot <- gg_scatter(SNV_eQTL_locus, pcutoff = FALSE, yzero=T,  labels = SNP_name, legend_pos = "right",
                            LD_scheme = c("#e5e5e5", "#e5e5e5", "#3e70b4", "#3f7d1d",
                                          "orange", "red", "red")) + 
  annotate("text",x=loci_end/10^6, y=max(-log10(SNV_eQTL_data$p)),label="eQTL_Whole Blood", hjust="right") +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank())

gene_plot <- gg_genetracks(SNV_eQTL_locus, highlight = "STIM1",
                           filter_gene_name = c("STIM1"),
                           filter_gene_biotype = c("protein_coding"))

wrap_plots(list(GWAS_plot, SNV_eQTL_plot, gene_plot), ncol = 1, heights = c(3,3,1))


## Fig 4f he-QTL ---------------------------------------------------------------
count_df <- fread("./input/Fig4f.txt")
color_vec <- readRDS("./input/tissue_color.RDS")

ggplot(count_df, aes(count, reorder(tissue, count), color=tissue)) +
  ylab('') + xlab('Number of he-QTLs') +
  geom_segment(aes(x=0, xend=count, y=tissue, yend=tissue)) +
  geom_point(shape=16, size=4) +
  scale_color_manual(values = color_vec)+
  theme_classic(base_size = 12) +
  theme(text=element_text(size=12),
        legend.position = "none",
        legend.title = element_blank(),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"))


## Fig 4g Z-score of he-QTL between GTOP and GTEx ------------------------------
count_df <- fread("./input/Fig4g.txt")

ggplot(data = plot_df, aes(x = type1, y = abs(value), fill=type1)) + 
  geom_violin() +
  geom_boxplot(width=0.1, outliers = FALSE) +
  theme_classic() +
  scale_fill_manual(values = c("GTEx"="#227377", "GTOP"="#8c3426")) +
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(x="", y="|Z-score|")

