#==============================================#
# GTOP-eQTL-correlation-SNV #
# Supp-Figure-21#
#==============================================#
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(data.table)
library(corrplot)

setwd("/path/to/GTOP_code/supp/supp_fig21")
# Supp.Fig.21a correlation of eQTL effects between pancreas -----------

dat.w <- fread("./input/Figure S21a.cov.txt")
rownames(dat.w) <- dat.w$replication_tissue
dat.w$replication_tissue <- NULL
corrplot(as.matrix(dat.w),method = "number",is.corr = F,col = COL1('Blues', 200)[50:200],
         type = "upper",diag = F)

scatterdf <- fread("./input/Figure S21a.scatter.txt")
df_tissue_pairs <- scatterdf %>% distinct(rep, dis)

scatter_list <- list()

for(i in 1:dim(df_tissue_pairs)[1]){
  cat(i,": ",df_tissue_pairs$dis[i],", ",df_tissue_pairs$rep[i],"\n")
  
  tissuePair <- paste0(df_tissue_pairs$dis[i],".",df_tissue_pairs$rep[i])
  
  dat <- scatterdf %>% dplyr::filter(rep==df_tissue_pairs$rep[i], dis==df_tissue_pairs$dis[i])
  p <- ggplot(dat,aes(x=Beta,y=slope)) + geom_point() + theme_pubr() + ggtitle(label = tissuePair)
  scatter_list[[i]] <- p
}

cowplot::plot_grid(plotlist = scatter_list,ncol = 4,align = "hv")


# Supp.Fig.21b rb of eQTL between tissues -----------

df <- readRDS("input/Figure S22b.rds")

ggplot(df, aes(x = tis1, y = tis2, fill = r_b)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = label), 
            size = 3, 
            color = "black") +
  scale_fill_gradient2(
    high = "#2171B5",
    mid = "#9ECAE1",
    low = "white",
    midpoint = median(df$r_b, na.rm = TRUE),
    name = "Rb"
  )+
  labs(
    x = "Discovery",
    y = "Replication",
  ) +
  theme_pubr() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(hjust = 1),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  ) +
  coord_fixed(ratio = 1)

