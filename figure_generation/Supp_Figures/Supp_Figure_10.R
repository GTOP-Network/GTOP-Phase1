#===============================================#
# MS peptide  #
# Supp-Figure-10 #
#===============================================#
library(ggplot2)
library(ggpubr)
library(tidyverse)

setwd("/path/to/GTOP_code/supp/supp_fig10/input")

# Supp.Fig.10a ------------------------------------------------------------

plot_df <- fread("supp10a.MS_peptide_length_dist.txt")

ggplot(plot_df)+
  geom_bar(aes(x=peptide_length, y=count), stat = "identity", fill="#2171a9")+
  xlab("Peptide length (aa)")+
  ylab("Number of peptides")+
  theme_pubr()


# Supp.Fig.10b ------------------------------------------------------------

library(dplyr)
dat<-fread("supp10b.MS_all_peptide_validate.annotated_transcript.txt")
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


# Supp.Fig.10c --------------------------------------------------------------

df <- fread("supp10c.MS_unique_peptide_support_novel_transcript.txt")
plot_df <- df %>% 
  pivot_longer(cols = -c("tissue"), names_to = "class",
               values_to = "number")
plot_df %>% 
  mutate(class=factor(class, levels=c("> 5", "2 - 5", "1"))) %>% 
  ggplot(.)+
  geom_bar(aes(x=tissue, y=number, fill=class), stat = "identity", width = 0.6)+
  scale_fill_manual(nam = "# uniquely supported peptides",values = c("1"="#cedbeb", "2 - 5"="#82aece", "> 5"="#4174ad"))+
  xlab("")+
  ylab("Number of novel transcripts")+
  theme_pubr()+
  theme(axis.text.x = element_text(angle = 55, hjust = 1))
