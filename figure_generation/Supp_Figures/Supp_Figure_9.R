#===============================================#
# LRS read QC #
# Supp-Figure-9                         #
#===============================================#
library(ggplot2)
library(ggpubr)
library(ComplexUpset)
library(tidyverse)
library(data.table)

setwd("/path/to/GTOP_code/supp/supp_fig9/input")


# Supp.Fig.9a Candidate transcripts --------------------------------------------

df <- read_tsv("supp9a.LR_candidate_isoform_tool_upset.txt",
  show_col_types = FALSE
)
df <- df %>%
  rename(
    Bambu   = bambu,
    FLAIR   = flair,
    `Iso-Seq` = isoseq
  )
tool_cols <- c("Bambu", "FLAIR", "Iso-Seq")

upset(
  df,
  intersect = tool_cols,
  base_annotations = list(
    "Intersection size" = intersection_size(
      counts  = FALSE
    )
  ),
  
  width_ratio = 0.3
)  +
  theme(
    text = element_text(size = 12, color = "black"),
    axis.text  = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    axis.line  = element_line(color = "black", linewidth = 0.6),
    panel.background = element_blank(),
    panel.grid = element_blank(), 
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 11),
    strip.background = element_blank(),
    strip.text = element_text(color = "black")
  )


# Supp.Fig.9c FLNC read counts --------------------------------------------

# Load data
df <- fread("supp9c.LR_isoform_QC_reads.txt.gz")
flnc_reads <- df[, rowSums(.SD), .SDcols = patterns("^")]
log_flnc_reads <- log10(flnc_reads)

original_threshold <- 10
log_threshold <- log10(original_threshold)  # exactly 1.0

# Create bins
max_log_val <- ceiling(max(log_flnc_reads, na.rm = TRUE) * 10) / 5
custom_bins <- seq(0, max_log_val + 0.2, by = 0.2)

bin_table <- tibble(
  bin_left  = custom_bins[-length(custom_bins)],
  bin_right = custom_bins[-1]
)

binned_counts <- tibble(log_flnc = log_flnc_reads) %>%
  mutate(
    bin = cut(
      log_flnc,
      breaks = custom_bins,
      include.lowest = TRUE,
      right = FALSE
    )
  ) %>%
  count(bin)

plot_df <- bin_table %>%
  mutate(
    bin = factor(
      paste0("[", bin_left, ",", bin_right, ")"),
      levels = levels(binned_counts$bin)
    )
  ) %>%
  left_join(binned_counts, by = "bin") %>%
  mutate(
    n = ifelse(is.na(n), 0, n),
    pass = ifelse(bin_left >= log_threshold, "passed", "filtered")
  )
plot_df$pass <- factor(plot_df$pass, levels = c("filtered", "passed"))
legend_labels <- c(
  "filtered" = paste0("Filtered out (< ", original_threshold, ")"),
  "passed"   = paste0("Passed filters (≥ ", original_threshold, ")")
)

# Plot
ggplot(plot_df) +
  geom_rect(
    aes(
      xmin = bin_left,
      xmax = bin_right,
      ymin = 0,
      ymax = n,
      fill = pass
    ),
    color = "black",
    linewidth = 0.4
  ) +
  geom_vline(
    xintercept = log_threshold,
    linetype = "dashed",
    linewidth = 1,
    color = "red"
  ) +
  scale_fill_manual(
    values = c(
      "filtered" = "#FFA07A",
      "passed"   = "#4682B4"
    ),
    labels = legend_labels,
    drop = FALSE
  ) +
  scale_x_continuous(
    limits = c(0,log10(max(flnc_reads))+0.2),
    breaks = seq(0, max_log_val + 1, by = 1),
    expand = c(0, 0)
  ) +
  labs(
    x = "Supported FLNC read count (log10 scale)",
    y = "Number of transcripts",
    fill = NULL
  ) +
  theme_pubr() +
  theme(
    legend.position = "top"
  )


# Supp.Fig.9d -------------------------------------------------------------


plot_df <- fread("supp9d.LR_isoform_QC_support_samples.txt.gz")
binwidth <- (max(plot_df$support_samples) - min(plot_df$support_samples)) / 20
ggplot(plot_df, aes(x = support_samples)) +
  geom_histogram(binwidth = binwidth, 
                 fill = "#7f8080", center = min(plot_df$support_samples) + binwidth/2,
                 color = "black") +
  xlab("Number of supporting samples") +
  ylab("Frequency") +
  theme_pubr()



# Supp.Fig.9e -------------------------------------------------------------


plot_df <- fread("supp9e.LR_read_QC_intra_priming.txt")


ggplot(plot_df) +
  geom_bar(
    aes(x = A_number, y = frequency, fill = type),
    stat = "identity",
    color = "black",
    size = 0.2, 
    width = 1     
  ) +
  scale_fill_manual(values = setNames(plot_df$color, plot_df$type)) +
  xlab("# A in 20bp downstream TTS") +
  ylab("Number of transcripts") +
  theme_pubr() +
  theme(legend.position = "none")


# Supp.Fig.9f -------------------------------------------------------------


plot_df <- fread("supp9f.LR_read_QC_pass_category.Passed filters.txt")
color <- plot_df %>% distinct(type, color)

plot_df$category <- factor(plot_df$category, levels = plot_df$category)

ggplot(plot_df)+
  geom_bar(aes(x=category, y=count, fill=type), stat = "identity", color = "black")+
  scale_fill_manual(values=setNames(color$color, color$type))+
  xlab("")+
  ylab("Transcript count")+
  theme_pubr()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")


# Supp.Fig.9g -------------------------------------------------------------


plot_df <- fread("supp9g.LR_read_QC_pass_category.Filtered out.txt")
color <- plot_df %>% distinct(type, color)
plot_df$category <- factor(plot_df$category, levels = plot_df$category)

ggplot(plot_df)+
  geom_bar(aes(x=category, y=count, fill=type), stat = "identity", color = "black")+
  scale_fill_manual(values=setNames(color$color, color$type))+
  xlab("")+
  ylab("Transcript count")+
  theme_pubr()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")


# Supp.Fig.9h -------------------------------------------------------------

df <- read_tsv(
  "supp9h.LR_final_isoform_tool_with_cate_upset.txt",
  show_col_types = FALSE
)

df <- df %>%
  rename(
    Bambu   = bambu,
    FLAIR   = flair,
    `Iso-Seq` = isoseq
  )
tool_cols <- c("Bambu", "FLAIR", "Iso-Seq")

sorted_cate <- c('FSM', 'ISM', 'NIC', 'NNC', 'Fusion', 'Antisense', 'Genic', 'Intergenic')

df <- df %>%
  mutate(
    struct_class = factor(
      struct_class,
      levels = c(
        intersect(sorted_cate, unique(struct_class)),
        setdiff(unique(struct_class), sorted_cate)
      )
    )
  )

fill_colors <- c(
  "FSM"       = "#226ca0",
  "ISM"       = "#f39800",
  "NIC"       = "#65a066",
  "NNC"       = "#00a29a",
  "Fusion"    = "#595757",
  "Antisense" = "#e94651",
  "Genic"     = "#7861a8",
  "Intergenic"= "#c17272"
)


upset(
  df,tool_cols,
  base_annotations=list(
    'Intersection size'=intersection_size(
      counts=FALSE,
      mapping=aes(fill=struct_class)
    ) + scale_fill_manual(values=fill_colors)
  ),
  width_ratio=0.3
)  +
  theme(
    text = element_text(size = 12, color = "black"),
    axis.text  = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    axis.line  = element_line(color = "black", linewidth = 0.6),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 11),
    strip.background = element_blank(),
    strip.text = element_text(color = "black")
  )


# Supp.Fig.9i -------------------------------------------------------------

plot_scatter_with_fit <- function(file, x_col){
  
  df <- read.table(file,
                   header = TRUE,
                   sep = "\t",
                   check.names = FALSE)
  
  y_col <- "Number of tissue-specific transcripts"
  group_col <- "annotation"
  
  df <- df %>%
    select(all_of(c(x_col, y_col, group_col))) %>%
    na.omit()
  
  stats <- df %>%
    group_by(.data[[group_col]]) %>%
    summarise(
      p = cor.test(
        .data[[x_col]],
        .data[[y_col]],
        method = "spearman",
        exact = FALSE
      )$p.value
    )
  
  stats$label <- paste0(stats[[group_col]],
                        ": p=",
                        sprintf("%.2f", stats$p))
  novel_annotated_COLOR <- c(
    "Annotated" = "#cf928f",
    "Novel" = "#ad3b2b"
  )
  
  ggplot(df,
         aes(x = .data[[x_col]],
             y = .data[[y_col]],
             color = .data[[group_col]])) +
    geom_point(size = 2.8) +
    geom_smooth(method = "lm",
                se = FALSE,
                linetype = "dashed",
                linewidth = 0.7) +
    scale_color_manual(values = novel_annotated_COLOR) +
    annotate("text",
             x = min(df[[x_col]]),
             y = max(df[[y_col]]),
             label = paste(stats$label, collapse = "\n"),
             hjust = 0,
             vjust = 1,
             size = 4) +
    labs(
      x = x_col,
      y = y_col,
      color = "Annotation"
    ) +
    theme_classic(base_size = 12) +
    theme(legend.position = "top")
}


plot_scatter_with_fit(
  file = "supp9i-k.LR_tissue_specific.txt",
  x_col = "Sample size"
)


# Supp.Fig.9j -------------------------------------------------------------

plot_scatter_with_fit(
  file = "supp9i-k.LR_tissue_specific.txt",
  x_col = "Median number of reads"
)


# Supp.Fig.9k -------------------------------------------------------------

plot_scatter_with_fit(
  file = "supp9i-k.LR_tissue_specific.txt",
  x_col = "Median RIN"
)


# Supp.Fig.9l -------------------------------------------------------------

library(tidyr)
dat<-fread("supp9l.LR_transcript_tissue_specificity_10samples.txt")

dat_long <- pivot_longer(
  dat,
  cols = c("Annotated", "Novel"),
  names_to = "Type",
  values_to = "Proportion"
)

ggplot(dat_long, aes(x = index, y = Proportion, color = Type)) +
  geom_point(size = 2) +           
  geom_line(linewidth = 1) +  
  scale_y_continuous(
    breaks = seq(0, 0.50, by = 0.05)
  ) +
  scale_color_manual(values = c("Annotated" = "#c3968f", "Novel" = "#9d3929")) +  
  labs(x = "Tissue number", y = "Proportion of transcripts", color = "Type") +
  theme_classic(base_size = 13)



