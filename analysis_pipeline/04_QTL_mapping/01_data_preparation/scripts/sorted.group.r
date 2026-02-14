library(data.table)
library(optparse)



option_list <- list(
  make_option(c("-p", "--phenotype"), 
              type = "character", 
              default = NULL,
              help = "Path to input BED file (required)",
              metavar = "FILE"),
  
  make_option(c("-g", "--group"), 
              type = "character", 
              default = NULL,
              help = "Path to group file, two columns without header: phenotype_id, group (required)",
              metavar = "FILE"),
  
  make_option(c("-w", "--wkpath"), 
              type = "character", 
              default = NULL,
              help = "Output directory path",
              metavar = "FILE"),
  
  make_option(c("-o", "--output"), 
              type = "character", 
              default = "sorted_output",
              help = "Output file prefix [default: %default]",
              metavar = "PREFIX"))
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)


setwd(opt$wkpath)
phenotype_bed <- fread(opt$phenotype, header = TRUE)
group_df <- fread(opt$group, header = FALSE, col.names = c("phenotype_id", "group"))


colnames(phenotype_bed)[1:4] <- c("#chr", "start", "end", "phenotype_id")




phenotype_bed[, chr_num := as.integer(fcase(
  `#chr` == "chr1", 1L,
  `#chr` == "chr2", 2L,
  `#chr` == "chr3", 3L,
  `#chr` == "chr4", 4L,
  `#chr` == "chr5", 5L,
  `#chr` == "chr6", 6L,
  `#chr` == "chr7", 7L,
  `#chr` == "chr8", 8L,
  `#chr` == "chr9", 9L,
  `#chr` == "chr10", 10L,
  `#chr` == "chr11", 11L,
  `#chr` == "chr12", 12L,
  `#chr` == "chr13", 13L,
  `#chr` == "chr14", 14L,
  `#chr` == "chr15", 15L,
  `#chr` == "chr16", 16L,
  `#chr` == "chr17", 17L,
  `#chr` == "chr18", 18L,
  `#chr` == "chr19", 19L,
  `#chr` == "chr20", 20L,
  `#chr` == "chr21", 21L,
  `#chr` == "chr22", 22L,
  `#chr` == "chrX", 23L,
  `#chr` == "chrM", 24L
))]


setorder(phenotype_bed, chr_num, start)


phenotype_bed <- merge(phenotype_bed, group_df, by = "phenotype_id", all.x = TRUE, sort = FALSE)


phenotype_bed[is.na(group), group := phenotype_id]
phenotype_bed[, group_index := .GRP, by = .(chr_num, group)]
setorder(phenotype_bed, chr_num, group_index, start)
group_sorted <- phenotype_bed[, .(phenotype_id, group)]
fwrite(group_sorted, "phenotype_groups.sorted.txt", sep = "\t", col.names = FALSE)


cols_to_remove <- c("chr_num", "group_index", "group")
output_cols <- setdiff(names(phenotype_bed), cols_to_remove)


output_bed <- phenotype_bed[, .SD, .SDcols = output_cols]
setcolorder(output_bed, c("#chr", "start", "end", "phenotype_id"))

fwrite(output_bed, "phenotype.sorted.bed", sep = "\t", quote = FALSE)
