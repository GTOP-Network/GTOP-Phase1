
# SNP Analysis Pipeline: Allele Frequency Comparisons
# Date: 2025-06-21

library(data.table)
library(dplyr)
library(Biostrings)
library(ggplot2)
library(ggpubr)
library(waddR)

# SECTION 1: SNP FREQUENCY ANALYSIS (1000 Genomes) ---------------

#' Load and process GTOP SNP frequency data
#' Note: Original dataset contains ~7.4 million SNPs
process_gtop_snps <- function() {
  # Read small variant rsID mapping
  CHRPOS2rsid <- fread("2025-06-11-specific_xQTL/input/GTOP/SNV/GTOP_chrpos_rsid.txt",
                       col.names = c("chr_pos_ref_alt", "variant_id"))
  
  # Load allele frequencies from GTOP dataset
  GTOP_freq <- fread("2025-06-11-specific_xQTL/input/GTOP/SNV/GTOP_SNV.afreq")
  GTOP_freq <- merge(GTOP_freq[, .(variant_id = ID, af_GTOP = ALT_FREQS)], 
                     CHRPOS2rsid, by = "variant_id")[, .(chr_pos_ref_alt, variant_id, af_GTOP)]
  
  return(GTOP_freq)
}


#' Load and process 1000 Genomes population frequencies
process_1kg_snps <- function() {
  # populations <- c("AFR", "EUR", "SAS", "EAS", "AMR")
  # freq_list <- rbindlist(lapply(populations, function(pop) {
  #   file_path <- sprintf("/media/bora_A/zhangt/src/data/1000G/five_ancestry_groups/%s/%s.afreq", pop, pop)
  #   fread(file_path)[, Population := pop]
  # }))
  
  # Reshape to wide format (one row per SNP)
  OKG_freq <- fread("/media/bora_A/zhangt/src/data/1000G/five_ancestry_groups/all.afreq")
  colnames(OKG_freq)[2] <- "SNVID"
  colnames(OKG_freq)[5] <- "af_AFR"
  colnames(OKG_freq)[11] <- "af_AMR"
  colnames(OKG_freq)[17] <- "af_EAS"
  colnames(OKG_freq)[23] <- "af_EUR"
  colnames(OKG_freq)[29] <- "af_SAS"
  
  OKG_freq <- OKG_freq[, .(SNVID, af_AFR, af_AMR, af_EAS, af_EUR, af_SAS)] %>% 
    .[, chr_pos_ref_alt := paste0("chr", gsub(":", "_", SNVID))] %>%
    .[, c("chr", "pos") := tstrsplit(chr_pos_ref_alt, "_", keep = 1:2)] %>% 
    .[, .(chr_pos_ref_alt, chr, pos, af_AFR, af_EUR, af_SAS, af_EAS, af_AMR)]
  OKG_freq$pos <- as.numeric(OKG_freq$pos)
  
  OKG_freq <- OKG_freq[order(chr, as.numeric(pos))]
  
  return(OKG_freq)
}


# Main SNP processing pipeline -------------------------------------------------
SNV_freq_1KG <- process_1kg_snps()
fwrite(SNV_freq_1KG, "2025-06-11-specific_xQTL/output/data/processing/1KG_SNV_information.txt.gz", sep = "\t")


GTOP_freq <- process_gtop_snps() # 7374672
GTOP_freq[, c("chr", "pos", "ref", "alt") := 
            tstrsplit(chr_pos_ref_alt, "_", keep = 1:4)]
fwrite(GTOP_freq, "2025-06-11-specific_xQTL/output/data/processing/GTOP_SNV_information.txt.gz", sep = "\t")


## reload data
SNV_freq_1KG <- fread("2025-06-11-specific_xQTL/output/data/processing/1KG_SNV_information.txt.gz")
GTOP_freq <- fread("2025-06-11-specific_xQTL/output/data/processing/GTOP_SNV_information.txt.gz")
GTOP_freq$compare1KG <- "No"
GTOP_freq$compare1KG[GTOP_freq$chr_pos_ref_alt %in% SNV_freq_1KG$chr_pos_ref_alt] <- "Yes"
table(GTOP_freq$compare1KG)

fwrite(GTOP_freq, "2025-06-11-specific_xQTL/output/data/processing/GTOP_SNV_information.txt.gz")


ggvenn::ggvenn(list("GTOP"=GTOP_freq$chr_pos_ref_alt, "1KGP"=SNV_freq_1KG$chr_pos_ref_alt), show_percentage = F)

# Process 1000 Genomes data
GTOP_1KG_res <- merge(GTOP_freq, SNV_freq_1KG, by = c("chr_pos_ref_alt", "chr", "pos"))
GTOP_1KG_res$af_GTOP <- as.numeric(GTOP_1KG_res$af_GTOP)
GTOP_1KG_res$af_AFR <- as.numeric(GTOP_1KG_res$af_AFR)
GTOP_1KG_res$af_EUR <- as.numeric(GTOP_1KG_res$af_EUR)
GTOP_1KG_res$af_SAS <- as.numeric(GTOP_1KG_res$af_SAS)
GTOP_1KG_res$af_EAS <- as.numeric(GTOP_1KG_res$af_EAS)
GTOP_1KG_res$af_AMR <- as.numeric(GTOP_1KG_res$af_AMR)

GTOP_1KG_res$Type <- "C"
index1 <- (GTOP_1KG_res$af_AFR < 0.01 | GTOP_1KG_res$af_AFR > 0.99) |
  (GTOP_1KG_res$af_EUR < 0.01 | GTOP_1KG_res$af_EUR > 0.99) |
  (GTOP_1KG_res$af_SAS < 0.01 | GTOP_1KG_res$af_SAS > 0.99) |
  (GTOP_1KG_res$af_AMR < 0.01 | GTOP_1KG_res$af_AMR > 0.99)

index2 <- (GTOP_1KG_res$af_AFR < 0.01 | GTOP_1KG_res$af_AFR > 0.99) &
  (GTOP_1KG_res$af_EUR < 0.01 | GTOP_1KG_res$af_EUR > 0.99) &
  (GTOP_1KG_res$af_SAS < 0.01 | GTOP_1KG_res$af_SAS > 0.99) &
  (GTOP_1KG_res$af_AMR < 0.01 | GTOP_1KG_res$af_AMR > 0.99)
GTOP_1KG_res$Type[index1] <- "R1"
GTOP_1KG_res$Type[index2] <- "R2"
table(GTOP_1KG_res$Type)

table(abs(GTOP_1KG_res$af_GTOP-GTOP_1KG_res$af_EAS) < 0.2)

GTOP_1KG_res <- GTOP_1KG_res[abs(GTOP_1KG_res$af_GTOP-GTOP_1KG_res$af_EAS) < 0.2, ]
fwrite(GTOP_1KG_res, "2025-06-11-specific_xQTL/output/data/processing/GTOP_1KG_SNV_freq.txt.gz", sep = "\t")
