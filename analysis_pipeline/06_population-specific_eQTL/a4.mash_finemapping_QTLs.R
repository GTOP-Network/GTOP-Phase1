
library(mashr)

"%&%" = function(a, b) { paste0(a, b) }
readFile = function(f){
  if(grepl(".RDS", f, ignore.case=T)){
    DF = readRDS(f)
  } else{
    DF = fread(f)
  }
  return(DF)
}


############ strong pairs
raw.zval_data.strong = as.data.frame(readFile("2025-09-30-mash/output/SNV_eQTL/Strong2/strong_pairs.MashR_input.txt.gz"))
rownames(raw.zval_data.strong) = raw.zval_data.strong$pair_id
raw.zval_data.strong$pair_id = NULL

raw.zval_data.strong[is.na(raw.zval_data.strong)] = 0
zval_data.strong = as.matrix(raw.zval_data.strong)

############ all pairs of a random subset of 1000000 tests
zval_data.random.subset = as.data.frame(readFile("2025-09-30-mash/output/SNV_eQTL/Random/MashR.random_subset_1000000.RDS"))
rownames(zval_data.random.subset) = zval_data.random.subset$pair_id
zval_data.random.subset$pair_id = NULL
zval_data.random.subset = as.matrix(zval_data.random.subset)

zval_data.random.subset <- zval_data.random.subset[, colnames(zval_data.strong)]
if (sum(colnames(zval_data.random.subset) != colnames(zval_data.strong))>0){
  stop("Different colnames between random subset and strong subset.")
}

## Correlation structure
data.temp = mash_set_data(Bhat=zval_data.random.subset,alpha=1,zero_Bhat_Shat_reset=1e6)
Vhat = estimate_null_correlation_simple(data.temp)
rm(data.temp)

data.strong = mash_set_data(Bhat=zval_data.strong,alpha=1,V=Vhat,zero_Bhat_Shat_reset=1e6)

## Compute posterior summaries
m.r <- readRDS("2025-09-30-mash/output/SNV_eQTL/top_pairs/m.r.RDS")

finemappint.m.s <- mash(data.strong, g=get_fitted_g(m.r), fixg=TRUE)

saveRDS(finemappint.m.s, "2025-09-30-mash/output/SNV_eQTL/top_pairs2/m.s_zval.RDS")
saveRDS(get_lfsr(finemappint.m.s), "2025-09-30-mash/output/SNV_eQTL/top_pairs2/lfsr_m.s.RDS")
