library(optparse)

option_list <- list(
	make_option(c("-a","--annotation"),type="character",default="NA",action="store",help="specify a tissue"),
	make_option(c("-t","--tissue"),type="character",default="NA",action="store",help="specify a tissue")
)

opt <- parse_args(OptionParser(option_list=option_list,usage="usage: %prog [options]"))

library(data.table)
library(dplyr)
library(magrittr)


# functions 
rm_version <- function(x){
	return(strsplit(x,split=".",fixed=T)[[1]][1])
}
get_tss_pos <- function(x){
	strand=x[4]
	if(strand=="+"){
		return(x[2])
	}else{
		return(x[3])
	}
}
# load annotation
df_ref <- fread(opt$annotation,header=F,sep="\t")
names(df_ref) <- c("chromosome","pos0","pos1","strand","gene","genename","genetype")

df_ref$TSS <- apply(df_ref,1,get_tss_pos)
# load gene list
genelist_file <- paste0(opt$tissue,".snv_egenes.FDR.05.txt")
df_genes <- fread(genelist_file,header=T,sep="\t")

df_gene_out <- df_ref %>% filter(gene %in% df_genes$Gene) %>% select(gene,chromosome,TSS)
df_gene_out$gene <- sapply(df_gene_out$gene,rm_version)

fwrite(df_gene_out,file=paste0("./input/",opt$tissue,".egene_list.txt"),quote=F,sep="\t",row.names=F,col.names=F)
