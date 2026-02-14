library(optparse)

option_list <- list(
	make_option(c("-r","--ref"),type="character",default="NA",action="store",help="specify a tissue"),
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
# load reference
df_ref <- fread("./input/reference/GTOP.SNP_list.txt",header=F,sep="\t")
names(df_ref) <- c("chrom","bp","SNP","A2","A1")

# load gene list
genelist_file <- paste0("./input/",opt$tissue,".egene_list.txt")
df_genes <- fread(genelist_file,header=F,sep="\t")
names(df_genes) <- c("Gene","chrom","pos")

#selected_genes <- sapply(unique(df_genes$Gene),rm_version)
selected_genes <- unique(df_genes$Gene)

for(i in 1:length(selected_genes)){
	cat(i,selected_genes[i]," \n")
	infile <- paste0("./output/",opt$tissue,"/SuSiEx/examples/GTOP_",selected_genes[i],".txt")
	if(file.exists(infile)){
		df_gene <- fread(infile,header=T,sep="\t")
		df_ref.sub <- df_ref %>% filter(SNP %in% df_gene$SNP)
		df_gene <- merge(df_gene,df_ref.sub,by="SNP")
		df_gene %<>% select(chrom,bp,A2,A1,SNP,pvalue,beta,se)
		output_file <- paste0("./output/",opt$tissue,"/SuSiEx/examples/GTOP_",selected_genes[i],".sumstat.txt")
		fwrite(df_gene,file=output_file,quote=F,sep="\t",row.names=F)
	}else{
		cat(infile," not exists!!\n")
	}
}

