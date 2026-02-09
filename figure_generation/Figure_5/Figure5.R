#==============================================#
# LRS WGS QC #
# Figure-5#
#==============================================#
library(data.table)
library(dplyr)
library(magrittr)
library(ggplot2)
library(ggpubr)

setwd("/media/london_A/mengxin/GTOP_code/fig-5")

# Fig.5a ------------------------------------------------------------------
# load data
df_sv_eqtl.r2 <- readRDS("./input/dat_sv.RDS")
df_tr_eqtl.r2 <- readRDS("./input/dat_tr.RDS")
df_tissue_color <- readRDS("./input/tissue_colors.RDS")

# poorly tagged eSV: LD R2<0.8
df_sv.poor <- df_sv_eqtl.r2 %>% filter(R2<0.8) # 21.9%
df_tr.poor <- df_tr_eqtl.r2 %>% filter(R2<0.8) # 33.8%

# no. of poorly tagged SV-eQTL and TR-eQTL
median(table(df_sv.poor$Tissue))# :
median(table(df_tr.poor$Tissue))

df_sv_eqtl.r2$R2_pseu <- df_sv_eqtl.r2$R2
df_sv_eqtl.r2$R2_pseu[df_sv_eqtl.r2$R2<0.2] <- 0.2
df_tr_eqtl.r2$R2_pseu <- df_tr_eqtl.r2$R2
df_tr_eqtl.r2$R2_pseu[df_tr_eqtl.r2$R2<0.2] <- 0.2


p1 <- ggplot(df_sv_eqtl.r2,aes(x=R2_pseu,fill=Tissue)) + geom_histogram(position = "stack") + theme_pubr() + 
  scale_fill_manual(breaks = df_tissue_color$Tissue,values = paste0("#",df_tissue_color$Tissue_Color_Code)) + 
  scale_x_continuous(breaks = c(0.2,0.4,0.6,0.8,1.0));p1
p2 <- ggplot(df_tr_eqtl.r2,aes(x=R2_pseu,fill=Tissue)) + geom_histogram(position = "stack") + theme_pubr() + 
  scale_fill_manual(breaks = df_tissue_color$Tissue,values = paste0("#",df_tissue_color$Tissue_Color_Code)) + 
  theme(legend.position = "none") + scale_x_continuous(breaks = c(0.2,0.4,0.6,0.8,1.0));p2
cowplot::plot_grid(p1+theme(legend.position = "none"),p2,ncol=1,align = "v")


# Fig.5b ------------------------------------------------------------------


df_cs.nr <- readRDS("./input/joint_finemap.res.RDS")

df_var_count <- as.data.frame(as.matrix.data.frame(table(df_cs.nr$Tissue,df_cs.nr$varType)))
df_var_count$Tissue <- unique(df_cs.nr$Tissue)
df_var_count <- df_var_count[order(-df_var_count$V1),]

df_plot_count <- as.data.frame(table(df_cs.nr$Tissue,df_cs.nr$varType))
df_plot_count$Var1 <- as.character(df_plot_count$Var1)
df_plot_count$Var2 <- as.character(df_plot_count$Var2)

df_plot_count$Var1 <- factor(df_plot_count$Var1,levels = df_var_count$Tissue)
df_plot_count$Var2 <- factor(df_plot_count$Var2,levels = c("SNV","TR","SV"))

p3 <- ggplot(df_plot_count,aes(x=Var1,y=log2(Freq))) + geom_bar(stat = "identity",width=.6,aes(fill = Var2),position = position_dodge()) + theme_pubr() + 
  scale_fill_manual(breaks = c("SNV","TR","SV"),values = c("#8090B4","#963628","#227E84")) + 
  theme(axis.text.x = element_text(angle = 90,hjust = 1, vjust = .5));p3


# Fig.5c ------------------------------------------------------------------


df_cs.nr <- readRDS("./input/non_redundant_cs.RDS")
df_prp <- readRDS("./input/prp_svtr_lead.RDS")

df_cs_count <- as.data.frame(table(df_cs.nr$Tissue))
names(df_cs_count) <- c("Tissue","Count")
df_cs_count <- df_cs_count[order(-df_cs_count$Count),]

library(ggplot2)
library(ggpubr)
df_prp$Var1 <- factor(df_prp$Var1,levels = df_cs_count$Tissue)
df_prp$Var2 <- factor(df_prp$Var2,levels = c(0,1,2))
p4 <- ggplot(df_prp,aes(x=Var1,y=Freq,fill=Var2)) + geom_bar(stat = "identity",position = "stack") + theme_pubr() + 
  scale_fill_manual(breaks = c(0,1,2),values = c("#7D8BAD","#7D8BAD","#7D8BAD")) + 
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = .5)) + 
  xlab("Tissues") + ylab("# of credible sets");p1

library(reshape2)
df_prp.w <- dcast(df_prp,Var1 ~ Var2, value.var = "Freq")
names(df_prp.w) <- c("Tissue","SNV","SV_TR","SV_TR_lead")
df_prp.w$total <- apply(df_prp.w[,-1],1,sum)

df_prp.w %<>% mutate(PC_snv = SNV/total, PC_sv_tr=SV_TR/total, PC_sv_trLead=SV_TR_lead/total) %>% select(Tissue,PC_snv,PC_sv_tr,PC_sv_trLead)
#write.table(df_prp.w,file = )
df_prp.wl <- melt(df_prp.w,id.vars = "Tissue")
df_prp.wl$Tissue <- as.character(df_prp.wl$Tissue)
df_prp.wl$Tissue <- factor(df_prp.wl$Tissue, levels = df_cs_count$Tissue)
p5 <- ggplot(df_prp.wl,aes(x=Tissue,y=value,fill=variable)) + geom_bar(stat = "identity",position = "stack") + theme_pubr() + 
  scale_fill_manual(breaks = c("PC_snv","PC_sv_tr","PC_sv_trLead"),values = c("#bdbdbd","#fc9272","#de2d26")) + 
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = .5)) + 
  xlab("Tissues") + ylab("# of credible sets");p2


cowplot::plot_grid(p4,p5,ncol = 1,align = "v")


# Fig.5d ------------------------------------------------------------------


df_comb <- readRDS("./input/Fig5d_data.RDS")

df_comb$cs_class <- factor(df_comb$cs_class,levels = c("SNV_TR","SNV_only","TR_only","SNV_SV_TR","SNV_SV","SV_only","SV_TR"))
df_comb.snv <- df_comb %>% filter(cs_class =="SNV_only")
df_comb.sv <- df_comb %>% filter(cs_class %in% c("SNV_SV_TR","SNV_SV","SV_only","SV_TR"))
df_comb.tr <- df_comb %>% filter(cs_class %in% c("SNV_SV_TR","SNV_TR","TR_only","SV_TR"))
df_comb.sv %<>% group_by(Var1) %>% mutate(Count=sum(Freq)) %>% select(Var1,Count) %>% distinct()
df_comb.sv$cs_class <- "cs_sv"
df_comb.tr %<>% group_by(Var1) %>% mutate(Count=sum(Freq)) %>% select(Var1,Count) %>% distinct()
df_comb.tr$cs_class <- "cs_tr"
names(df_comb.snv) <- c("Var1","Count","cs_class")
df_plot <- rbind(df_comb.snv,df_comb.sv,df_comb.tr)
df_plot$Var1 <- as.character(df_plot$Var1)
df_plot$cs_class <- as.character(df_plot$cs_class)
df_plot$Var1[df_plot$Var1=="SNV" & df_plot$cs_class=="SNV_only"] <- "not_lead"
df_plot$Var1[df_plot$Var1=="SV" & df_plot$cs_class=="cs_sv"] <- "lead"
df_plot$Var1[df_plot$Var1=="TR" & df_plot$cs_class=="cs_tr"] <- "lead"

df_plot$Var1[df_plot$Var1!="lead" & df_plot$cs_class=="cs_sv"] <- "not_lead"
df_plot$Var1[df_plot$Var1!="lead" & df_plot$cs_class=="cs_tr"] <- "not_lead"
df_plot$Var1 <- factor(df_plot$Var1,levels = c("not_lead","lead"))


p6 <- ggplot(df_plot,aes(x=cs_class,y=Count,fill=Var1)) + geom_bar(stat = "identity",position = "stack",width=.85) + theme_pubr() + 
  scale_fill_manual(breaks = c("not_lead","lead"),values = c("#bdbdbd","#FC9272")) +
  xlab("CS groups") + ylab("# of credible sets");p6


# Fig.5e ------------------------------------------------------------------


df_m <- readRDS("./input/SVTR_enrich_causal_CS.compare_to_SNV.RDS")
df_m$PIP_num <- df_m$PIP
df_m$PIP <- factor(df_m$PIP)

# load tissue color
df_tissue_color <- readRDS("./input/tissue_colors.RDS")

p7 <- ggplot(df_m,aes(x=PIP,y=OR)) + 
  geom_boxplot(outlier.shape = NA,aes(fill = PIP_num),alpha=.8) +
  geom_jitter(aes(color=Tissue),width=.2,size=1) + theme_pubr() +
  scale_color_manual(breaks = df_tissue_color$Tissue,values = paste0("#",df_tissue_color$Tissue_Color_Code)) +
  scale_fill_gradient(low = "#fff7ec",high="#b30000");p7


# conduct linear regression
res <- lm(OR ~ PIP_num, data=df_m)

summary(res) # Adjusted R2 = 0.706, P=6.477e-16
