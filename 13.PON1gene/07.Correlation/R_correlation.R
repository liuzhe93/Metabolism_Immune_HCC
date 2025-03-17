setwd("F:/Ashley_Projects/Metabolism_Immune/13.PON1gene/07.Correlation/")
rm(list=ls())

merged_data_sel<-read.csv("F:/Ashley_Projects/Metabolism_Immune/01.IdentificationOfKeyMetabolicPathways/11_ssGSEA_corebiology/data_merged.csv", header = T)
dim(merged_data_sel)
#[1] 390  22

dataExp<-read.csv("F:/Ashley_Projects/Metabolism_Immune/11.Immunotherapy/02.Drug/geneExp.csv", header = T, row.names = 1)
dim(dataExp)
#[1] 12925   390
dataExp_t<-t(dataExp)
dim(dataExp_t)
#[1]   390 12925
dataExp_t<-as.data.frame(dataExp_t)
dataExp_t$X<-row.names(dataExp_t)
dataExp_t_sel_PON1<-dataExp_t[, c("X", "PON1")]


data_all<-cbind(merged_data_sel[,5:ncol(merged_data_sel)], dataExp_t_sel_PON1$PON1)
dim(data_all)
#[1] 390  22
rownames(data_all)<-merged_data_sel$X
colnames(data_all)[ncol(data_all)]<-"PON1"

data_all_format<-cbind(data_all[, c(19, 1:18)])

M = cor(data_all_format, use = "everything", method = "spearman")
library(corrplot)
testRes = cor.mtest(data_all_format,  use = "everything", method = "spearman", conf.level = 0.95)

## add significant level stars
pdf("Cor_18corebio.pdf", height = 8, width = 8)
corrplot(M, p.mat = testRes$p, method = 'circle', type = 'lower', diag = T, insig='blank', pch.col = 'grey20',
         order = 'original')
dev.off()




