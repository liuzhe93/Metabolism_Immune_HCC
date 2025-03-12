setwd("F:/Ashley_Projects/Metabolism_Immune/01.IdentificationOfKeyMetabolicPathways/09_ssGSEA_immunecells")
rm(list=ls())

dataExp<-read.table("F:/Ashley_Projects/Metabolism_Immune/01.IdentificationOfKeyMetabolicPathways/08_EstimateScores/dataexp.txt", header = T, row.names = 1)

genelist<-read.csv("ImmunecellMetagenes.csv")
#28 immune cell types and metagenes are collected from https://doi.org/10.1016/j.celrep.2016.12.019
gene_set<-genelist[,1:2]
head(gene_set)
#  Metagene        Cell.type
#1   ADAM28 Activated B cell
#2    CD180 Activated B cell
#3    CD79B Activated B cell
#4      BLK Activated B cell
#5     CD19 Activated B cell
#6    MS4A1 Activated B cell
list<-split(as.matrix(gene_set)[,1], gene_set[,2])
list

library("GSVA")
ssgsea_par <- ssgseaParam(as.matrix(dataExp), list)
gsva_matrix <- gsva(ssgsea_par)
dim(gsva_matrix)
#[1]  28 390
gsva_matrix<-as.data.frame(gsva_matrix)
gsva_matrix_t<-t(gsva_matrix)
dim(gsva_matrix_t)
#[1] 390  28

scores_aa<-read.csv("F:/Ashley_Projects/Metabolism_Immune/01.IdentificationOfKeyMetabolicPathways/01_Aminoacid/Normalize.csv", row.names = 1)
scores_aa_df<-as.data.frame(t(scores_aa))
colnames(scores_aa_df)<-"amino_acid"
head(scores_aa_df)

scores_lp<-read.csv("F:/Ashley_Projects/Metabolism_Immune/01.IdentificationOfKeyMetabolicPathways/04_Lipid/Normalize.csv", row.names = 1)
scores_lp_df<-as.data.frame(t(scores_lp))
colnames(scores_lp_df)<-"lipid"
head(scores_lp_df)

scores_nt<-read.csv("F:/Ashley_Projects/Metabolism_Immune/01.IdentificationOfKeyMetabolicPathways/05_Nucleotide/Normalize.csv", row.names = 1)
scores_nt_df<-as.data.frame(t(scores_nt))
colnames(scores_nt_df)<-"nucleotide"
head(scores_nt_df)

dim(scores_aa_df)
#[1] 390   1
dim(scores_lp_df)
#[1] 390   1
dim(scores_nt_df)
#[1] 390   1

merged_scores_met<-cbind(scores_aa_df, scores_lp_df, scores_nt_df)
merged_scores_met$sample<-row.names(merged_scores_met)
head(merged_scores_met)

data_merged<-cbind(merged_scores_met, gsva_matrix_t)
data_merged_sel<-data_merged[,-4]
write.csv(data_merged_sel, "data_merged.csv", quote = F)

M = cor(data_merged_sel, use = "everything", method = "spearman")
library(corrplot)
testRes = cor.mtest(data_merged_sel,  use = "everything", method = "spearman", conf.level = 0.95)

## add significant level stars
pdf("Cor_28immune.pdf", height = 8, width = 8)
corrplot(M, p.mat = testRes$p, method = 'color', diag = T, type = 'lower',
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,
         insig = 'label_sig', pch.col = 'grey20', order = 'original')
dev.off()

#0.05-->*
#0.01-->**
#0.001-->***
