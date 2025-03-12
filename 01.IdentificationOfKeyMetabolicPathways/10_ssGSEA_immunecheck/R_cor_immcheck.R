setwd("F:/Ashley_Projects/Metabolism_Immune/01.IdentificationOfKeyMetabolicPathways/10_ssGSEA_immunecheck")
rm(list=ls())

#https://www.nature.com/articles/nature25501
#Supplementary Table 8

dataExp<-read.table("F:/Ashley_Projects/Metabolism_Immune/01.IdentificationOfKeyMetabolicPathways/08_EstimateScores/dataexp.txt", header = T, row.names = 1)
icpg<-c("IDO1", "CD274", "HAVCR2", "PDCD1", "CTLA4", "LAG3", "PDCD1LG2")
dataExp_sel<-dataExp[icpg,]
dim(dataExp_sel)
#[1]   7 390
dataExp_sel_t<-t(dataExp_sel)

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

merged_scores_met<-cbind(scores_aa_df, scores_lp_df, scores_nt_df)
head(merged_scores_met)

data_merged<-cbind(merged_scores_met, dataExp_sel_t)
write.csv(data_merged, "data_merged.csv", quote = F)


library(corrplot)
M = cor(data_merged, use = "everything", method = "spearman")
testRes = cor.mtest(data_merged,  use = "everything", method = "spearman", conf.level = 0.95)

pdf("Cor_ICPGs.pdf")
corrplot(M, p.mat = testRes$p, method = 'circle', type = 'full', insig='blank',
         addCoef.col ='black', number.cex = 0.8, order = 'original', diag=T)
dev.off()
#0.05-->*
#0.01-->**
#0.001-->***
write.csv(M, "cor.csv", quote = F)
write.csv(testRes$p, "pvalues.csv", quote = F)


