setwd("F:/Ashley_Projects/Metabolism_Immune/04.MRGsClusters/05_Heatmap/")
rm(list=ls())

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

sample_cluster<-read.csv("F:/Ashley_Projects/Metabolism_Immune/04.MRGsClusters/04_MRG_clusters/Sample_Cluster.csv", header = T) 
dim(sample_cluster)
#[1] 390   2
head(sample_cluster)

merged<-cbind(sample_cluster, scores_aa_df$amino_acid, scores_lp_df$lipid, scores_nt_df$nucleotide)
rownames(merged)<-merged$sample_name
merged_new<-merged[,-2]
head(merged_new)
write.csv(merged_new, "SampleID_MetabolicScore.csv", quote = F)
merged_srt<-merged_new[order(merged_new[,1]),]
dim(merged_srt)
#[1] 390   4
head(merged_srt)
table(merged_srt$cluster)
# 1   2   3 
#88 193 109

merged_srt<-merged_srt[,-1]
head(merged_srt)

colnames(merged_srt)<-c("amino_acid","lipid","nucleotide")
merged_t<-t(merged_srt)
dim(merged_t)
#[1]   3 390
merged_t[1:3,1:3]
mat<-merged_t

library("pheatmap")
packageVersion("pheatmap")
#[1] ‘1.0.12’
library("ggplot2")
library("ggplotify")

annotation_col = data.frame(Cluster = c(rep("1",88),rep("2",193),rep("3",109)))
rownames(annotation_col) = colnames(mat)
ann_colors = list(Cluster = c("1" = "blue", "2" = "red", "3"="green"))

heatmap_input <-mat

pdf(file="heatmap_C123.pdf",width = 6,height = 3)
pheatmap(heatmap_input,scale = "row", clustering_distance_rows = "correlation", cluster_rows = F, 
         cluster_col = F, color = colorRampPalette(c("blue", "white", "red"))(50), 
         annotation_col = annotation_col, gaps_col = c(88, 281),angle_col = "45", 
         annotation_colors = ann_colors, show_rownames = T, show_colnames = F, main = "Title")
#fontsize = 1
dev.off()

