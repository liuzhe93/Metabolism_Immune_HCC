setwd("F:/Ashley_Projects/Metabolism_Immune/05.GC_TME/01_EstimateScores")
rm(list=ls())

library("estimate")

mydata<-read.csv("F:/Ashley_Projects/Metabolism_Immune/01.IdentificationOfKeyMetabolicPathways/08_EstimateScores/scores_estimate.csv")
head(mydata)
sample_cluster<-read.csv("F:/Ashley_Projects/Metabolism_Immune/04.MRGsClusters/05_Heatmap/SampleID_MetabolicScore.csv", header = T) 
head(sample_cluster)

merged_data<-cbind(sample_cluster, mydata)
dim(merged_data)
#[1] 390   11
head(merged_data)


library("tidyverse")
library("gapminder")
library("ggsci")
library("ggprism")
library("rstatix")
library("ggpubr")
library("reshape2")
library("knitr")
library("rstatix")


###############################################################################################
StromalScore<-merged_data[, c("X", "cluster", "StromalScore")]
StromalScore$cluster<-as.factor(StromalScore$cluster)
df<-StromalScore
df_p_val1 <- df %>% 
  wilcox_test(StromalScore  ~ cluster) %>%
  adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  add_significance(p.col = "p.adj") %>% 
  add_xy_position(x = "Metabolism Cluster", dodge = 0.8) 
pdf("StromalScore_boxplot.pdf")
p <- ggboxplot(df, x = "cluster", y = "StromalScore",fill = "cluster", 
               color = 'cluster', palette = c("blue", "red","green"), add = "boxplot", 
               add.params = list(fill = "white")) + 
  labs(x = "Metabolism cluster", y = "Stromal Score", fill = "cluster") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p + stat_compare_means(aes(group = cluster), label = "p.signif", paired = FALSE)
dev.off()

###############################################################################################
ImmuneScore<-merged_data[, c("X", "cluster", "ImmuneScore")]
ImmuneScore$cluster<-as.factor(ImmuneScore$cluster)
df<-ImmuneScore
df_p_val1 <- df %>% 
  wilcox_test(ImmuneScore  ~ cluster) %>%
  adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  add_significance(p.col = "p.adj") %>% 
  add_xy_position(x = "Metabolism Cluster", dodge = 0.8) 
pdf("ImmuneScore_boxplot.pdf")
p <- ggboxplot(df, x = "cluster", y = "ImmuneScore",fill = "cluster", 
               color = 'cluster', palette = c("blue", "red","green"), add = "boxplot", 
               add.params = list(fill = "white")) + 
  labs(x = "Metabolism cluster", y = "Immune Score", fill = "cluster") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p + stat_compare_means(aes(group = cluster), label = "p.signif", paired = FALSE)
dev.off()

###############################################################################################
TumourPurity<-merged_data[, c("X", "cluster", "TumourPurity")]
TumourPurity$cluster<-as.factor(TumourPurity$cluster)
df<-TumourPurity
df_p_val1 <- df %>% 
  wilcox_test(TumourPurity  ~ cluster) %>%
  adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  add_significance(p.col = "p.adj") %>% 
  add_xy_position(x = "Metabolism Cluster", dodge = 0.8) 
pdf("TumourPurity_boxplot.pdf")
p <- ggboxplot(df, x = "cluster", y = "TumourPurity",fill = "cluster", 
               color = 'cluster', palette = c("blue", "red","green"), add = "boxplot", 
               add.params = list(fill = "white")) + 
  labs(x = "Metabolism cluster", y = "Tumour Purity", fill = "cluster") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p + stat_compare_means(aes(group = cluster), label = "p.signif", paired = FALSE)
dev.off()







