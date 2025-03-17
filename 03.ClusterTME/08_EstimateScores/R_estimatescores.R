setwd("F:/Ashley_Projects/Metabolism_Immune/03.ClusterTME/08_EstimateScores")
rm(list=ls())

library("estimate")

mydata<-read.csv("F:/Ashley_Projects/Metabolism_Immune/01.IdentificationOfKeyMetabolicPathways/08_EstimateScores/scores_estimate.csv")
head(mydata)
sample_cluster<-read.csv("F:/Ashley_Projects/Metabolism_Immune/03.ClusterTME/01_cluster/Sample_Cluster.csv", header = T) 
head(sample_cluster)

merged_data<-cbind(sample_cluster, mydata)
dim(merged_data)
#[1] 390   8
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
StromalScore<-merged_data[, c("ID", "sample_subtypes", "StromalScore")]
StromalScore$sample_subtypes<-as.factor(StromalScore$sample_subtypes)
df<-StromalScore
df_p_val1 <- df %>% 
  wilcox_test(StromalScore  ~ sample_subtypes) %>%
  adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  add_significance(p.col = "p.adj") %>% 
  add_xy_position(x = "Metabolism Cluster", dodge = 0.8) 
pdf("StromalScore_boxplot.pdf")
p <- ggboxplot(df, x = "sample_subtypes", y = "StromalScore",fill = "sample_subtypes", 
               color = 'sample_subtypes', palette = c("#00AFBB", "#E7B800"), add = "boxplot", 
               add.params = list(fill = "white")) + 
  labs(x = "Metabolism cluster", y = "Stromal Score", fill = "sample_subtypes") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p + stat_compare_means(aes(group = sample_subtypes), label = "p.signif", paired = FALSE)
dev.off()

###############################################################################################
ImmuneScore<-merged_data[, c("ID", "sample_subtypes", "ImmuneScore")]
ImmuneScore$sample_subtypes<-as.factor(ImmuneScore$sample_subtypes)
df<-ImmuneScore
df_p_val1 <- df %>% 
  wilcox_test(ImmuneScore  ~ sample_subtypes) %>%
  adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  add_significance(p.col = "p.adj") %>% 
  add_xy_position(x = "Metabolism Cluster", dodge = 0.8) 
pdf("ImmuneScore_boxplot.pdf")
p <- ggboxplot(df, x = "sample_subtypes", y = "ImmuneScore",fill = "sample_subtypes", 
               color = 'sample_subtypes', palette = c("#00AFBB", "#E7B800"), add = "boxplot", 
               add.params = list(fill = "white")) + 
  labs(x = "Metabolism cluster", y = "Immune Score", fill = "sample_subtypes") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p + stat_compare_means(aes(group = sample_subtypes), label = "p.signif", paired = FALSE)
dev.off()

###############################################################################################
TumourPurity<-merged_data[, c("ID", "sample_subtypes", "TumourPurity")]
TumourPurity$sample_subtypes<-as.factor(TumourPurity$sample_subtypes)
df<-TumourPurity
df_p_val1 <- df %>% 
  wilcox_test(TumourPurity  ~ sample_subtypes) %>%
  adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  add_significance(p.col = "p.adj") %>% 
  add_xy_position(x = "Metabolism Cluster", dodge = 0.8) 
pdf("TumourPurity_boxplot.pdf")
p <- ggboxplot(df, x = "sample_subtypes", y = "TumourPurity",fill = "sample_subtypes", 
               color = 'sample_subtypes', palette = c("#00AFBB", "#E7B800"), add = "boxplot", 
               add.params = list(fill = "white")) + 
  labs(x = "Metabolism cluster", y = "Tumour Purity", fill = "sample_subtypes") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p + stat_compare_means(aes(group = sample_subtypes), label = "p.signif", paired = FALSE)
dev.off()




