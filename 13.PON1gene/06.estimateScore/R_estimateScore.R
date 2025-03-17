setwd("F:/Ashley_Projects/Metabolism_Immune/13.PON1gene/06.estimateScore/")
rm(list=ls())

concat<-read.csv("F:/Ashley_Projects/Metabolism_Immune/13.PON1gene/03.MC/geneExp_class.csv", header = T, row.names = 1)
head(concat)

library("estimate")

mydata<-read.csv("F:/Ashley_Projects/Metabolism_Immune/01.IdentificationOfKeyMetabolicPathways/08_EstimateScores/scores_estimate.csv")
head(mydata)

merged_data<-merge(concat, mydata, by.x = "ID", by.y = "sample_name")

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
StromalScore<-merged_data[, c("ID", "PON1_class", "StromalScore")]
StromalScore$PON1_class<-as.factor(StromalScore$PON1_class)
df<-StromalScore
df_p_val1 <- df %>% 
  wilcox_test(StromalScore  ~ PON1_class) %>%
  adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  add_significance(p.col = "p.adj") %>% 
  add_xy_position(x = "Metabolism Cluster", dodge = 0.8) 
pdf("StromalScore_boxplot.pdf")
p <- ggboxplot(df, x = "PON1_class", y = "StromalScore",fill = "PON1_class", 
               color = 'PON1_class', palette = c("#00AFBB", "#E7B800"), add = "boxplot", 
               add.params = list(fill = "white")) + 
  labs(x = "Metabolism cluster", y = "Stromal Score", fill = "PON1_class") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p + stat_compare_means(aes(group = PON1_class), label = "p.signif", paired = FALSE)
dev.off()

###############################################################################################
ImmuneScore<-merged_data[, c("ID", "PON1_class", "ImmuneScore")]
ImmuneScore$PON1_class<-as.factor(ImmuneScore$PON1_class)
df<-ImmuneScore
df_p_val1 <- df %>% 
  wilcox_test(ImmuneScore  ~ PON1_class) %>%
  adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  add_significance(p.col = "p.adj") %>% 
  add_xy_position(x = "Metabolism Cluster", dodge = 0.8) 
pdf("ImmuneScore_boxplot.pdf")
p <- ggboxplot(df, x = "PON1_class", y = "ImmuneScore",fill = "PON1_class", 
               color = 'PON1_class', palette = c("#00AFBB", "#E7B800"), add = "boxplot", 
               add.params = list(fill = "white")) + 
  labs(x = "Metabolism cluster", y = "Immune Score", fill = "PON1_class") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p + stat_compare_means(aes(group = PON1_class), label = "p.signif", paired = FALSE)
dev.off()

###############################################################################################
TumourPurity<-merged_data[, c("ID", "PON1_class", "TumourPurity")]
TumourPurity$PON1_class<-as.factor(TumourPurity$PON1_class)
df<-TumourPurity
df_p_val1 <- df %>% 
  wilcox_test(TumourPurity  ~ PON1_class) %>%
  adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  add_significance(p.col = "p.adj") %>% 
  add_xy_position(x = "Metabolism Cluster", dodge = 0.8) 
pdf("TumourPurity_boxplot.pdf")
p <- ggboxplot(df, x = "PON1_class", y = "TumourPurity",fill = "PON1_class", 
               color = 'PON1_class', palette = c("#00AFBB", "#E7B800"), add = "boxplot", 
               add.params = list(fill = "white")) + 
  labs(x = "Metabolism cluster", y = "Tumour Purity", fill = "PON1_class") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p + stat_compare_means(aes(group = PON1_class), label = "p.signif", paired = FALSE)
dev.off()




