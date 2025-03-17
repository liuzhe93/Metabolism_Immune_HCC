setwd("F:/Ashley_Projects/Metabolism_Immune/13.PON1gene/04.GC/")
rm(list=ls())

concat<-read.csv("F:/Ashley_Projects/Metabolism_Immune/13_PON1gene/03.MC/geneExp_class.csv", header = T, row.names = 1)
head(concat)

library("tidyverse")
library("gapminder")
library("ggsci")
library("ggprism")
library("rstatix")
library("ggpubr")
library("reshape2")
library("knitr")
library("rstatix")

exp_pon1<-concat[, c("ID", "GC", "PON1_exp")]
exp_pon1$GC<-as.factor(exp_pon1$GC)
df<-exp_pon1
df_p_val1 <- df %>% 
  wilcox_test(PON1_exp  ~ GC) %>%
  adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  add_significance(p.col = "p.adj") %>% 
  add_xy_position(x = "Gene Cluster", dodge = 0.8) 
pdf("GCexp_boxplot.pdf")
p <- ggboxplot(df, x = "GC", y = "PON1_exp",fill = "GC", 
               color = 'GC', palette = c("blue", "red","green"), add = "boxplot", 
               add.params = list(fill = "white")) + 
  labs(x = "Gene cluster", y = "The expression of PON1", fill = "GC") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p + stat_compare_means(aes(group = GC), label = "p.signif", paired = FALSE)
dev.off()





