setwd("F:/Ashley_Projects/Metabolism_Immune/13.PON1gene/05.RL/")
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

exp_pon1<-concat[, c("ID", "RL", "PON1_exp")]
exp_pon1$RL<-as.factor(exp_pon1$RL)
df<-exp_pon1
df_p_val1 <- df %>% 
  wilcox_test(PON1_exp  ~ RL) %>%
  adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  add_significance(p.col = "p.adj") %>% 
  add_xy_position(x = "Gene Cluster", dodge = 0.8) 
pdf("RLexp_boxplot.pdf")
p <- ggboxplot(df, x = "RL", y = "PON1_exp",fill = "RL", 
               color = 'RL', palette = c("#00AFBB", "#E7B800"), add = "boxplot", 
               add.params = list(fill = "white")) + 
  labs(x = "Gene cluster", y = "The expression of PON1", fill = "RL") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p + stat_compare_means(aes(group = RL), label = "p.signif", paired = FALSE)
dev.off()





