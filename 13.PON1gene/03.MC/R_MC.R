setwd("F:/Ashley_Projects/Metabolism_Immune/13.PON1gene/03.MC/")
rm(list=ls())

dataExp<-read.csv("F:/Ashley_Projects/Metabolism_Immune/11.Immunotherapy/02.Drug/geneExp.csv", header = T, row.names = 1)
dim(dataExp)
#[1] 12925   390

dataExp_t<-t(dataExp)
dim(dataExp_t)
#[1]   390 12925
dataExp_t<-as.data.frame(dataExp_t)
dataExp_t$Id<-row.names(dataExp_t)
dataExp_t$Id<-substring(dataExp_t$Id, 1, 12)
dataExp_t$Id<-gsub("\\.", "-", dataExp_t$Id)

category_info<-read.csv("F:/Ashley_Projects/Metabolism_Immune/10.AlluvialPlot/AP_data.csv", header = T)
dim(category_info)
#[1] 577   6
head(category_info)

merged_data<-merge(category_info, dataExp_t, by.x = "ID", by.y = "Id")
merged_data_sel<-subset(merged_data, select = c("ID", "MC", "GC.x", "RL", "PON1"))
colnames(merged_data_sel)[3]<-"GC"
head(merged_data_sel)

geneExp<-read.csv("F:/Ashley_Projects/Metabolism_Immune/13.PON1gene/01.internal/res.csv", header = T)
dim(geneExp)
#[1] 355   5
head(geneExp)

concat<-merge(merged_data_sel, geneExp, by.x = "ID", by.y = "Id")
colnames(concat)[5]<-"PON1_exp"
colnames(concat)[9]<-"PON1_class"
concat$X<-NULL

concat$MC<-ifelse(concat$MC=="C1", "MC1", "MC2")
concat$GC<-ifelse(concat$GC=="G1", "GC1", ifelse(concat$GC=="G2", "GC2", "GC3"))
concat$PON1_class<-ifelse(concat$PON1_class=="high", "High", "Low")
head(concat)

write.csv(concat, "geneExp_class.csv", quote = F)

library("tidyverse")
library("gapminder")
library("ggsci")
library("ggprism")
library("rstatix")
library("ggpubr")
library("reshape2")
library("knitr")
library("rstatix")

exp_pon1<-concat[, c("ID", "MC", "PON1_exp")]
exp_pon1$MC<-as.factor(exp_pon1$MC)
df<-exp_pon1
df_p_val1 <- df %>% 
  wilcox_test(PON1_exp  ~ MC) %>%
  adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  add_significance(p.col = "p.adj") %>% 
  add_xy_position(x = "Metabolism Cluster", dodge = 0.8) 
pdf("MCexp_boxplot.pdf")
p <- ggboxplot(df, x = "MC", y = "PON1_exp",fill = "MC", 
               color = 'MC', palette = c("#00AFBB", "#E7B800"), add = "boxplot", 
               add.params = list(fill = "white")) + 
  labs(x = "Metabolism cluster", y = "The expression of PON1", fill = "MC") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p + stat_compare_means(aes(group = MC), label = "p.signif", paired = FALSE)
dev.off()


