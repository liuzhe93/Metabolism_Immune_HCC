setwd("F:/Ashley_Projects/Metabolism_Immune/05.GC_TME/03_Boxplot_28immunecells")
rm(list=ls())

mydata<-read.csv("F:/Ashley_Projects/Metabolism_Immune/05.GC_TME/02_Heatmap28immunecells/Cluster_immunecells.csv", row.names = 1)

mydata_sel<-mydata[c(1,10:37),]
mydata_t<-t(mydata_sel)
mydata_t<-as.data.frame(mydata_t)
mydata_t$cluster<-as.factor(mydata_t$cluster)
mydata_t$X<-row.names(mydata_t)

library("reshape2")
library("knitr")
mydata_transf<-melt(mydata_t, id.vars = c("X", "cluster"), variable.names = "ImmuneCells", 
                    value.name = "MetabolicScore")
head(mydata_transf)
library("rstatix")
myt_test<-t_test(group_by(mydata_transf,variable), MetabolicScore~cluster)
myt_test<-add_significance(myt_test, "p")
myt_test

mywilcox_test<-wilcox_test(group_by(mydata_transf, variable), MetabolicScore~cluster)
mywilcox_test<-add_significance(mywilcox_test, "p")
mywilcox_test

library("ggplot2")
library("ggpubr")
pdf("Immunecells_sig.pdf", height = 6, width = 8)
p <- ggboxplot(mydata_transf, x = "variable", y = "MetabolicScore",fill = "cluster", 
               color = 'cluster', palette = c("blue", "red", "green"), add = "boxplot", 
               add.params = list(fill = "white")) + 
  labs(x = "28 immune cell infiltration", y = "Relative enrichment score", fill = "sample_subtypes") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p + stat_compare_means(aes(group = cluster), label = "p.signif", paired = FALSE)
dev.off()



