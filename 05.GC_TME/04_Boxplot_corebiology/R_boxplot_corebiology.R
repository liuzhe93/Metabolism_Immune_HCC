setwd("F:/Ashley_Projects/Metabolism_Immune/05.GC_TME/04_Boxplot_corebiology")
rm(list=ls())

mydata<-read.csv("F:/Ashley_Projects/Metabolism_Immune/01.IdentificationOfKeyMetabolicPathways/11_ssGSEA_corebiology/data_merged.csv")
head(mydata)
sample_cluster<-read.csv("F:/Ashley_Projects/Metabolism_Immune/04.MRGsClusters/04_MRG_clusters/Sample_Cluster.csv", header = T) 
head(sample_cluster)

merged_data<-cbind(sample_cluster, mydata)
dim(merged_data)
#[1] 390   24
head(merged_data)
merged_data<-merged_data[ , -which(colnames(merged_data) %in% c("X", "amino_acid", "lipid", "nucleotide"))]
dim(merged_data)
#[1] 390  20
head(merged_data)

library("reshape2")
library("knitr")
mydata_transf<-melt(merged_data, id.vars = c("sample_name", "cluster"), variable.names = "corebiology", 
                    value.name = "EnrichmentScore")
head(mydata_transf)
library("rstatix")
myt_test<-t_test(group_by(mydata_transf,variable), EnrichmentScore~cluster)
myt_test<-add_significance(myt_test, "p")
myt_test

mywilcox_test<-wilcox_test(group_by(mydata_transf, variable), EnrichmentScore~cluster)
mywilcox_test<-add_significance(mywilcox_test, "p")
mywilcox_test

library("ggplot2")
library("ggpubr")
mydata_transf$cluster<-as.factor(mydata_transf$cluster)

pdf("Corebiology_sig.pdf", height = 6, width = 8)
p <- ggboxplot(mydata_transf, x = "variable", y = "EnrichmentScore",fill = "cluster", 
               color = 'cluster', palette = c("blue", "red", "green"), add = "boxplot", 
               add.params = list(fill = "white")) + 
  labs(x = "18 core biological pathway", y = "Relative enrichment score", fill = "cluster") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p + stat_compare_means(aes(group = cluster), label = "p.signif", paired = FALSE)
dev.off()



