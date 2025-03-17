setwd("F:/Ashley_Projects/Metabolism_Immune/03.ClusterTME/09_boxplot_corebiology")
rm(list=ls())


mydata<-read.csv("F:/Ashley_Projects/Metabolism_Immune/01.IdentificationOfKeyMetabolicPathways/11_ssGSEA_corebiology/data_merged.csv")
head(mydata)
sample_cluster<-read.csv("F:/Ashley_Projects/Metabolism_Immune/03.ClusterTME/01_cluster/Sample_Cluster.csv", header = T) 
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
mydata_transf<-melt(merged_data, id.vars = c("ID", "sample_subtypes"), variable.names = "corebiology", 
                    value.name = "EnrichmentScore")
head(mydata_transf)
library("rstatix")
myt_test<-t_test(group_by(mydata_transf,variable), EnrichmentScore~sample_subtypes)
myt_test<-add_significance(myt_test, "p")
myt_test
# A tibble: 18 × 10
#variable                      .y.             group1 group2    n1    n2 statistic    df        p p.signif
#<fct>                         <chr>           <chr>  <chr>  <int> <int>     <dbl> <dbl>    <dbl> <chr>   
#1 Angiogenesis                  EnrichmentScore 1      2        126   264    -7.11   239. 1.34e-11 ****    
#2 Antigen.processing.machinery  EnrichmentScore 1      2        126   264    -0.151  201. 8.8 e- 1 ns      
#3 CD8.T.effector                EnrichmentScore 1      2        126   264    -0.206  190. 8.37e- 1 ns      
#4 Cell.cycle                    EnrichmentScore 1      2        126   264    12.6    221. 1.16e-27 ****    
#5 Cell.cycle.regulators         EnrichmentScore 1      2        126   264     0.903  180. 3.68e- 1 ns      
#6 DNA.damage.repair             EnrichmentScore 1      2        126   264     8.08   260. 2.47e-14 ****    
#7 DNA.replication               EnrichmentScore 1      2        126   264    10.2    285. 6.27e-21 ****    
#8 EMT_1                         EnrichmentScore 1      2        126   264    -0.667  200. 5.06e- 1 ns      
#9 EMT_2                         EnrichmentScore 1      2        126   264    -2.03   210. 4.4 e- 2 *       
#10 EMT_3                        EnrichmentScore 1      2        126   264     3.01   267. 2.86e- 3 **      
#11 Fanconi.anemia               EnrichmentScore 1      2        126   264     8.59   267. 7.41e-16 ****    
#12 FGFR3.related.genes          EnrichmentScore 1      2        126   264     7.70   191. 7.15e-13 ****    
#13 Homologous.recombination     EnrichmentScore 1      2        126   264     6.50   211. 5.88e-10 ****    
#14 Immune.checkpoint            EnrichmentScore 1      2        126   264     4.46   158. 1.58e- 5 ****    
#15 Mismatch.repair              EnrichmentScore 1      2        126   264     8.04   270. 2.84e-14 ****    
#16 Nucleotide.excision.repair   EnrichmentScore 1      2        126   264     8.45   263. 2.04e-15 ****    
#17 Pan.F.TBRS                   EnrichmentScore 1      2        126   264     0.230  198. 8.18e- 1 ns      
#18 WNT.target                   EnrichmentScore 1      2        126   264     1.92   218. 5.63e- 2 ns 
mywilcox_test<-wilcox_test(group_by(mydata_transf, variable), EnrichmentScore~sample_subtypes)
mywilcox_test<-add_significance(mywilcox_test, "p")
mywilcox_test
# A tibble: 18 × 9
#variable                        .y.             group1 group2    n1    n2 statistic        p p.signif
#<fct>                           <chr>           <chr>  <chr>  <int> <int>     <dbl>    <dbl> <chr>   
#1  Angiogenesis                 EnrichmentScore 1      2        126   264     9334. 2.38e-12 ****    
#2  Antigen.processing.machinery EnrichmentScore 1      2        126   264    16654  9.84e- 1 ns      
#3  CD8.T.effector               EnrichmentScore 1      2        126   264    15030. 1.24e- 1 ns      
#4  Cell.cycle                   EnrichmentScore 1      2        126   264    27782  9.18e-27 ****    
#5  Cell.cycle.regulators        EnrichmentScore 1      2        126   264    17776  2.72e- 1 ns      
#6  DNA.damage.repair            EnrichmentScore 1      2        126   264    24205  3.5 e-13 ****    
#7  DNA.replication              EnrichmentScore 1      2        126   264    25839  9.30e-19 ****    
#8  EMT_1                        EnrichmentScore 1      2        126   264    15360  2.22e- 1 ns      
#9  EMT_2                        EnrichmentScore 1      2        126   264    14530  4.35e- 2 *       
#10 EMT_3                        EnrichmentScore 1      2        126   264    19600. 4.36e- 3 **      
#11 Fanconi.anemia               EnrichmentScore 1      2        126   264    24590  2.11e-14 ****    
#12 FGFR3.related.genes          EnrichmentScore 1      2        126   264    25277  1.01e-16 ****    
#13 Homologous.recombination     EnrichmentScore 1      2        126   264    23135  4.21e-10 ****    
#14 Immune.checkpoint            EnrichmentScore 1      2        126   264    20877  4.56e- 5 ****    
#15 Mismatch.repair              EnrichmentScore 1      2        126   264    24087  8.05e-13 ****    
#16 Nucleotide.excision.repair   EnrichmentScore 1      2        126   264    24548  2.89e-14 ****    
#17 Pan.F.TBRS                   EnrichmentScore 1      2        126   264    16986  7.34e- 1 ns      
#18 WNT.target                   EnrichmentScore 1      2        126   264    19094  1.81e- 2 *
  
  
library("ggplot2")
library("ggpubr")
mydata_transf$sample_subtypes<-as.factor(mydata_transf$sample_subtypes)

pdf("Corebiology_sig.pdf", height = 6, width = 8)
p <- ggboxplot(mydata_transf, x = "variable", y = "EnrichmentScore",fill = "sample_subtypes", 
               color = 'sample_subtypes', palette = c("#00AFBB", "#E7B800"), add = "boxplot", 
               add.params = list(fill = "white")) + 
  labs(x = "18 core biological pathway", y = "Relative enrichment score", fill = "sample_subtypes") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p + stat_compare_means(aes(group = sample_subtypes), label = "p.signif", paired = FALSE)
dev.off()



