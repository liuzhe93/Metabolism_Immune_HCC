setwd("F:/Ashley_Projects/Metabolism_Immune/04.MRGsClusters/06_Metabolic_activity/")
rm(list=ls())

mydata<-read.csv("F:/Ashley_Projects/Metabolism_Immune/04.MRGsClusters/05_Heatmap/SampleID_MetabolicScore.csv")
head(mydata)
mydata$cluster<-as.factor(mydata$cluster)

library("reshape2")
library("knitr")
mydata_transf<-melt(mydata, id.vars = c("X", "cluster"), variable.names = "Pathway", 
                    value.name = "MetabolicScore")
mydata_transf$variable<-gsub("scores_aa_df.", "", mydata_transf$variable)
mydata_transf$variable<-gsub("scores_lp_df.", "", mydata_transf$variable)
mydata_transf$variable<-gsub("scores_nt_df.", "", mydata_transf$variable)

head(mydata_transf)
library("rstatix")
myt_test<-t_test(group_by(mydata_transf,variable), MetabolicScore~cluster)
myt_test<-add_significance(myt_test, "p")
myt_test
## A tibble: 9 × 12
#  variable   .y.          group1 group2    n1    n2 statistic    df        p    p.adj p.adj.signif p.signif
#  <chr>      <chr>        <chr>  <chr>  <int> <int>     <dbl> <dbl>    <dbl>    <dbl> <chr>        <chr>   
#  1 amino_acid MetabolicSc… 1      2         88   193    -14.0   101. 2.51e-25 7.53e-25 ****         ****    
#  2 amino_acid MetabolicSc… 1      3         88   109     -8.51  123. 4.88e-14 4.88e-14 ****         ****    
#  3 amino_acid MetabolicSc… 2      3        193   109      9.55  189. 7.04e-18 1.41e-17 ****         ****    
#  4 lipid      MetabolicSc… 1      2         88   193    -20.1   117. 1.22e-39 3.66e-39 ****         ****    
#  5 lipid      MetabolicSc… 1      3         88   109    -11.8   126. 3.72e-22 3.72e-22 ****         ****    
#  6 lipid      MetabolicSc… 2      3        193   109     13.6   247. 1.15e-31 2.3 e-31 ****         ****    
#  7 nucleotide MetabolicSc… 1      2         88   193     -4.67  126. 7.76e- 6 2.33e- 5 ****         ****    
#  8 nucleotide MetabolicSc… 1      3         88   109     -3.34  172. 1   e- 3 2   e- 3 **           ***     
#  9 nucleotide MetabolicSc… 2      3        193   109      1.04  187. 2.99e- 1 2.99e- 1 ns           ns 

mywilcox_test<-wilcox_test(group_by(mydata_transf, variable), MetabolicScore~cluster)
mywilcox_test<-add_significance(mywilcox_test, "p")
mywilcox_test
## A tibble: 9 × 11
# variable   .y.            group1 group2    n1    n2 statistic        p    p.adj p.adj.signif p.signif
# <chr>      <chr>          <chr>  <chr>  <int> <int>     <dbl>    <dbl>    <dbl> <chr>        <chr>   
#  1 amino_acid MetabolicScore 1      2         88   193       902 3.03e-33 9.09e-33 ****         ****    
#  2 amino_acid MetabolicScore 1      3         88   109      1760 2.34e-14 2.34e-14 ****         ****    
#  3 amino_acid MetabolicScore 2      3        193   109     16739 1.41e-17 2.82e-17 ****         ****    
#  4 lipid      MetabolicScore 1      2         88   193       102 3.04e-40 9.12e-40 ****         ****    
#  5 lipid      MetabolicScore 1      3         88   109       773 4.93e-24 4.93e-24 ****         ****    
#  6 lipid      MetabolicScore 2      3        193   109     18353 6   e-27 1.20e-26 ****         ****    
#  7 nucleotide MetabolicScore 1      2         88   193      5880 3.57e- 5 1.07e- 4 ***          ****    
#  8 nucleotide MetabolicScore 1      3         88   109      3658 4   e- 3 9   e- 3 **           **      
#  9 nucleotide MetabolicScore 2      3        193   109     11210 3.43e- 1 3.43e- 1 ns           ns
  
library("ggplot2")
library("ggpubr")
pdf("ThreeMetaPath_sig.pdf", height = 6, width = 8)
p <- ggboxplot(mydata_transf, x = "variable", y = "MetabolicScore",fill = "cluster", 
               color = 'cluster', palette = c("1" = "blue", "2" = "red", "3"="green"), add = "boxplot", 
               add.params = list(fill = "white")) + 
  labs(x = "Three key metabolic pathways", y = "Relative pathway activity", fill = "cluster") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p + stat_compare_means(aes(group = cluster), label = "p.signif", paired = FALSE)
dev.off()





  
