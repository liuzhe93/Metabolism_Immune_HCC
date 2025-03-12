setwd("F:/Ashley_Projects/Metabolism_Immune/03.ClusterTME/03_metabolic_activity/")
rm(list=ls())

mydata<-read.csv("../02_heatmap/SampleID_MetabolicScore.csv")
head(mydata)
mydata$sample_subtypes<-as.factor(mydata$sample_subtypes)

library("reshape2")
library("knitr")
mydata_transf<-melt(mydata, id.vars = c("X", "sample_subtypes"), variable.names = "Pathway", 
                    value.name = "MetabolicScore")
mydata_transf$variable<-gsub("scores_aa_df.", "", mydata_transf$variable)
head(mydata_transf)
library("rstatix")
myt_test<-t_test(group_by(mydata_transf,variable), MetabolicScore~sample_subtypes)
myt_test<-add_significance(myt_test, "p")
myt_test
## A tibble: 3 × 10
#  variable                  .y.            group1 group2    n1    n2  statistic  df        p    p.signif
#  <chr>                     <chr>          <chr>  <chr>  <int> <int>     <dbl>  <dbl>    <dbl>  <chr>   
#  1 amino_acid              MetabolicScore 1      2        126   264    -13.3    156.  1.62e-27 ****    
#  2 scores_lp_df.lipid      MetabolicScore 1      2        126   264    -18.2    178.  1.84e-42 ****    
#  3 scores_nt_df.nucleotide MetabolicScore 1      2        126   264     -4.14   187.  5.19e- 5 **** 
mywilcox_test<-wilcox_test(group_by(mydata_transf, variable), MetabolicScore~sample_subtypes)
mywilcox_test<-add_significance(mywilcox_test, "p")
mywilcox_test
## A tibble: 3 × 9
#  variable                  .y.            group1 group2    n1    n2 statistic        p p.signif
#  <chr>                     <chr>          <chr>  <chr>  <int> <int>     <dbl>    <dbl> <chr>   
#  1 amino_acid              MetabolicScore 1      2        126   264      3575    4.44e-36 ****    
#  2 scores_lp_df.lipid      MetabolicScore 1      2        126   264      1359    1.01e-48 ****    
#  3 scores_nt_df.nucleotide MetabolicScore 1      2        126   264     12890    3.26e- 4 ***

library("ggplot2")
library("ggpubr")
pdf("ThreeMetaPath_sig.pdf", height = 6, width = 8)
p <- ggboxplot(mydata_transf, x = "variable", y = "MetabolicScore",fill = "sample_subtypes", 
               color = 'sample_subtypes', palette = c("#00AFBB", "#E7B800"), add = "boxplot", 
               add.params = list(fill = "white")) + 
  labs(x = "Three key metabolic pathways", y = "Relative pathway activity", fill = "sample_subtypes") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p + stat_compare_means(aes(group = sample_subtypes), label = "p.signif", paired = FALSE)
dev.off()
