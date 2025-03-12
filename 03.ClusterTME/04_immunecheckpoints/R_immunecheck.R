setwd("F:/Ashley_Projects/Metabolism_Immune/03.ClusterTME/04_immunecheckpoints")
rm(list=ls())

icpg<-c("IDO1", "CD274", "HAVCR2", "PDCD1", "CTLA4", "LAG3", "PDCD1LG2")
geneExp<-read.csv("F:/Ashley_Projects/Metabolism_Immune/01.IdentificationOfKeyMetabolicPathways/10_ssGSEA_immunecheck/data_merged.csv")


mydata<-read.csv("../02_heatmap/SampleID_MetabolicScore.csv")

merged_data<-cbind(mydata[,1:2], geneExp[,5:11])

merged_data$sample_subtypes<-as.factor(merged_data$sample_subtypes)

#############################gene IDO1
q1<-quantile(merged_data$IDO1, 0.01)
q99<-quantile(merged_data$IDO1, 0.99)
for(i in 1:390){
  temp<-merged_data$IDO1[i]
  if(temp<q1){
    merged_data$IDO1[i]<-q1
  }else if(temp>q99){
    merged_data$IDO1[i]<-q99
  }else {
    
  }
}
merged_data[,merged_data$IDO1<q1]$IDO1<-q1
merged_data[,merged_data$IDO1>q99]$IDO1<-q99

#############################gene CD274
q1<-quantile(merged_data$CD274, 0.01)
q99<-quantile(merged_data$CD274, 0.99)
for(i in 1:390){
  temp<-merged_data$CD274[i]
  if(temp<q1){
    merged_data$CD274[i]<-q1
  }else if(temp>q99){
    merged_data$CD274[i]<-q99
  }else {
    
  }
}
merged_data[,merged_data$CD274<q1]$CD274<-q1
merged_data[,merged_data$CD274>q99]$CD274<-q99

#############################gene HAVCR2
q1<-quantile(merged_data$HAVCR2, 0.01)
q99<-quantile(merged_data$HAVCR2, 0.99)
for(i in 1:390){
  temp<-merged_data$HAVCR2[i]
  if(temp<q1){
    merged_data$HAVCR2[i]<-q1
  }else if(temp>q99){
    merged_data$HAVCR2[i]<-q99
  }else {
    
  }
}
merged_data[,merged_data$HAVCR2<q1]$HAVCR2<-q1
merged_data[,merged_data$HAVCR2>q99]$HAVCR2<-q99

#############################gene PDCD1
q1<-quantile(merged_data$PDCD1, 0.01)
q99<-quantile(merged_data$PDCD1, 0.99)
for(i in 1:390){
  temp<-merged_data$PDCD1[i]
  if(temp<q1){
    merged_data$PDCD1[i]<-q1
  }else if(temp>q99){
    merged_data$PDCD1[i]<-q99
  }else {
    
  }
}
merged_data[,merged_data$PDCD1<q1]$PDCD1<-q1
merged_data[,merged_data$PDCD1>q99]$PDCD1<-q99

#############################gene CTLA4
q1<-quantile(merged_data$CTLA4, 0.01)
q99<-quantile(merged_data$CTLA4, 0.99)
for(i in 1:390){
  temp<-merged_data$CTLA4[i]
  if(temp<q1){
    merged_data$CTLA4[i]<-q1
  }else if(temp>q99){
    merged_data$CTLA4[i]<-q99
  }else {
    
  }
}
merged_data[,merged_data$CTLA4<q1]$CTLA4<-q1
merged_data[,merged_data$CTLA4>q99]$CTLA4<-q99

#############################gene LAG3
q1<-quantile(merged_data$LAG3, 0.01)
q99<-quantile(merged_data$LAG3, 0.99)
for(i in 1:390){
  temp<-merged_data$LAG3[i]
  if(temp<q1){
    merged_data$LAG3[i]<-q1
  }else if(temp>q99){
    merged_data$LAG3[i]<-q99
  }else {
    
  }
}
merged_data[,merged_data$LAG3<q1]$LAG3<-q1
merged_data[,merged_data$LAG3>q99]$LAG3<-q99

#############################gene PDCD1LG2
q1<-quantile(merged_data$PDCD1LG2, 0.01)
q99<-quantile(merged_data$PDCD1LG2, 0.99)
for(i in 1:390){
  temp<-merged_data$PDCD1LG2[i]
  if(temp<q1){
    merged_data$PDCD1LG2[i]<-q1
  }else if(temp>q99){
    merged_data$PDCD1LG2[i]<-q99
  }else {
    
  }
}
merged_data[,merged_data$PDCD1LG2<q1]$PDCD1LG2<-q1
merged_data[,merged_data$PDCD1LG2>q99]$PDCD1LG2<-q99


library("reshape2")
library("knitr")
mydata_transf<-melt(merged_data, id.vars = c("X", "sample_subtypes"), variable.names = "ImmuneCheckPoints", 
                    value.name = "GeneExp")
head(mydata_transf)
library("rstatix")
myt_test<-t_test(group_by(mydata_transf,variable), GeneExp~sample_subtypes)
myt_test<-add_significance(myt_test, "p")
myt_test
# A tibble: 7 × 10
# variable  .y.     group1 group2    n1    n2 statistic    df       p p.signif
# <fct>     <chr>   <chr>  <chr>  <int> <int>     <dbl> <dbl>   <dbl> <chr>   
#1 IDO1     GeneExp 1      2        126   264     1.93   170. 0.0553  ns      
#2 CD274    GeneExp 1      2        126   264     0.345  221. 0.73    ns      
#3 HAVCR2   GeneExp 1      2        126   264     3.08   178. 0.00239 **      
#4 PDCD1    GeneExp 1      2        126   264     3.33   129. 0.00113 **      
#5 CTLA4    GeneExp 1      2        126   264     3.22   133. 0.00161 **      
#6 LAG3     GeneExp 1      2        126   264     2.98   145. 0.00336 **      
#7 PDCD1LG2 GeneExp 1      2        126   264    -3.38   263  0.00083 ***


mywilcox_test<-wilcox_test(group_by(mydata_transf, variable), GeneExp~sample_subtypes)
mywilcox_test<-add_significance(mywilcox_test, "p")
mywilcox_test
# A tibble: 7 × 9
#variable  .y.     group1 group2    n1    n2 statistic          p p.signif
#<fct>     <chr>   <chr>  <chr>  <int> <int>     <dbl>      <dbl> <chr>   
#1 IDO1     GeneExp 1      2        126   264    17685  0.0286     *       
#2 CD274    GeneExp 1      2        126   264    16775  0.552      ns      
#3 HAVCR2   GeneExp 1      2        126   264    19470. 0.000624   ***     
#4 PDCD1    GeneExp 1      2        126   264    18876. 0.00000311 ****    
#5 CTLA4    GeneExp 1      2        126   264    18093  0.0000134  ****    
#6 LAG3     GeneExp 1      2        126   264    18278. 0.00647    **      
#7 PDCD1LG2 GeneExp 1      2        126   264    15939  0.0204     *  

head(mydata_transf)
mydata_transf$ExpTran<-log2(mydata_transf$GeneExp+1)

pdf("SevenImmuneCheck_sig.pdf", height = 6, width = 8)

p <- ggboxplot(mydata_transf, x = "variable", y = "GeneExp",fill = "sample_subtypes", 
               color = 'sample_subtypes', palette = c("#00AFBB", "#E7B800"), add = "boxplot", 
               add.params = list(fill = "white")) + 
  labs(x = "Immune check point genes", y = "The gene expression levels", fill = "sample_subtypes") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p + stat_compare_means(aes(group = sample_subtypes), label = "p.signif", paired = FALSE)

dev.off()

pdf("SevenImmuneCheck_sig_log.pdf", height = 6, width = 8)
p <- ggboxplot(mydata_transf, x = "variable", y = "ExpTran",fill = "sample_subtypes", 
               color = 'sample_subtypes', palette = c("#00AFBB", "#E7B800"), add = "boxplot", 
               add.params = list(fill = "white")) + 
  labs(x = "Immune check point genes", y = "The gene expression levels", fill = "sample_subtypes") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p + stat_compare_means(aes(group = sample_subtypes), label = "p.signif", paired = FALSE)
dev.off()



cluster1<-filter(mydata_transf, sample_subtypes == "1")
cluster2<-filter(mydata_transf, sample_subtypes == "2")

cluster1_IDO1<-filter(cluster1, variable == "IDO1")
cluster2_IDO1<-filter(cluster2, variable == "IDO1")

