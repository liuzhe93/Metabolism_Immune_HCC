setwd("F:/Ashley_Projects/Metabolism_Immune/07.MR_ImmInfil")
rm(list=ls())

dataExp<-read.table("F:/Ashley_Projects/Metabolism_Immune/06.MR_relatedRiskModel/05_EstimateScores/dataExp.txt", 
                    header = T, row.names = 1)
genelist<-read.csv("F:/Ashley_Projects/Metabolism_Immune/01.IdentificationOfKeyMetabolicPathways/09_ssGSEA_immunecells/ImmunecellMetagenes.csv")
#28 immune cell types and metagenes are collected from https://doi.org/10.1016/j.celrep.2016.12.019
gene_set<-genelist[, 1:2]
head(gene_set)
#  Metagene        Cell.type
#1   ADAM28 Activated B cell
#2    CD180 Activated B cell
#3    CD79B Activated B cell
#4      BLK Activated B cell
#5     CD19 Activated B cell
#6    MS4A1 Activated B cell
list<-split(as.matrix(gene_set)[,1], gene_set[,2])
list

library("GSVA")
ssgsea_par <- ssgseaParam(as.matrix(dataExp), list)
gsva_matrix <- gsva(ssgsea_par)
dim(gsva_matrix)
#[1]  28 70
gsva_matrix<-as.data.frame(gsva_matrix)
gsva_matrix_t<-t(gsva_matrix)
dim(gsva_matrix_t)
#[1] 70  28
gsva_matrix_t[1:4,1:4]
#      Activated B cell Activated CD4 T cell Activated CD8 T cell Activated dendritic cell
#L001P        0.7777076            0.7827752             1.087999                0.8707364
#L003P        0.7462330            0.8594134             1.093309                0.9230732
#L004P        0.7665463            0.8788208             1.093421                0.8802251
#L005P        0.6899906            0.8524650             1.038292                0.8967997
write.csv(gsva_matrix_t, "data_merged.csv", quote = F)
gsva_matrix_t<-as.data.frame(gsva_matrix_t)
gsva_matrix_t$X<-row.names(gsva_matrix_t)

mydata<-read.csv("F:/Ashley_Projects/Metabolism_Immune/06.MR_relatedRiskModel/05_EstimateScores/Merged_data.csv")
dim(mydata)
#[1] 70 20
mydata[1:4,15:20]
#  SEPTIN11      ZYX   ZZEF1     ZZZ3 riskScore riskCLass
#1  7.49662 19.31170 4.04251  5.01862 0.6973069       low
#2 24.30710 48.29770 4.36698  9.99004 0.4608506       low
#3 21.12940  8.37228 2.25808 11.12160 6.8012530      high
#4 36.29560 29.28320 2.27734 17.21810 6.6360938       low
merged_data<-merge(mydata, gsva_matrix_t, by = "X")
row.names(merged_data)<-merged_data$X
merged_data$sample_name<-NULL
merged_data$riskCLass<-as.factor(merged_data$riskCLass)
merged_data_sel<-merged_data[, c(1, 19, 20:47)]


library("reshape2")
library("knitr")
mydata_transf<-melt(merged_data_sel, id.vars = c("X", "riskCLass"), variable.names = "ImmuneCells", 
                    value.name = "EnrichmentScore")
head(mydata_transf)
library("rstatix")
myt_test<-t_test(group_by(mydata_transf,variable), EnrichmentScore~riskCLass)
myt_test<-add_significance(myt_test, "p")
myt_test
## A tibble: 28 × 10
#variable                       .y.             group1 group2    n1    n2 statistic    df      p p.signif
#<fct>                          <chr>           <chr>  <chr>  <int> <int>     <dbl> <dbl>  <dbl> <chr>   
#1 Activated B cell               EnrichmentScore high   low        7    63    -1.22   7.39 0.259  ns      
#2 Activated CD4 T cell           EnrichmentScore high   low        7    63     2.97  10.7  0.013  *       
#3 Activated CD8 T cell           EnrichmentScore high   low        7    63    -0.407  7.00 0.696  ns      
#4 Activated dendritic cell       EnrichmentScore high   low        7    63    -0.102 12.8  0.92   ns      
#5 CD56bright natural killer cell EnrichmentScore high   low        7    63    -0.196  6.72 0.85   ns      
#6 CD56dim natural killer cell    EnrichmentScore high   low        7    63    -2.64   6.93 0.0337 *       
#7 Central memory CD4 T cell      EnrichmentScore high   low        7    63    -2.25  20.0  0.0357 *       
#8 Central memory CD8 T cell      EnrichmentScore high   low        7    63     0.797  7.71 0.449  ns      
#9 Effector memeory CD4 T cell    EnrichmentScore high   low        7    63     0.581  7.72 0.578  ns      
#10 Effector memeory CD8 T cell    EnrichmentScore high   low        7    63    -1.68   7.08 0.136  ns      
## ℹ 18 more rows
## ℹ Use `print(n = ...)` to see more rows

mywilcox_test<-wilcox_test(group_by(mydata_transf, variable), EnrichmentScore~riskCLass)
mywilcox_test<-add_significance(mywilcox_test, "p")
mywilcox_test
## A tibble: 28 × 9
#variable                       .y.             group1 group2    n1    n2 statistic      p p.signif
#<fct>                          <chr>           <chr>  <chr>  <int> <int>     <dbl>  <dbl> <chr>   
#1 Activated B cell               EnrichmentScore high   low        7    63       148 0.159  ns      
#2 Activated CD4 T cell           EnrichmentScore high   low        7    63       327 0.038  *       
#3 Activated CD8 T cell           EnrichmentScore high   low        7    63       203 0.739  ns      
#4 Activated dendritic cell       EnrichmentScore high   low        7    63       206 0.784  ns      
#5 CD56bright natural killer cell EnrichmentScore high   low        7    63       225 0.938  ns      
#6 CD56dim natural killer cell    EnrichmentScore high   low        7    63        97 0.016  *       
#7 Central memory CD4 T cell      EnrichmentScore high   low        7    63       154 0.196  ns      
#8 Central memory CD8 T cell      EnrichmentScore high   low        7    63       258 0.469  ns      
#9 Effector memeory CD4 T cell    EnrichmentScore high   low        7    63       241 0.695  ns      
#10 Effector memeory CD8 T cell    EnrichmentScore high   low        7    63       131 0.0814 ns      
## ℹ 18 more rows
## ℹ Use `print(n = ...)` to see more rows

library("ggplot2")
library("ggpubr")
pdf("Immunecells_sig.pdf", height = 6, width = 8)
p <- ggboxplot(mydata_transf, x = "variable", y = "EnrichmentScore",fill = "riskCLass", 
               color = 'riskCLass', palette = c("#00AFBB", "#E7B800"), add = "boxplot", 
               add.params = list(fill = "white")) + 
  labs(x = "28 immune cell infiltration", y = "Relative enrichment score", fill = "riskCLass") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p + stat_compare_means(aes(group = riskCLass), label = "p.signif", paired = FALSE)
dev.off()






