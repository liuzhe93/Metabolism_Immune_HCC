setwd("F:/Ashley_Projects/Metabolism_Immune/03.ClusterTME/07_boxplot_28immunecells")
rm(list=ls())

mydata<-read.csv("../06_heatmap_28immunecells/Cluster_immunecells.csv", row.names = 1)
mydata_sel<-mydata[c(1,10:37),]
mydata_t<-t(mydata_sel)
mydata_t<-as.data.frame(mydata_t)
mydata_t$sample_subtypes<-as.factor(mydata_t$sample_subtypes)
mydata_t$X<-row.names(mydata_t)

library("reshape2")
library("knitr")
mydata_transf<-melt(mydata_t, id.vars = c("X", "sample_subtypes"), variable.names = "ImmuneCells", 
                    value.name = "MetabolicScore")
head(mydata_transf)
library("rstatix")
myt_test<-t_test(group_by(mydata_transf,variable), MetabolicScore~sample_subtypes)
myt_test<-add_significance(myt_test, "p")
myt_test
# A tibble: 28 × 10
#   variable                    .y.            group1 group2    n1    n2 statistic    df        p p.signif
#   <fct>                       <chr>          <chr>  <chr>  <int> <int>     <dbl> <dbl>    <dbl> <chr>   
#1  Activated B cell            MetabolicScore 1      2        120   235      1.04  161. 3   e- 1 ns      
#2  Activated CD4 T cell        MetabolicScore 1      2        120   235      9.02  246. 5.58e-17 ****    
#3  Activated CD8 T cell        MetabolicScore 1      2        120   235      2.39  176. 1.77e- 2 *       
#4  Central memory CD4 T cell   MetabolicScore 1      2        120   235      6.50  200. 6.21e-10 ****    
#5  Central memory CD8 T cell   MetabolicScore 1      2        120   235     -1.64  187. 1.03e- 1 ns      
#6  Effector memeory CD4 T cell MetabolicScore 1      2        120   235      5.59  191. 7.90e- 8 ****    
#7  Effector memeory CD8 T cell MetabolicScore 1      2        120   235     -3.44  186. 7.21e- 4 ***     
#8  Gamma delta T cell          MetabolicScore 1      2        120   235      1.25  194. 2.15e- 1 ns      
#9  Immature  B cell            MetabolicScore 1      2        120   235     -1.92  205. 5.62e- 2 ns      
#10 Memory B cell               MetabolicScore 1      2        120   235     -2.25  189. 2.53e- 2 *   
# ℹ 18 more rows
# ℹ Use `print(n = ...)` to see more rows
mywilcox_test<-wilcox_test(group_by(mydata_transf, variable), MetabolicScore~sample_subtypes)
mywilcox_test<-add_significance(mywilcox_test, "p")
mywilcox_test
# A tibble: 28 × 9
#variable                    .y.            group1 group2    n1    n2 statistic        p p.signif
#<fct>                       <chr>          <chr>  <chr>  <int> <int>     <dbl>    <dbl> <chr>   
#1  Activated B cell            MetabolicScore 1      2        120   235     13702 6.64e- 1 ns      
#2  Activated CD4 T cell        MetabolicScore 1      2        120   235     21308 3.27e-15 ****    
#3  Activated CD8 T cell        MetabolicScore 1      2        120   235     15284 1.96e- 1 ns      
#4  Central memory CD4 T cell   MetabolicScore 1      2        120   235     19575 2.16e- 9 ****    
#5  Central memory CD8 T cell   MetabolicScore 1      2        120   235     12943 2.06e- 1 ns      
#6  Effector memeory CD4 T cell MetabolicScore 1      2        120   235     19338 1.03e- 8 ****    
#7  Effector memeory CD8 T cell MetabolicScore 1      2        120   235     10895 4.59e- 4 ***     
#8  Gamma delta T cell          MetabolicScore 1      2        120   235     15415 1.51e- 1 ns      
#9  Immature  B cell            MetabolicScore 1      2        120   235     12527 8.56e- 2 ns      
#10 Memory B cell               MetabolicScore 1      2        120   235     11695 8.57e- 3 **      
# ℹ 18 more rows
# ℹ Use `print(n = ...)` to see more rows

library("ggplot2")
library("ggpubr")
pdf("Immunecells_sig.pdf", height = 6, width = 8)
p <- ggboxplot(mydata_transf, x = "variable", y = "MetabolicScore",fill = "sample_subtypes", 
               color = 'sample_subtypes', palette = c("#00AFBB", "#E7B800"), add = "boxplot", 
               add.params = list(fill = "white")) + 
  labs(x = "28 immune cell infiltration", y = "Relative enrichment score", fill = "sample_subtypes") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p + stat_compare_means(aes(group = sample_subtypes), label = "p.signif", paired = FALSE)
dev.off()



