setwd("F:/Ashley_Projects/Metabolism_Immune/11.Immunotherapy/01.TIDE/GC")
rm(list=ls())

LIHC_Expr<-read.csv("F:/Ashley_Projects/Metabolism_Immune/04.MRGsClusters/04_MRG_clusters/geneExp.csv", row.names = 1)
LIHC_Expr[1:5,1:5]


cluster<-read.csv("F:/Ashley_Projects/Metabolism_Immune/04.MRGsClusters/04_MRG_clusters/Sample_Cluster.csv")
clinical<-read.table("F:/Ashley_Projects/Metabolism_Immune/DataCollection/clinical_rmNA.txt", header = T)
cluster$Id<-substr(cluster$sample_name, 1, 12)
cluster$GC<-paste0("C", cluster$cluster)
cluster<-cluster[, c("Id", "GC")]
dim(cluster)
dim(clinical)
clin_info<-merge(clinical, cluster, by = "Id")
dim(clin_info)
#[1] 355  10

dt <- as.data.frame(t(LIHC_Expr)) 
dt$Id <- substr(row.names(dt), 1, 12)
dt$Id <- gsub("\\.", "-", dt$Id)

merged_data<-merge(clin_info, dt, by = "Id")
df<-merged_data
df$id2 <- paste(df$GC.x, df$Id, sep = '_') # 将风险分组和id串联  
for(i in 1:nrow(df)){
  df$id2[i]<-paste(df$id2[i], i, sep = "_")
}
rownames(df) <- df$id2
df <- df[, c(11:12935)]
df2 <- t(df)
df2[1:4,1:4]
Expr <- t(apply(df2, 1, function(x){x-(mean(x))})) # 均值标准化  
Expr[1:6,1:6]
write.table(Expr, file = 'TIDE.txt', sep = "\t", quote = F, row.names = T) # 矩阵保存到本地
#http://tide.dfci.harvard.edu/login/

result <- read.csv('TIDE_res.csv')
colnames(result)
result$Risk<-result$Patient
for(i in 1:nrow(result)){
  temp<-substr(result[i,1], 1, 2)
  result$Risk[i]<-temp
}
result$Risk <- factor(result$Risk, levels = c('C1','C2','C3'))  
head(result)
library("reshape2")
library("knitr")
mydata<-subset(result, select = c("Patient", "CAF", "Dysfunction", "Exclusion", "MDSC", "TAM.M2", "TIDE", "Risk"))
head(mydata)
mydata_transf<-melt(mydata, id.vars = c("Patient", "Risk"), variable.names = "TIDE_related_scores", 
                    value.name = "Value")
head(mydata_transf)
library("rstatix")
myt_test<-t_test(group_by(mydata_transf,variable), Value~Risk)
myt_test<-add_significance(myt_test, "p")
myt_test
# A tibble: 18 × 12
# variable       .y.   group1 group2    n1    n2 statistic    df        p    p.adj p.adj.signif p.signif
# <fct>          <chr> <chr>  <chr>  <int> <int>     <dbl> <dbl>    <dbl>    <dbl> <chr>        <chr>   
#  1 CAF          Value C1     C2        87   222   -2.17    158. 3.2 e- 2 9.5 e- 2 ns           *       
#  2 CAF          Value C1     C3        87   120   -1.78    204. 7.6 e- 2 1.52e- 1 ns           ns      
#  3 CAF          Value C2     C3       222   120   -0.0802  197. 9.36e- 1 9.36e- 1 ns           ns      
#  4 Dysfunction  Value C1     C2        87   222   -6.44    177. 1.10e- 9 3.30e- 9 ****         ****    
#  5 Dysfunction  Value C1     C3        87   120   -3.09    197. 2   e- 3 5   e- 3 **           **      
#  6 Dysfunction  Value C2     C3       222   120    3.03    245. 3   e- 3 5   e- 3 **           **      
#  7 Exclusion    Value C1     C2        87   222    9.89    128. 1.82e-17 5.46e-17 ****         ****    
#  8 Exclusion    Value C1     C3        87   120    5.51    161. 1.36e- 7 2.72e- 7 ****         ****    
#  9 Exclusion    Value C2     C3       222   120   -5.07    232. 8.27e- 7 8.27e- 7 ****         ****    
#  10 MDSC        Value C1     C2        87   222   12.6     144. 3.92e-25 1.18e-24 ****         ****    
#  11 MDSC        Value C1     C3        87   120    7.70    177. 9.21e-13 1.84e-12 ****         ****    
#  12 MDSC        Value C2     C3       222   120   -4.83    238. 2.46e- 6 2.46e- 6 ****         ****    
#  13 TAM.M2      Value C1     C2        87   222    5.46    166. 1.69e- 7 5.07e- 7 ****         ****    
#  14 TAM.M2      Value C1     C3        87   120    3.47    180. 6.54e- 4 1   e- 3 **           ***     
#  15 TAM.M2      Value C2     C3       222   120   -1.99    266. 4.8 e- 2 4.8 e- 2 *            *       
#  16 TIDE        Value C1     C2        87   222   10.3     130. 1.65e-18 4.95e-18 ****         ****    
#  17 TIDE        Value C1     C3        87   120    5.79    164. 3.5 e- 8 7   e- 8 ****         ****    
#  18 TIDE        Value C2     C3       222   120   -5.07    232. 8.27e- 7 8.27e- 7 ****         ****


library("ggplot2")
library("ggpubr")
mydata_transf$Risk<-as.factor(mydata_transf$Risk)

pdf("GeneCluster_sig.pdf", height = 6, width = 8)
p <- ggboxplot(mydata_transf, x = "variable", y = "Value",fill = "Risk", 
               color = 'Risk', palette = c("blue", "red", "green"), add = "boxplot", 
               add.params = list(fill = "white")) + 
  labs(x = "TIDE-related scores", y = "Value", fill = "Risk") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p + stat_compare_means(aes(group = Risk), label = "p.signif", paired = FALSE)
dev.off()





library("tidyr")
library("dplyr")
mydata<-subset(result, select = c("Patient", "CAF", "Dysfunction", "Exclusion", "MDSC", "TAM.M2", "TIDE", "Risk"))
head(mydata)

tmp_mydata <- mydata %>% mutate(Risk = factor(Risk))
library("ggplot2")
library("ggunchained")
library("reshape2")
pdf("TIDE_vilplot.pdf")
tmp_mydata %>% 
  mutate(Risk = factor(Risk)) %>% 
  pivot_longer(-c(Patient,Risk), names_to = "TIDE_related_scores",values_to = "value") %>% 
  ggplot(aes(TIDE_related_scores,value,fill=Risk))+
  geom_split_violin(draw_quantiles = c(0.25, 0.5, 0.75))+
  theme_bw()+
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45,hjust = 1)
  )
dev.off()


