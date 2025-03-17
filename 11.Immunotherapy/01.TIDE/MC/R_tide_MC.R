setwd("F:/Ashley_Projects/Metabolism_Immune/11.Immunotherapy/01.TIDE/MC/")
rm(list=ls())

LIHC_Expr<-read.csv("F:/Ashley_Projects/Metabolism_Immune/04.MRGsClusters/04_MRG_clusters/geneExp.csv", row.names = 1)
LIHC_Expr[1:5,1:5]


cluster<-read.csv("F:/Ashley_Projects/Metabolism_Immune/03.ClusterTME/01_cluster/Sample_Cluster.csv")
clinical<-read.table("F:/Ashley_Projects/Metabolism_Immune/DataCollection/clinical_rmNA.txt", header = T)
cluster$Id<-substr(cluster$ID, 1, 12)
cluster$MC<-paste0("C", cluster$sample_subtypes)
cluster<-cluster[, c("Id", "MC")]
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
df$id2 <- paste(df$MC, df$Id, sep = '_') # 将风险分组和id串联  
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
library("stringr")

result$Risk <- ifelse(  
  str_sub(result$Patient, 1, 2) == 'C1', 'C1', 'C2'  
)  
result$Risk <- factor(result$Risk, levels = c('C1','C2'))  
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
# A tibble: 6 × 10
#   variable    .y.   group1 group2    n1    n2 statistic    df        p p.signif
#   <fct>       <chr> <chr>  <chr>  <int> <int>     <dbl> <dbl>    <dbl> <chr>   
#  1 CAF         Value C1     C2       136   293     -2.67  278. 7.95e- 3 **      
#  2 Dysfunction Value C1     C2       136   293     -5.10  268. 6.45e- 7 ****    
#  3 Exclusion   Value C1     C2       136   293      7.79  220. 2.58e-13 ****    
#  4 MDSC        Value C1     C2       136   293      9.80  220. 4.73e-19 ****    
#  5 TAM.M2      Value C1     C2       136   293      4.92  264. 1.54e- 6 ****    
#  6 TIDE        Value C1     C2       136   293      8.01  222. 6.40e-14 ****


library("ggplot2")
library("ggpubr")
mydata_transf$Risk<-as.factor(mydata_transf$Risk)

pdf("MetabolismCluster_sig.pdf", height = 6, width = 8)
p <- ggboxplot(mydata_transf, x = "variable", y = "Value",fill = "Risk", 
               color = 'Risk', palette = c("#00AFBB", "#E7B800"), add = "boxplot", 
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


