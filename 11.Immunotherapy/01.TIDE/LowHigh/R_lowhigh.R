setwd("F:/Ashley_Projects/Metabolism_Immune/11.Immunotherapy/01.TIDE/LowHigh")
rm(list=ls())


rt<-read.table("F:/Ashley_Projects/Metabolism_Immune/06.MR_relatedRiskModel/01_RiskModel/risk.txt")
rt=rt[,c(1,2,10)]
head(rt)
library("survminer")
res.cut<-surv_cutpoint(rt,
                       time = "futime",
                       event = "fustat",
                       variables = "riskScore")
summary(res.cut)
#          cutpoint statistic
#riskScore 1.387188  7.648931
res.cat<-surv_categorize(res.cut)
head(res.cat)
table(res.cat$riskScore)
#high  low 
#  73  282

identical(row.names(rt), row.names(res.cat))
#[1] TRUE
sample_cluster<-cbind(rownames(res.cat), rt$riskScore, res.cat$riskScore)
colnames(sample_cluster)<-c("SampleName", "riskScore", "riskClass")
write.csv(sample_cluster, "Sample_Cluster.csv", row.names = F, quote = F)


rm(list=ls())
LIHC_Expr<-read.csv("F:/Ashley_Projects/Metabolism_Immune/04.MRGsClusters/04_MRG_clusters/geneExp.csv", row.names = 1)
LIHC_Expr[1:5,1:5]


cluster<-read.csv("Sample_Cluster.csv")
clinical<-read.table("F:/Ashley_Projects/Metabolism_Immune/DataCollection/clinical_rmNA.txt", header = T)
cluster$Id<-substr(cluster$SampleName, 1, 12)
cluster$RL<-ifelse(cluster$riskClass=="low", "Low", "High")
cluster<-cluster[, c("Id", "RL")]
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
df$id2 <- paste(df$RL, df$Id, sep = '_') # 将风险分组和id串联  
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
  temp<-substr(result[i,1], 1, 3)
  result$Risk[i]<-temp
}
result$Risk<-ifelse(result$Risk=="Low", "Low", "High")
result$Risk <- factor(result$Risk, levels = c('Low','High'))  
head(result)
library("reshape2")
library("knitr")
mydata<-subset(result, select = c("Risk", "Responder"))

RL_SP<-table(mydata$Risk, mydata$Responder)
RL_NO<-addmargins(RL_SP)
RL_NO_matrix<-as.matrix(RL_NO)[1:2,1:2]
RL_NO_df<-as.data.frame.array(RL_NO_matrix)
RL_NO_df
RL_NO_df$RL<-row.names(RL_NO_df)
library(tidyverse)
RL_count <- RL_NO_df %>% pivot_longer(cols=c(False:True),
                                      names_to = 'Responder',
                                      values_to = 'freq')
RL_count

library(plyr)
RL_count_transfer = ddply(RL_count,'RL',transform,percent_con=freq/sum(freq)*100)
RL_count_transfer

pdf("RL_percentage.pdf")
p3 <- ggplot(RL_count_transfer,aes(x=RL,y=percent_con,fill=Responder))+
  geom_bar(stat = 'identity',width = 0.5,colour='black')
p4 <- p3+labs(x='RL',y='Percentage')+
  theme(axis.title = element_text(size=12),
        axis.text = element_text(size=11))+
  scale_y_continuous(breaks=seq(0,100,25),
                     labels=c('0','25%','50%','75%','100%'))

p4
dev.off()






