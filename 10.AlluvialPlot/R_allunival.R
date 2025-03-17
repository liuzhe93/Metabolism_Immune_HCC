setwd("F:/Ashley_Projects/Metabolism_Immune/10.AlluvialPlot")
rm(list=ls())

library("ggalluvial")

metabolism_cluster<-read.csv("F:/Ashley_Projects/Metabolism_Immune/03.ClusterTME/01_cluster/Sample_Cluster.csv")
head(metabolism_cluster)
#                            ID sample_subtypes
#1 TCGA-FV-A3I0-01A-11R-A22L-07               1
#2 TCGA-BD-A3ER-01A-11R-A213-07               2
#3 TCGA-CC-5261-01A-01R-A131-07               1
#4 TCGA-DD-AAVZ-01A-11R-A41C-07               2
#5 TCGA-DD-AADN-01A-11R-A41C-07               1
#6 TCGA-DD-A1EB-01A-11R-A131-07               2
gene_cluster<-read.csv("F:/Ashley_Projects/Metabolism_Immune/04.MRGsClusters/04_MRG_clusters/Sample_Cluster.csv")
head(gene_cluster)
#  cluster                  sample_name
#1       1 TCGA-FV-A3I0-01A-11R-A22L-07
#2       2 TCGA-BD-A3ER-01A-11R-A213-07
#3       1 TCGA-CC-5261-01A-01R-A131-07
#4       3 TCGA-DD-AAVZ-01A-11R-A41C-07
#5       1 TCGA-DD-AADN-01A-11R-A41C-07
#6       2 TCGA-DD-A1EB-01A-11R-A131-07

risk_level<-read.csv("F:/Ashley_Projects/Metabolism_Immune/06.MR_relatedRiskModel/01_RiskModel/Sample_Cluster.csv")
head(risk_level)
#      SampleName fustat riskClass
#1 TCGA-2Y-A9GT_1      1       low
#2 TCGA-2Y-A9GU_2      0       low
#3 TCGA-2Y-A9GV_3      1       low
#4 TCGA-2Y-A9GW_4      1       low
#5 TCGA-2Y-A9GY_5      1       low
#6 TCGA-2Y-A9GZ_6      1       low

MC_data<-metabolism_cluster
MC_data$ID<-substr(MC_data$ID, 1, 12)
colnames(MC_data)<-c("ID", "Cluster")
for(i in 1:nrow(MC_data)){
  if(MC_data$Cluster[i] == 1){
    MC_data$Cluster[i]="C1"
  }else{
    MC_data$Cluster[i]="C2"
  }
}
head(MC_data)


GC_data<-gene_cluster
GC_data$sample_name<-substr(GC_data$sample_name, 1, 12)
colnames(GC_data)<-c("Cluster", "ID")
for(i in 1:nrow(GC_data)){
  if(GC_data$Cluster[i] == 1){
    GC_data$Cluster[i]="G1"
  }else if(GC_data$Cluster[i] == 2){
    GC_data$Cluster[i]="G2"
  }else{
    GC_data$Cluster[i]="G3"
  }
}
head(GC_data)

RL_data<-risk_level
RL_data$SampleName<-substr(RL_data$SampleName, 1, 12)
colnames(RL_data)<-c("ID", "Survived", "Cluster")
for(i in 1:nrow(RL_data)){
  if(RL_data$Cluster[i] == "low"){
    RL_data$Cluster[i]="Low"
  }else{
    RL_data$Cluster[i]="High"
  }
}
for(i in 1:nrow(RL_data)){
  if(RL_data$Survived[i] == "1"){
    RL_data$Survived[i]="No"
  }else{
    RL_data$Survived[i]="Yes"
  }
}
head(RL_data)

merged_data1<-merge(MC_data, GC_data, by = "ID")
merged_data2<-merge(merged_data1, RL_data, by = "ID")

merged_data2$Freq<-1
colnames(merged_data2)[c(2, 3, 5)]<-c("MC", "GC", "RL")
# 直接使用原始数据（宽数据）绘图；
# 绘制alluvium，width调整结点处条带水平宽度；knot.pos调整曲率；reverse调整着色顺序；
df<-merged_data2
write.csv(df, "AP_data.csv", quote = F, row.names = F)

pdf("AlluvivalDiagram.pdf")
mydata<-read.csv("AP_data_processed.csv", header = T)
ggplot(mydata, aes(y = Freq, axis1 = MC, axis2 = GC, axis3 = RL)) +
  geom_alluvium(aes(fill = Survived), width = 0, knot.pos = 0, reverse = FALSE) + 
  guides(fill = FALSE) +geom_stratum(width= 1/8, reverse = FALSE) +geom_text(stat = "stratum", infer.label =TRUE, reverse = FALSE) +
  scale_x_continuous(breaks = 1:3, labels =c("MC", "GC", "RL")) +
  ggtitle("AlluvialDiagram")
dev.off()
#+coord_flip() 

#data(vaccinations)
#levels(vaccinations$response) <- rev(levels(vaccinations$response))
#ggplot(vaccinations,
#       aes(x = survey, stratum = response, alluvium = subject,
#           y = freq,
#           fill = response, label = response)) +
#  scale_x_discrete(expand = c(.1, .1)) +
#  geom_flow() +
#  geom_stratum(alpha = .5) +
#  geom_text(stat = "stratum", size = 3) +
#  theme(legend.position = "none") +
#  ggtitle("vaccination survey responses at three points in time")



#ggplot(data = df, aes(axis1 = MC, axis2 = GC, axis3 = RL, weight = Freq)) +
#  scale_x_discrete(limits = c("MC", "GC", "RL"), expand = c(.1, .05)) +
#  geom_alluvium(aes(fill = Survived)) +
#  geom_stratum() + geom_text(stat = "stratum", label.strata = TRUE) +
#  theme_minimal() +
#  ggtitle("Plot",
#          "stratified by different clusters and survival")






