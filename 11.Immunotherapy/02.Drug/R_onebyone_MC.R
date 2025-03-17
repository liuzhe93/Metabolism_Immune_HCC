setwd("F:/Ashley_Projects/Metabolism_Immune/11.Immunotherapy/02.Drug")
rm(list=ls())

metabolism_cluster<-read.csv("F:/Ashley_Projects/Metabolism_Immune/03.ClusterTME/01_cluster/Sample_Cluster.csv")
head(metabolism_cluster)

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

ic50<-read.csv("drugSensitivity.csv")
dim(ic50)
#[1] 390  61

for(i in 1:nrow(ic50)){
  ic50$sample_name[i]<-substring(ic50$X[i], 1, 12)
  ic50$sample_name[i]<-gsub("\\.", "-", ic50$sample_name[i])
}

merged<-merge(MC_data, ic50, by.x = "ID", by.y = "sample_name")
dim(merged)
#[1] 480  63

for(i in 1:nrow(merged)){
  row.names(merged)[i]<-paste(merged$ID[i], i, sep = "_")
}


merged_format<-cbind(merged[, 2], merged[, 4:ncol(merged)])
colnames(merged_format)[1]<-"MC"


####################################drug1： Sorafenib################################################
drug_sen_Sorafenib<-cbind(merged_format$MC,as.numeric(merged_format$Sorafenib))
rownames(drug_sen_Sorafenib)<-row.names(merged_format)
colnames(drug_sen_Sorafenib)<-c("Risk","IC50")
head(drug_sen_Sorafenib)
drug_sen_Sorafenib<-as.data.frame(drug_sen_Sorafenib)
drug_sen_Sorafenib$IC50<-as.numeric(drug_sen_Sorafenib$IC50)
library(ggpubr)
pdf("Sorafenib.pdf")
p <- ggboxplot(drug_sen_Sorafenib, x="Risk",
               y = "IC50", color = "Risk",
               palette = "jco", add = "jitter")
# 添加p值
p + stat_compare_means()
dev.off()

####################################drug2： Erlotinib################################################
drug_sen_Erlotinib<-cbind(merged_format$MC,as.numeric(merged_format$Erlotinib))
rownames(drug_sen_Erlotinib)<-row.names(merged_format)
colnames(drug_sen_Erlotinib)<-c("Risk","IC50")
head(drug_sen_Erlotinib)
drug_sen_Erlotinib<-as.data.frame(drug_sen_Erlotinib)
drug_sen_Erlotinib$IC50<-as.numeric(drug_sen_Erlotinib$IC50)
library(ggpubr)
pdf("Erlotinib.pdf")
p <- ggboxplot(drug_sen_Erlotinib, x="Risk",
               y = "IC50", color = "Risk",
               palette = "jco", add = "jitter")
# 添加p值
p + stat_compare_means()
dev.off()

####################################drug3： Rapamycin################################################
drug_sen_Rapamycin<-cbind(merged_format$MC,as.numeric(merged_format$Rapamycin))
rownames(drug_sen_Rapamycin)<-row.names(merged_format)
colnames(drug_sen_Rapamycin)<-c("Risk","IC50")
head(drug_sen_Rapamycin)
drug_sen_Rapamycin<-as.data.frame(drug_sen_Rapamycin)
drug_sen_Rapamycin$IC50<-as.numeric(drug_sen_Rapamycin$IC50)
library(ggpubr)
pdf("Rapamycin.pdf")
p <- ggboxplot(drug_sen_Rapamycin, x="Risk",
               y = "IC50", color = "Risk",
               palette = "jco", add = "jitter")
# 添加p值
p + stat_compare_means()
dev.off()

####################################drug4： Paclitaxel################################################
drug_sen_Paclitaxel<-cbind(merged_format$MC,as.numeric(merged_format$Paclitaxel))
rownames(drug_sen_Paclitaxel)<-row.names(merged_format)
colnames(drug_sen_Paclitaxel)<-c("Risk","IC50")
head(drug_sen_Paclitaxel)
drug_sen_Paclitaxel<-as.data.frame(drug_sen_Paclitaxel)
drug_sen_Paclitaxel$IC50<-as.numeric(drug_sen_Paclitaxel$IC50)
library(ggpubr)
pdf("Paclitaxel.pdf")
p <- ggboxplot(drug_sen_Paclitaxel, x="Risk",
               y = "IC50", color = "Risk",
               palette = "jco", add = "jitter")
# 添加p值
p + stat_compare_means()
dev.off()

####################################drug5： Cyclopamine################################################
drug_sen_Cyclopamine<-cbind(merged_format$MC,as.numeric(merged_format$Cyclopamine))
rownames(drug_sen_Cyclopamine)<-row.names(merged_format)
colnames(drug_sen_Cyclopamine)<-c("Risk","IC50")
head(drug_sen_Cyclopamine)
drug_sen_Cyclopamine<-as.data.frame(drug_sen_Cyclopamine)
drug_sen_Cyclopamine$IC50<-as.numeric(drug_sen_Cyclopamine$IC50)
library(ggpubr)
pdf("Cyclopamine.pdf")
p <- ggboxplot(drug_sen_Cyclopamine, x="Risk",
               y = "IC50", color = "Risk",
               palette = "jco", add = "jitter")
# 添加p值
p + stat_compare_means()
dev.off()

####################################drug6： AZ628################################################
drug_sen_AZ628<-cbind(merged_format$MC,as.numeric(merged_format$AZ628))
rownames(drug_sen_AZ628)<-row.names(merged_format)
colnames(drug_sen_AZ628)<-c("Risk","IC50")
head(drug_sen_AZ628)
drug_sen_AZ628<-as.data.frame(drug_sen_AZ628)
drug_sen_AZ628$IC50<-as.numeric(drug_sen_AZ628$IC50)
library(ggpubr)
pdf("AZ628.pdf")
p <- ggboxplot(drug_sen_AZ628, x="Risk",
               y = "IC50", color = "Risk",
               palette = "jco", add = "jitter")
# 添加p值
p + stat_compare_means()
dev.off()

####################################drug7： Imatinib################################################
drug_sen_Imatinib<-cbind(merged_format$MC,as.numeric(merged_format$Imatinib))
rownames(drug_sen_Imatinib)<-row.names(merged_format)
colnames(drug_sen_Imatinib)<-c("Risk","IC50")
head(drug_sen_Imatinib)
drug_sen_Imatinib<-as.data.frame(drug_sen_Imatinib)
drug_sen_Imatinib$IC50<-as.numeric(drug_sen_Imatinib$IC50)
library(ggpubr)
pdf("Imatinib.pdf")
p <- ggboxplot(drug_sen_Imatinib, x="Risk",
               y = "IC50", color = "Risk",
               palette = "jco", add = "jitter")
# 添加p值
p + stat_compare_means()
dev.off()

####################################drug8： CMK################################################
drug_sen_CMK<-cbind(merged_format$MC,as.numeric(merged_format$CMK))
rownames(drug_sen_CMK)<-row.names(merged_format)
colnames(drug_sen_CMK)<-c("Risk","IC50")
head(drug_sen_CMK)
drug_sen_CMK<-as.data.frame(drug_sen_CMK)
drug_sen_CMK$IC50<-as.numeric(drug_sen_CMK$IC50)
library(ggpubr)
pdf("CMK.pdf")
p <- ggboxplot(drug_sen_CMK, x="Risk",
               y = "IC50", color = "Risk",
               palette = "jco", add = "jitter")
# 添加p值
p + stat_compare_means()
dev.off()

####################################drug9： Pyrimethamine################################################
drug_sen_Pyrimethamine<-cbind(merged_format$MC,as.numeric(merged_format$Pyrimethamine))
rownames(drug_sen_Pyrimethamine)<-row.names(merged_format)
colnames(drug_sen_Pyrimethamine)<-c("Risk","IC50")
head(drug_sen_Pyrimethamine)
drug_sen_Pyrimethamine<-as.data.frame(drug_sen_Pyrimethamine)
drug_sen_Pyrimethamine$IC50<-as.numeric(drug_sen_Pyrimethamine$IC50)
library(ggpubr)
pdf("Pyrimethamine.pdf")
p <- ggboxplot(drug_sen_Pyrimethamine, x="Risk",
               y = "IC50", color = "Risk",
               palette = "jco", add = "jitter")
# 添加p值
p + stat_compare_means()
dev.off()

####################################drug10： GW843682X################################################
drug_sen_GW843682X<-cbind(merged_format$MC,as.numeric(merged_format$GW843682X))
rownames(drug_sen_GW843682X)<-row.names(merged_format)
colnames(drug_sen_GW843682X)<-c("Risk","IC50")
head(drug_sen_GW843682X)
drug_sen_GW843682X<-as.data.frame(drug_sen_GW843682X)
drug_sen_GW843682X$IC50<-as.numeric(drug_sen_GW843682X$IC50)
library(ggpubr)
pdf("GW843682X.pdf")
p <- ggboxplot(drug_sen_GW843682X, x="Risk",
               y = "IC50", color = "Risk",
               palette = "jco", add = "jitter")
# 添加p值
p + stat_compare_means()
dev.off()

####################################drug11： Parthenolide################################################
drug_sen_Parthenolide<-cbind(merged_format$MC,as.numeric(merged_format$Parthenolide))
rownames(drug_sen_Parthenolide)<-row.names(merged_format)
colnames(drug_sen_Parthenolide)<-c("Risk","IC50")
head(drug_sen_Parthenolide)
drug_sen_Parthenolide<-as.data.frame(drug_sen_Parthenolide)
drug_sen_Parthenolide$IC50<-as.numeric(drug_sen_Parthenolide$IC50)
library(ggpubr)
pdf("Parthenolide.pdf")
p <- ggboxplot(drug_sen_Parthenolide, x="Risk",
               y = "IC50", color = "Risk",
               palette = "jco", add = "jitter")
# 添加p值
p + stat_compare_means()
dev.off()

####################################drug12： Roscovitine################################################
drug_sen_Roscovitine<-cbind(merged_format$MC,as.numeric(merged_format$Roscovitine))
rownames(drug_sen_Roscovitine)<-row.names(merged_format)
colnames(drug_sen_Roscovitine)<-c("Risk","IC50")
head(drug_sen_Roscovitine)
drug_sen_Roscovitine<-as.data.frame(drug_sen_Roscovitine)
drug_sen_Roscovitine$IC50<-as.numeric(drug_sen_Roscovitine$IC50)
library(ggpubr)
pdf("Roscovitine.pdf")
p <- ggboxplot(drug_sen_Roscovitine, x="Risk",
               y = "IC50", color = "Risk",
               palette = "jco", add = "jitter")
# 添加p值
p + stat_compare_means()
dev.off()

####################################drug13： Salubrinal################################################
drug_sen_Salubrinal<-cbind(merged_format$MC,as.numeric(merged_format$Salubrinal))
rownames(drug_sen_Salubrinal)<-row.names(merged_format)
colnames(drug_sen_Salubrinal)<-c("Risk","IC50")
head(drug_sen_Salubrinal)
drug_sen_Salubrinal<-as.data.frame(drug_sen_Salubrinal)
drug_sen_Salubrinal$IC50<-as.numeric(drug_sen_Salubrinal$IC50)
library(ggpubr)
pdf("Salubrinal.pdf")
p <- ggboxplot(drug_sen_Salubrinal, x="Risk",
               y = "IC50", color = "Risk",
               palette = "jco", add = "jitter")
# 添加p值
p + stat_compare_means()
dev.off()

####################################drug14： Lapatinib################################################
drug_sen_Lapatinib<-cbind(merged_format$MC,as.numeric(merged_format$Lapatinib))
rownames(drug_sen_Lapatinib)<-row.names(merged_format)
colnames(drug_sen_Lapatinib)<-c("Risk","IC50")
head(drug_sen_Lapatinib)
drug_sen_Lapatinib<-as.data.frame(drug_sen_Lapatinib)
drug_sen_Lapatinib$IC50<-as.numeric(drug_sen_Lapatinib$IC50)
library(ggpubr)
pdf("Lapatinib.pdf")
p <- ggboxplot(drug_sen_Lapatinib, x="Risk",
               y = "IC50", color = "Risk",
               palette = "jco", add = "jitter")
# 添加p值
p + stat_compare_means()
dev.off()

####################################drug15： Doxorubicin################################################
drug_sen_Doxorubicin<-cbind(merged_format$MC,as.numeric(merged_format$Doxorubicin))
rownames(drug_sen_Doxorubicin)<-row.names(merged_format)
colnames(drug_sen_Doxorubicin)<-c("Risk","IC50")
head(drug_sen_Doxorubicin)
drug_sen_Doxorubicin<-as.data.frame(drug_sen_Doxorubicin)
drug_sen_Doxorubicin$IC50<-as.numeric(drug_sen_Doxorubicin$IC50)
library(ggpubr)
pdf("Doxorubicin.pdf")
p <- ggboxplot(drug_sen_Doxorubicin, x="Risk",
               y = "IC50", color = "Risk",
               palette = "jco", add = "jitter")
# 添加p值
p + stat_compare_means()
dev.off()

####################################drug16： Etoposide################################################
drug_sen_Etoposide<-cbind(merged_format$MC,as.numeric(merged_format$Etoposide))
rownames(drug_sen_Etoposide)<-row.names(merged_format)
colnames(drug_sen_Etoposide)<-c("Risk","IC50")
head(drug_sen_Etoposide)
drug_sen_Etoposide<-as.data.frame(drug_sen_Etoposide)
drug_sen_Etoposide$IC50<-as.numeric(drug_sen_Etoposide$IC50)
library(ggpubr)
pdf("Etoposide.pdf")
p <- ggboxplot(drug_sen_Etoposide, x="Risk",
               y = "IC50", color = "Risk",
               palette = "jco", add = "jitter")
# 添加p值
p + stat_compare_means()
dev.off()

####################################drug17： Gemcitabine################################################
drug_sen_Gemcitabine<-cbind(merged_format$MC,as.numeric(merged_format$Gemcitabine))
rownames(drug_sen_Gemcitabine)<-row.names(merged_format)
colnames(drug_sen_Gemcitabine)<-c("Risk","IC50")
head(drug_sen_Gemcitabine)
drug_sen_Gemcitabine<-as.data.frame(drug_sen_Gemcitabine)
drug_sen_Gemcitabine$IC50<-as.numeric(drug_sen_Gemcitabine$IC50)
library(ggpubr)
pdf("Gemcitabine.pdf")
p <- ggboxplot(drug_sen_Gemcitabine, x="Risk",
               y = "IC50", color = "Risk",
               palette = "jco", add = "jitter")
# 添加p值
p + stat_compare_means()
dev.off()

####################################drug18： Vinorelbine################################################
drug_sen_Vinorelbine<-cbind(merged_format$MC,as.numeric(merged_format$Vinorelbine))
rownames(drug_sen_Vinorelbine)<-row.names(merged_format)
colnames(drug_sen_Vinorelbine)<-c("Risk","IC50")
head(drug_sen_Vinorelbine)
drug_sen_Vinorelbine<-as.data.frame(drug_sen_Vinorelbine)
drug_sen_Vinorelbine$IC50<-as.numeric(drug_sen_Vinorelbine$IC50)
library(ggpubr)
pdf("Vinorelbine.pdf")
p <- ggboxplot(drug_sen_Vinorelbine, x="Risk",
               y = "IC50", color = "Risk",
               palette = "jco", add = "jitter")
# 添加p值
p + stat_compare_means()
dev.off()

####################################drug19： Bicalutamide################################################
drug_sen_Bicalutamide<-cbind(merged_format$MC,as.numeric(merged_format$Bicalutamide))
rownames(drug_sen_Bicalutamide)<-row.names(merged_format)
colnames(drug_sen_Bicalutamide)<-c("Risk","IC50")
head(drug_sen_Bicalutamide)
drug_sen_Bicalutamide<-as.data.frame(drug_sen_Bicalutamide)
drug_sen_Bicalutamide$IC50<-as.numeric(drug_sen_Bicalutamide$IC50)
library(ggpubr)
pdf("Bicalutamide.pdf")
p <- ggboxplot(drug_sen_Bicalutamide, x="Risk",
               y = "IC50", color = "Risk",
               palette = "jco", add = "jitter")
# 添加p值
p + stat_compare_means()
dev.off()

####################################drug20： QS11################################################
drug_sen_QS11<-cbind(merged_format$MC,as.numeric(merged_format$QS11))
rownames(drug_sen_QS11)<-row.names(merged_format)
colnames(drug_sen_QS11)<-c("Risk","IC50")
head(drug_sen_QS11)
drug_sen_QS11<-as.data.frame(drug_sen_QS11)
drug_sen_QS11$IC50<-as.numeric(drug_sen_QS11$IC50)
library(ggpubr)
pdf("QS11.pdf")
p <- ggboxplot(drug_sen_QS11, x="Risk",
               y = "IC50", color = "Risk",
               palette = "jco", add = "jitter")
# 添加p值
p + stat_compare_means()
dev.off()

####################################drug21： Midostaurin################################################
drug_sen_Midostaurin<-cbind(merged_format$MC,as.numeric(merged_format$Midostaurin))
rownames(drug_sen_Midostaurin)<-row.names(merged_format)
colnames(drug_sen_Midostaurin)<-c("Risk","IC50")
head(drug_sen_Midostaurin)
drug_sen_Midostaurin<-as.data.frame(drug_sen_Midostaurin)
drug_sen_Midostaurin$IC50<-as.numeric(drug_sen_Midostaurin$IC50)
library(ggpubr)
pdf("Midostaurin.pdf")
p <- ggboxplot(drug_sen_Midostaurin, x="Risk",
               y = "IC50", color = "Risk",
               palette = "jco", add = "jitter")
# 添加p值
p + stat_compare_means()
dev.off()

####################################drug22： AZD6482################################################
drug_sen_AZD6482<-cbind(merged_format$MC,as.numeric(merged_format$AZD6482))
rownames(drug_sen_AZD6482)<-row.names(merged_format)
colnames(drug_sen_AZD6482)<-c("Risk","IC50")
head(drug_sen_AZD6482)
drug_sen_AZD6482<-as.data.frame(drug_sen_AZD6482)
drug_sen_AZD6482$IC50<-as.numeric(drug_sen_AZD6482$IC50)
library(ggpubr)
pdf("AZD6482.pdf")
p <- ggboxplot(drug_sen_AZD6482, x="Risk",
               y = "IC50", color = "Risk",
               palette = "jco", add = "jitter")
# 添加p值
p + stat_compare_means()
dev.off()

####################################drug23： DMOG################################################
drug_sen_DMOG<-cbind(merged_format$MC,as.numeric(merged_format$DMOG))
rownames(drug_sen_DMOG)<-row.names(merged_format)
colnames(drug_sen_DMOG)<-c("Risk","IC50")
head(drug_sen_DMOG)
drug_sen_DMOG<-as.data.frame(drug_sen_DMOG)
drug_sen_DMOG$IC50<-as.numeric(drug_sen_DMOG$IC50)
library(ggpubr)
pdf("DMOG.pdf")
p <- ggboxplot(drug_sen_DMOG, x="Risk",
               y = "IC50", color = "Risk",
               palette = "jco", add = "jitter")
# 添加p值
p + stat_compare_means()
dev.off()

####################################drug24： Embelin################################################
drug_sen_Embelin<-cbind(merged_format$MC,as.numeric(merged_format$Embelin))
rownames(drug_sen_Embelin)<-row.names(merged_format)
colnames(drug_sen_Embelin)<-c("Risk","IC50")
head(drug_sen_Embelin)
drug_sen_Embelin<-as.data.frame(drug_sen_Embelin)
drug_sen_Embelin$IC50<-as.numeric(drug_sen_Embelin$IC50)
library(ggpubr)
pdf("Embelin.pdf")
p <- ggboxplot(drug_sen_Embelin, x="Risk",
               y = "IC50", color = "Risk",
               palette = "jco", add = "jitter")
# 添加p值
p + stat_compare_means()
dev.off()

####################################drug25： FH535################################################
drug_sen_FH535<-cbind(merged_format$MC,as.numeric(merged_format$FH535))
rownames(drug_sen_FH535)<-row.names(merged_format)
colnames(drug_sen_FH535)<-c("Risk","IC50")
head(drug_sen_FH535)
drug_sen_FH535<-as.data.frame(drug_sen_FH535)
drug_sen_FH535$IC50<-as.numeric(drug_sen_FH535$IC50)
library(ggpubr)
pdf("FH535.pdf")
p <- ggboxplot(drug_sen_FH535, x="Risk",
               y = "IC50", color = "Risk",
               palette = "jco", add = "jitter")
# 添加p值
p + stat_compare_means()
dev.off()

####################################drug26： Thapsigargin################################################
drug_sen_Thapsigargin<-cbind(merged_format$MC,as.numeric(merged_format$Thapsigargin))
rownames(drug_sen_Thapsigargin)<-row.names(merged_format)
colnames(drug_sen_Thapsigargin)<-c("Risk","IC50")
head(drug_sen_Thapsigargin)
drug_sen_Thapsigargin<-as.data.frame(drug_sen_Thapsigargin)
drug_sen_Thapsigargin$IC50<-as.numeric(drug_sen_Thapsigargin$IC50)
library(ggpubr)
pdf("Thapsigargin.pdf")
p <- ggboxplot(drug_sen_Thapsigargin, x="Risk",
               y = "IC50", color = "Risk",
               palette = "jco", add = "jitter")
# 添加p值
p + stat_compare_means()
dev.off()

####################################drug27： Bexarotene################################################
drug_sen_Bexarotene<-cbind(merged_format$MC,as.numeric(merged_format$Bexarotene))
rownames(drug_sen_Bexarotene)<-row.names(merged_format)
colnames(drug_sen_Bexarotene)<-c("Risk","IC50")
head(drug_sen_Bexarotene)
drug_sen_Bexarotene<-as.data.frame(drug_sen_Bexarotene)
drug_sen_Bexarotene$IC50<-as.numeric(drug_sen_Bexarotene$IC50)
library(ggpubr)
pdf("Bexarotene.pdf")
p <- ggboxplot(drug_sen_Bexarotene, x="Risk",
               y = "IC50", color = "Risk",
               palette = "jco", add = "jitter")
# 添加p值
p + stat_compare_means()
dev.off()

####################################drug28： Bleomycin################################################
drug_sen_Bleomycin<-cbind(merged_format$MC,as.numeric(merged_format$Bleomycin))
rownames(drug_sen_Bleomycin)<-row.names(merged_format)
colnames(drug_sen_Bleomycin)<-c("Risk","IC50")
head(drug_sen_Bleomycin)
drug_sen_Bleomycin<-as.data.frame(drug_sen_Bleomycin)
drug_sen_Bleomycin$IC50<-as.numeric(drug_sen_Bleomycin$IC50)
library(ggpubr)
pdf("Bleomycin.pdf")
p <- ggboxplot(drug_sen_Bleomycin, x="Risk",
               y = "IC50", color = "Risk",
               palette = "jco", add = "jitter")
# 添加p值
p + stat_compare_means()
dev.off()

####################################drug29： Pazopanib################################################
drug_sen_Pazopanib<-cbind(merged_format$MC,as.numeric(merged_format$Pazopanib))
rownames(drug_sen_Pazopanib)<-row.names(merged_format)
colnames(drug_sen_Pazopanib)<-c("Risk","IC50")
head(drug_sen_Pazopanib)
drug_sen_Pazopanib<-as.data.frame(drug_sen_Pazopanib)
drug_sen_Pazopanib$IC50<-as.numeric(drug_sen_Pazopanib$IC50)
library(ggpubr)
pdf("Pazopanib.pdf")
p <- ggboxplot(drug_sen_Pazopanib, x="Risk",
               y = "IC50", color = "Risk",
               palette = "jco", add = "jitter")
# 添加p值
p + stat_compare_means()
dev.off()

####################################drug30： Tipifarnib################################################
drug_sen_Tipifarnib<-cbind(merged_format$MC,as.numeric(merged_format$Tipifarnib))
rownames(drug_sen_Tipifarnib)<-row.names(merged_format)
colnames(drug_sen_Tipifarnib)<-c("Risk","IC50")
head(drug_sen_Tipifarnib)
drug_sen_Tipifarnib<-as.data.frame(drug_sen_Tipifarnib)
drug_sen_Tipifarnib$IC50<-as.numeric(drug_sen_Tipifarnib$IC50)
library(ggpubr)
pdf("Tipifarnib.pdf")
p <- ggboxplot(drug_sen_Tipifarnib, x="Risk",
               y = "IC50", color = "Risk",
               palette = "jco", add = "jitter")
# 添加p值
p + stat_compare_means()
dev.off()

####################################drug31： AS601245################################################
drug_sen_AS601245<-cbind(merged_format$MC,as.numeric(merged_format$AS601245))
rownames(drug_sen_AS601245)<-row.names(merged_format)
colnames(drug_sen_AS601245)<-c("Risk","IC50")
head(drug_sen_AS601245)
drug_sen_AS601245<-as.data.frame(drug_sen_AS601245)
drug_sen_AS601245$IC50<-as.numeric(drug_sen_AS601245$IC50)
library(ggpubr)
pdf("AS601245.pdf")
p <- ggboxplot(drug_sen_AS601245, x="Risk",
               y = "IC50", color = "Risk",
               palette = "jco", add = "jitter")
# 添加p值
p + stat_compare_means()
dev.off()

####################################drug32： AICAR################################################
drug_sen_AICAR<-cbind(merged_format$MC,as.numeric(merged_format$AICAR))
rownames(drug_sen_AICAR)<-row.names(merged_format)
colnames(drug_sen_AICAR)<-c("Risk","IC50")
head(drug_sen_AICAR)
drug_sen_AICAR<-as.data.frame(drug_sen_AICAR)
drug_sen_AICAR$IC50<-as.numeric(drug_sen_AICAR$IC50)
library(ggpubr)
pdf("AICAR.pdf")
p <- ggboxplot(drug_sen_AICAR, x="Risk",
               y = "IC50", color = "Risk",
               palette = "jco", add = "jitter")
# 添加p值
p + stat_compare_means()
dev.off()

####################################drug33： Camptothecin################################################
drug_sen_Camptothecin<-cbind(merged_format$MC,as.numeric(merged_format$Camptothecin))
rownames(drug_sen_Camptothecin)<-row.names(merged_format)
colnames(drug_sen_Camptothecin)<-c("Risk","IC50")
head(drug_sen_Camptothecin)
drug_sen_Camptothecin<-as.data.frame(drug_sen_Camptothecin)
drug_sen_Camptothecin$IC50<-as.numeric(drug_sen_Camptothecin$IC50)
library(ggpubr)
pdf("Camptothecin.pdf")
p <- ggboxplot(drug_sen_Camptothecin, x="Risk",
               y = "IC50", color = "Risk",
               palette = "jco", add = "jitter")
# 添加p值
p + stat_compare_means()
dev.off()

####################################drug34： Vinblastine################################################
drug_sen_Vinblastine<-cbind(merged_format$MC,as.numeric(merged_format$Vinblastine))
rownames(drug_sen_Vinblastine)<-row.names(merged_format)
colnames(drug_sen_Vinblastine)<-c("Risk","IC50")
head(drug_sen_Vinblastine)
drug_sen_Vinblastine<-as.data.frame(drug_sen_Vinblastine)
drug_sen_Vinblastine$IC50<-as.numeric(drug_sen_Vinblastine$IC50)
library(ggpubr)
pdf("Vinblastine.pdf")
p <- ggboxplot(drug_sen_Vinblastine, x="Risk",
               y = "IC50", color = "Risk",
               palette = "jco", add = "jitter")
# 添加p值
p + stat_compare_means()
dev.off()

####################################drug35： Cisplatin################################################
drug_sen_Cisplatin<-cbind(merged_format$MC,as.numeric(merged_format$Cisplatin))
rownames(drug_sen_Cisplatin)<-row.names(merged_format)
colnames(drug_sen_Cisplatin)<-c("Risk","IC50")
head(drug_sen_Cisplatin)
drug_sen_Cisplatin<-as.data.frame(drug_sen_Cisplatin)
drug_sen_Cisplatin$IC50<-as.numeric(drug_sen_Cisplatin$IC50)
library(ggpubr)
pdf("Cisplatin.pdf")
p <- ggboxplot(drug_sen_Cisplatin, x="Risk",
               y = "IC50", color = "Risk",
               palette = "jco", add = "jitter")
# 添加p值
p + stat_compare_means()
dev.off()

####################################drug36： Cytarabine################################################
drug_sen_Cytarabine<-cbind(merged_format$MC,as.numeric(merged_format$Cytarabine))
rownames(drug_sen_Cytarabine)<-row.names(merged_format)
colnames(drug_sen_Cytarabine)<-c("Risk","IC50")
head(drug_sen_Cytarabine)
drug_sen_Cytarabine<-as.data.frame(drug_sen_Cytarabine)
drug_sen_Cytarabine$IC50<-as.numeric(drug_sen_Cytarabine$IC50)
library(ggpubr)
pdf("Cytarabine.pdf")
p <- ggboxplot(drug_sen_Cytarabine, x="Risk",
               y = "IC50", color = "Risk",
               palette = "jco", add = "jitter")
# 添加p值
p + stat_compare_means()
dev.off()

####################################drug37： Docetaxel################################################
drug_sen_Docetaxel<-cbind(merged_format$MC,as.numeric(merged_format$Docetaxel))
rownames(drug_sen_Docetaxel)<-row.names(merged_format)
colnames(drug_sen_Docetaxel)<-c("Risk","IC50")
head(drug_sen_Docetaxel)
drug_sen_Docetaxel<-as.data.frame(drug_sen_Docetaxel)
drug_sen_Docetaxel$IC50<-as.numeric(drug_sen_Docetaxel$IC50)
library(ggpubr)
pdf("Docetaxel.pdf")
p <- ggboxplot(drug_sen_Docetaxel, x="Risk",
               y = "IC50", color = "Risk",
               palette = "jco", add = "jitter")
# 添加p值
p + stat_compare_means()
dev.off()

####################################drug38： Methotrexate################################################
drug_sen_Methotrexate<-cbind(merged_format$MC,as.numeric(merged_format$Methotrexate))
rownames(drug_sen_Methotrexate)<-row.names(merged_format)
colnames(drug_sen_Methotrexate)<-c("Risk","IC50")
head(drug_sen_Methotrexate)
drug_sen_Methotrexate<-as.data.frame(drug_sen_Methotrexate)
drug_sen_Methotrexate$IC50<-as.numeric(drug_sen_Methotrexate$IC50)
library(ggpubr)
pdf("Methotrexate.pdf")
p <- ggboxplot(drug_sen_Methotrexate, x="Risk",
               y = "IC50", color = "Risk",
               palette = "jco", add = "jitter")
# 添加p值
p + stat_compare_means()
dev.off()

####################################drug39： ATRA################################################
drug_sen_ATRA<-cbind(merged_format$MC,as.numeric(merged_format$ATRA))
rownames(drug_sen_ATRA)<-row.names(merged_format)
colnames(drug_sen_ATRA)<-c("Risk","IC50")
head(drug_sen_ATRA)
drug_sen_ATRA<-as.data.frame(drug_sen_ATRA)
drug_sen_ATRA$IC50<-as.numeric(drug_sen_ATRA$IC50)
library(ggpubr)
pdf("ATRA.pdf")
p <- ggboxplot(drug_sen_ATRA, x="Risk",
               y = "IC50", color = "Risk",
               palette = "jco", add = "jitter")
# 添加p值
p + stat_compare_means()
dev.off()

####################################drug40： Vorinostat################################################
drug_sen_Vorinostat<-cbind(merged_format$MC,as.numeric(merged_format$Vorinostat))
rownames(drug_sen_Vorinostat)<-row.names(merged_format)
colnames(drug_sen_Vorinostat)<-c("Risk","IC50")
head(drug_sen_Vorinostat)
drug_sen_Vorinostat<-as.data.frame(drug_sen_Vorinostat)
drug_sen_Vorinostat$IC50<-as.numeric(drug_sen_Vorinostat$IC50)
library(ggpubr)
pdf("Vorinostat.pdf")
p <- ggboxplot(drug_sen_Vorinostat, x="Risk",
               y = "IC50", color = "Risk",
               palette = "jco", add = "jitter")
# 添加p值
p + stat_compare_means()
dev.off()

####################################drug41： Temsirolimus################################################
drug_sen_Temsirolimus<-cbind(merged_format$MC,as.numeric(merged_format$Temsirolimus))
rownames(drug_sen_Temsirolimus)<-row.names(merged_format)
colnames(drug_sen_Temsirolimus)<-c("Risk","IC50")
head(drug_sen_Temsirolimus)
drug_sen_Temsirolimus<-as.data.frame(drug_sen_Temsirolimus)
drug_sen_Temsirolimus$IC50<-as.numeric(drug_sen_Temsirolimus$IC50)
library(ggpubr)
pdf("Temsirolimus.pdf")
p <- ggboxplot(drug_sen_Temsirolimus, x="Risk",
               y = "IC50", color = "Risk",
               palette = "jco", add = "jitter")
# 添加p值
p + stat_compare_means()
dev.off()

####################################drug42： Axitinib################################################
drug_sen_Axitinib<-cbind(merged_format$MC,as.numeric(merged_format$Axitinib))
rownames(drug_sen_Axitinib)<-row.names(merged_format)
colnames(drug_sen_Axitinib)<-c("Risk","IC50")
head(drug_sen_Axitinib)
drug_sen_Axitinib<-as.data.frame(drug_sen_Axitinib)
drug_sen_Axitinib$IC50<-as.numeric(drug_sen_Axitinib$IC50)
library(ggpubr)
pdf("Axitinib.pdf")
p <- ggboxplot(drug_sen_Axitinib, x="Risk",
               y = "IC50", color = "Risk",
               palette = "jco", add = "jitter")
# 添加p值
p + stat_compare_means()
dev.off()

####################################drug43： PLX4720################################################
drug_sen_PLX4720<-cbind(merged_format$MC,as.numeric(merged_format$PLX4720))
rownames(drug_sen_PLX4720)<-row.names(merged_format)
colnames(drug_sen_PLX4720)<-c("Risk","IC50")
head(drug_sen_PLX4720)
drug_sen_PLX4720<-as.data.frame(drug_sen_PLX4720)
drug_sen_PLX4720$IC50<-as.numeric(drug_sen_PLX4720$IC50)
library(ggpubr)
pdf("PLX4720.pdf")
p <- ggboxplot(drug_sen_PLX4720, x="Risk",
               y = "IC50", color = "Risk",
               palette = "jco", add = "jitter")
# 添加p值
p + stat_compare_means()
dev.off()

####################################drug44： PLX4720################################################
drug_sen_PLX4720<-cbind(merged_format$MC,as.numeric(merged_format$PLX4720))
rownames(drug_sen_PLX4720)<-row.names(merged_format)
colnames(drug_sen_PLX4720)<-c("Risk","IC50")
head(drug_sen_PLX4720)
drug_sen_PLX4720<-as.data.frame(drug_sen_PLX4720)
drug_sen_PLX4720$IC50<-as.numeric(drug_sen_PLX4720$IC50)
library(ggpubr)
pdf("PLX4720.pdf")
p <- ggboxplot(drug_sen_PLX4720, x="Risk",
               y = "IC50", color = "Risk",
               palette = "jco", add = "jitter")
# 添加p值
p + stat_compare_means()
dev.off()

####################################drug45： GDC0941################################################
drug_sen_GDC0941<-cbind(merged_format$MC,as.numeric(merged_format$GDC0941))
rownames(drug_sen_GDC0941)<-row.names(merged_format)
colnames(drug_sen_GDC0941)<-c("Risk","IC50")
head(drug_sen_GDC0941)
drug_sen_GDC0941<-as.data.frame(drug_sen_GDC0941)
drug_sen_GDC0941$IC50<-as.numeric(drug_sen_GDC0941$IC50)
library(ggpubr)
pdf("GDC0941.pdf")
p <- ggboxplot(drug_sen_GDC0941, x="Risk",
               y = "IC50", color = "Risk",
               palette = "jco", add = "jitter")
# 添加p值
p + stat_compare_means()
dev.off()

####################################drug46： CCT007093################################################
drug_sen_CCT007093<-cbind(merged_format$MC,as.numeric(merged_format$CCT007093))
rownames(drug_sen_CCT007093)<-row.names(merged_format)
colnames(drug_sen_CCT007093)<-c("Risk","IC50")
head(drug_sen_CCT007093)
drug_sen_CCT007093<-as.data.frame(drug_sen_CCT007093)
drug_sen_CCT007093$IC50<-as.numeric(drug_sen_CCT007093$IC50)
library(ggpubr)
pdf("CCT007093.pdf")
p <- ggboxplot(drug_sen_CCT007093, x="Risk",
               y = "IC50", color = "Risk",
               palette = "jco", add = "jitter")
# 添加p值
p + stat_compare_means()
dev.off()

####################################drug47： CCT018159################################################
drug_sen_CCT018159<-cbind(merged_format$MC,as.numeric(merged_format$CCT018159))
rownames(drug_sen_CCT018159)<-row.names(merged_format)
colnames(drug_sen_CCT018159)<-c("Risk","IC50")
head(drug_sen_CCT018159)
drug_sen_CCT018159<-as.data.frame(drug_sen_CCT018159)
drug_sen_CCT018159$IC50<-as.numeric(drug_sen_CCT018159$IC50)
library(ggpubr)
pdf("CCT018159.pdf")
p <- ggboxplot(drug_sen_CCT018159, x="Risk",
               y = "IC50", color = "Risk",
               palette = "jco", add = "jitter")
# 添加p值
p + stat_compare_means()
dev.off()


colnames(merged_format)
#[1] "MC"            "Erlotinib"     "Rapamycin"     "Sunitinib"     "Paclitaxel"    "Cyclopamine"  
#[7] "AZ628"         "Sorafenib"     "Imatinib"      "Dasatinib"     "CMK"           "Pyrimethamine"
#[13] "GW843682X"     "Parthenolide"  "Bortezomib"    "Roscovitine"   "Salubrinal"    "Lapatinib"    
#[19] "GSK269962A"    "Doxorubicin"   "Etoposide"     "Gemcitabine"   "Vinorelbine"   "Bicalutamide" 
#[25] "QS11"          "Midostaurin"   "AZD6482"       "DMOG"          "Shikonin"      "Embelin"      
#[31] "FH535"         "Thapsigargin"  "Bexarotene"    "Bleomycin"     "AUY922"        "Pazopanib"    
#[37] "Tipifarnib"    "AS601245"      "AICAR"         "Camptothecin"  "Vinblastine"   "Cisplatin"    
#[43] "Cytarabine"    "Docetaxel"     "Methotrexate"  "ATRA"          "Gefitinib"     "Vorinostat"   
#[49] "Nilotinib"     "RDEA119"       "Temsirolimus"  "Bosutinib"     "Lenalidomide"  "Axitinib"     
#[55] "Elesclomol"    "PLX4720"       "GDC0941"       "AZD8055"       "SB590885"      "CCT007093"    
#[61] "CCT018159" 


