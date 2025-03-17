setwd("F:/Ashley_Projects/Metabolism_Immune/11.Immunotherapy/01.TIDE/RiskScore")
rm(list=ls())


riskScore<-read.table("F:/Ashley_Projects/Metabolism_Immune/06.MR_relatedRiskModel/01_RiskModel/risk.txt")
head(riskScore)
MC_result <- read.csv('F:/Ashley_Projects/Metabolism_Immune/11.Immunotherapy/01.TIDE/MC/TIDE_res.csv')
head(MC_result)

riskScore$ID<-row.names(riskScore)
library("stringr")
riskScore$ID<-str_sub(riskScore$ID, 1, 12)

MC_result$ID<-MC_result$Patient
MC_result$ID<-str_sub(MC_result$ID, 4, 15)

merged_data<-merge(riskScore, MC_result, by = "ID")
mydata<-subset(merged_data, select = c("ID", "riskScore", "CAF", "Dysfunction", 
                                       "Exclusion", "MDSC", "TAM.M2", "TIDE"))

pdf("CAF_riskscore.pdf")
ggscatter(mydata, x = "riskScore", y = "CAF", add = "reg.line", conf.int = TRUE, 
          add.params = list(full = "lightgrey"))+
  stat_cor(method = "pearson")
dev.off()

pdf("Dysfunction_riskscore.pdf")
ggscatter(mydata, x = "riskScore", y = "Dysfunction", add = "reg.line", conf.int = TRUE, 
          add.params = list(full = "lightgrey"))+
  stat_cor(method = "pearson")
dev.off()

pdf("Exclusion_riskscore.pdf")
ggscatter(mydata, x = "riskScore", y = "Exclusion", add = "reg.line", conf.int = TRUE, 
          add.params = list(full = "lightgrey"))+
  stat_cor(method = "pearson")
dev.off()

pdf("MDSC_riskscore.pdf")
ggscatter(mydata, x = "riskScore", y = "MDSC", add = "reg.line", conf.int = TRUE, 
          add.params = list(full = "lightgrey"))+
  stat_cor(method = "pearson")
dev.off()

pdf("TAM.M2_riskscore.pdf")
ggscatter(mydata, x = "riskScore", y = "TAM.M2", add = "reg.line", conf.int = TRUE, 
          add.params = list(full = "lightgrey"))+
  stat_cor(method = "pearson")
dev.off()

pdf("TIDE_riskscore.pdf")
ggscatter(mydata, x = "riskScore", y = "TIDE", add = "reg.line", conf.int = TRUE, 
          add.params = list(full = "lightgrey"))+
  stat_cor(method = "pearson")
dev.off()


