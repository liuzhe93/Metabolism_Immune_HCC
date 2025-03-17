setwd("F:/Ashley_Projects/Metabolism_Immune/11.Immunotherapy/02.Drug/RiskScoreDrugIC50")
rm(list=ls())


ic50<-read.csv("F:/Ashley_Projects/Metabolism_Immune/11.Immunotherapy/02.Drug/drugSensitivity.csv")
riskScore<-read.table("F:/Ashley_Projects/Metabolism_Immune/06.MR_relatedRiskModel/01_RiskModel/risk.txt")

ic50_sel<-ic50[, c("X", "Sorafenib")]
ic50_sel$ID<-substr(ic50_sel$X, 1, 12)
ic50_sel$ID<-gsub("\\.", "-", ic50_sel$ID)
  
riskScore$SampleName<-row.names(riskScore)
riskScore_sel<-riskScore[, c("riskScore", "SampleName")]
riskScore_sel$ID<-substring(riskScore_sel$SampleName, 1, 12)
riskScore_sel$SampleName<-NULL
head(riskScore_sel)

merged_data<-merge(ic50_sel, riskScore_sel, by = "ID")
for(i in 1:nrow(merged_data)){
  temp<-paste(merged_data$ID[i], i, sep = "_")
  row.names(merged_data)[i]<-temp
}
merged_format<-merged_data[, c("Sorafenib", "riskScore")]

drug_sen_Sorafenib<-merged_format
colnames(drug_sen_Sorafenib)<-c("IC50", "Risk")
head(drug_sen_Sorafenib)
drug_sen_Sorafenib<-as.data.frame(drug_sen_Sorafenib)
drug_sen_Sorafenib$IC50<-as.numeric(drug_sen_Sorafenib$IC50)



pdf("CAF_riskscore.pdf")
ggscatter(drug_sen_Sorafenib, x = "Risk", y = "IC50", add = "reg.line", conf.int = TRUE, 
          add.params = list(full = "lightgrey"))+
  stat_cor(method = "pearson")
dev.off()



