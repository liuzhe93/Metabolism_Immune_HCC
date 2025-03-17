setwd("F:/Ashley_Projects/Metabolism_Immune/13.PON1gene/02.external/")
rm(list=ls())

geneExp<-read.table("F:/Ashley_Projects/Metabolism_Immune/06.MR_relatedRiskModel/02_ExternalDataSet_GSE124535/dataset/GSE124535_HCC.RNA-seq.35.samples.fpkm.txt", header = 1, row.names = 1)
dim(geneExp)
#[1] 63652    74
geneExp$Chr<-NULL
geneExp$Biotype<-NULL
attach(geneExp)
geneExp_uniq<-aggregate(geneExp, by = list(gene_symbol), FUN = mean)
geneExp_uniq$gene_symbol<-NULL
geneExp_uniq[1:28,1]<-c("DELEC1", "MTARC1", "SEPT1", "MARCHF10", "SEPTIN10", 
                        "MARCHF11", "SEPTIN11", "SEPTIN12", "SEPTIN14", "SELENOF",
                        "MARCHF2", "SEPTIN2", "MARCHF3", "SEPTIN3", "MARCHF4",
                        "SEPTIN4", "MARCHF5", "SEPTIN5", "5S_rRNA", "MARCHF6", 
                        "SEPTIN6", "MARCHF7", "SEPTIN7", "7SK", "MARCH8",
                        "SEPTIN8", "MARCHF9", "SEPTIN9")
geneExp_uniq[1:6,1:4]

library("openxlsx")
clicalInfo<-read.xlsx("F:/Ashley_Projects/Metabolism_Immune/06.MR_relatedRiskModel/02_ExternalDataSet_GSE124535/dataset/2017-04-04932E-Supplementary Table 1.xlsx", startRow = 2, colNames = TRUE)
dim(clicalInfo)
clicalInfo[1:6,1:4]
clinicalInfo_sel<-clicalInfo[clicalInfo$`RNA-Seqb`=="1",]
dim(clinicalInfo_sel)
dim(clicalInfo)
head(clinicalInfo_sel)
clinicalInfo_sel<-subset(clinicalInfo_sel, select = c("Case.No.", "Gender", "Age", "Cancer.recurrenc.e.b",
                                                      "Disease.free.survival.(m)"))
colnames(clinicalInfo_sel)<-c("Patient", "Gender", "Age", "CancerRecurrence", "DiseaseFreeSurvival_m")
head(clinicalInfo_sel)
row.names(geneExp_uniq)<-geneExp_uniq$Group.1
geneExp_uniq$Group.1<-NULL
geneExp_t<-t(geneExp_uniq)
dim(geneExp_t)
geneExp_t[1:6,1:5]
geneExp_t<-as.data.frame(geneExp_t)
geneExp_t$SampleID<-row.names(geneExp_t)
dim(geneExp_t)
geneExp_t$SampleID<-substr(geneExp_t$SampleID, 1, 4)
geneExp_t[1:6, 56610:56616]
merged_data<-merge(x = geneExp_t, y = clinicalInfo_sel, by.x = "SampleID", by.y = "Patient")
dim(merged_data)
merged_data[1:6, c(1:8,56614:56620)]
merged_data_sel<-merged_data[, c("SampleID","CancerRecurrence", "PON1", "DiseaseFreeSurvival_m")]
for(i in 1:nrow(merged_data_sel)){
  merged_data_sel$id[i]<-paste(merged_data_sel$SampleID[i], i, sep = "_")
}
head(merged_data_sel)
rownames(merged_data_sel)<-merged_data_sel$id
merged_data_sel$id<-NULL
merged_data_sel$SampleID<-NULL
rt<-cbind(merged_data_sel$DiseaseFreeSurvival_m, merged_data_sel$CancerRecurrence, merged_data_sel$PON1)
row.names(rt)<-row.names(merged_data_sel)
colnames(rt)<-c("futime", "fustat", "PON1")
head(rt)
rt<-as.data.frame(rt)

#Calculate the optimal cutoff
library(survival)
library(survminer)
res.cut<-surv_cutpoint(rt,time="futime",event="fustat",variables="PON1")
summary(res.cut)
#     cutpoint statistic
#PON1  373.898  1.829868

pdf("Cluster_cutoff.pdf")
plot(res.cut,"PON1",palette="npg")
dev.off()

res.cat<-surv_categorize(res.cut)

fit=survfit(Surv(futime, fustat) ~PON1, data = res.cat)

head(res.cat)
diff<-survdiff(Surv(futime, fustat)~PON1, data = res.cat)
pValue<-1-pchisq(diff$chisq, df=1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
surPlot=ggsurvplot(fit, 
                   data=res.cat,
                   conf.int=TRUE,
                   pval=pValue,
                   pval.size=5,
                   legend.labs=c("High", "Low"),
                   legend.title="PON1",
                   xlab="Time (months)",
                   break.time.by = 10,
                   risk.table.title="",
                   palette=c("red", "blue"),
                   risk.table=T,
                   risk.table.height=.25)
outFile="Cluster_survival.pdf"
pdf(file=outFile,onefile = FALSE,width = 8,height =5)
print(surPlot)
dev.off()


