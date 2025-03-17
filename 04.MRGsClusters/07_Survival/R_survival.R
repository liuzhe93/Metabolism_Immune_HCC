setwd("F:/Ashley_Projects/Metabolism_Immune/04.MRGsClusters/07_Survival")
rm(list=ls())

Sample_cluster<-read.csv("F:/Ashley_Projects/Metabolism_Immune/04.MRGsClusters/05_Heatmap/SampleID_MetabolicScore.csv")
class(Sample_cluster)
head(Sample_cluster)
sam_clu_sel<-subset(Sample_cluster, select = c("X", "cluster"))
dim(sam_clu_sel)
#[1] 390   2

sur_time<-read.table("F:/Ashley_Projects/Metabolism_Immune/DataCollection/clinical_rmNA.txt",header=T)
head(sur_time)

sam_clu_sel$Id<-substring(sam_clu_sel$X, 1, 12)
head(sam_clu_sel)

merged_data<-merge(sur_time, sam_clu_sel, by = "Id")
dim(sur_time)
#[1] 349   9
dim(sam_clu_sel)
#[1] 390   3
dim(merged_data)
#[1] 355  11

merged_data_sel<-subset(merged_data, select = c("Id", "futime", "fustat", "cluster"))
#Calculate the optimal cutoff
library(survival)
library(survminer)


res.cat<-merged_data_sel[, -1]

fit=survfit(Surv(futime, fustat) ~cluster, data = res.cat)

head(res.cat)
diff<-survdiff(Surv(futime, fustat)~cluster, data = res.cat)
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
                   legend.labs=c("3" ,"2", "1"),
                   legend.title="Cluster",
                   xlab="Time (years)",
                   break.time.by = 1,
                   risk.table.title="",
                   palette=c("green", "red", "blue"),
                   risk.table=T,
                   risk.table.height=.25)
outFile="Cluster_survival.pdf"
pdf(file=outFile,onefile = FALSE,width = 8,height =5)
print(surPlot)
dev.off()

