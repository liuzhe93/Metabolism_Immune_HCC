setwd("F:/Ashley_Projects/Metabolism_Immune/03.ClusterTME/05_survival")
rm(list=ls())

Sample_cluster<-read.csv("F:/Ashley_Projects/Metabolism_Immune/03.ClusterTME/02_heatmap/SampleID_MetabolicScore.csv")
class(Sample_cluster)
head(Sample_cluster)
sam_clu_sel<-subset(Sample_cluster, select = c("X", "sample_subtypes"))
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

merged_data_sel<-subset(merged_data, select = c("Id", "futime", "fustat", "sample_subtypes"))



#Calculate the optimal cutoff
library(survival)
library(survminer)
res.cut<-surv_cutpoint(merged_data_sel,time="futime",event="fustat",variables="sample_subtypes")
summary(res.cut)
#                cutpoint statistic
#sample_subtypes        1  2.416011

pdf("Cluster_cutoff.pdf")
plot(res.cut,"sample_subtypes",palette="npg")
dev.off()

res.cat<-surv_categorize(res.cut)

fit=survfit(Surv(futime, fustat) ~sample_subtypes, data = res.cat)

head(res.cat)
diff<-survdiff(Surv(futime, fustat)~sample_subtypes, data = res.cat)
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
                   legend.labs=c("2", "1"),
                   legend.title="Cluster",
                   xlab="Time (years)",
                   break.time.by = 1,
                   risk.table.title="",
                   palette=c("red", "blue"),
                   risk.table=T,
                   risk.table.height=.25)
outFile="Cluster_survival.pdf"
pdf(file=outFile,onefile = FALSE,width = 8,height =5)
print(surPlot)
dev.off()


