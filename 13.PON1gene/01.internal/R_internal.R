setwd("F:/Ashley_Projects/Metabolism_Immune/13.PON1gene/01.internal")
rm(list=ls())

dataExp<-read.csv("F:/Ashley_Projects/Metabolism_Immune/11.Immunotherapy/02.Drug/geneExp.csv", header = T, row.names = 1)
dim(dataExp)
#[1] 12925   390

dataExp_t<-t(dataExp)
dim(dataExp_t)
#[1]   390 12925
dataExp_t<-as.data.frame(dataExp_t)
dataExp_t$Id<-row.names(dataExp_t)
dataExp_t$Id<-substring(dataExp_t$Id, 1, 12)
dataExp_t$Id<-gsub("\\.", "-", dataExp_t$Id)

sur_time<-read.table("F:/Ashley_Projects/Metabolism_Immune/DataCollection/clinical_rmNA.txt",header=T)
dim(sur_time)
#[1] 349   9
head(sur_time)

merged_data<-merge(sur_time, dataExp_t, by = "Id")
merged_data_sel<-subset(merged_data, select = c("Id", "PON1", "futime", "fustat"))

#Calculate the optimal cutoff
library(survival)
library(survminer)
res.cut<-surv_cutpoint(merged_data_sel,time="futime",event="fustat",variables="PON1")
summary(res.cut)
#     cutpoint statistic
#PON1 6.108524  3.497298

pdf("Cluster_cutoff.pdf")
plot(res.cut,"PON1",palette="npg")
dev.off()

res.cat<-surv_categorize(res.cut)
test<-res.cat
identical(test$futime, merged_data_sel$futime)
#[1] TRUE
test$Id<-merged_data_sel$Id
write.csv(test, "res.csv", quote = F)

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


