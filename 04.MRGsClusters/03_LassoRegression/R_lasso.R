setwd("F:/Ashley_Projects/Metabolism_Immune/04.MRGsClusters/03_LassoRegression/")
rm(list=ls())

library("glmnet")
library("survival")

####################read the gene expression profile matrix###################
dataFilt_LIHC_final <- read.csv("F:/Ashley_Projects/Metabolism_Immune/DataCollection/TCGA_LIHC_final.csv", row.names = 1)
dim(dataFilt_LIHC_final)
#[1] 12925   390
####################read the significant DEGs list############################
degs<-read.csv("../02_FunctionAnnotation/limma_res_sig.csv", row.names = 1)
head(degs)
dim(degs)
#[1] 660   7

genes_sel<-row.names(degs)
geneExp_sel<-dataFilt_LIHC_final[genes_sel,]
dim(geneExp_sel)
#[1] 660 390
geneExp_sel_rmNA<-na.omit(geneExp_sel)
dim(geneExp_sel_rmNA)
#[1] 629 390
geneExp_sel_t<-t(geneExp_sel_rmNA)
class(geneExp_sel_t)
#[1] "matrix" "array" 
geneExp_sel_t<-as.data.frame(geneExp_sel_t)
geneExp_sel_t$Id<-substring(rownames(geneExp_sel_t), 1, 12)
geneExp_sel_t$Id<-gsub("\\.", "-", geneExp_sel_t$Id)
######################read the clinical information###########################
sur_time<-read.table("F:/Ashley_Projects/Metabolism_Immune/DataCollection/clinical_rmNA.txt",header=T)
dim(sur_time)
#[1] 349   9
head(sur_time)

merged_data<-merge(sur_time, geneExp_sel_t, by = "Id")
dim(merged_data)
#[1] 355 638
for(i in 1:355){
  merged_data$Id[i]<-paste(merged_data$Id[i],i, sep = "_")
}
row.names(merged_data)<-merged_data$Id

rt<-merged_data
rt$futime<-rt$futime*365
rt$futime[rt$futime<=0]=1
rt$futime<-rt$futime/365

write.csv(rt,"uniSigExpID.csv",quote=F,row.names = F)


x<-as.matrix(rt[,c(10:ncol(rt))])
y<-data.matrix(Surv(rt$futime,rt$fustat))
fit<-glmnet(x,y,family="cox",maxit=1000)

pdf("lambda.pdf")
plot(fit,xvar="lambda",label=TRUE)
dev.off()

cvfit<-cv.glmnet(x,y,family="cox",maxit=1000)
pdf("cvfit.pdf")
plot(cvfit)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")
dev.off()

coef<-coef(fit,s=cvfit$lambda.min)
index<-which(coef!=0)
actCoef<-coef[index]
lassoGene<-row.names(coef)[index]
lassoGene<-c("futime","fustat",lassoGene)
lassoSigExp<-rt[,lassoGene]
lassoSigExp<-cbind(id=row.names(lassoSigExp),lassoSigExp)
write.table(lassoSigExp,file="lassoSigExp.txt",sep="\t",row.names = F,quote = F)

###################################multiCox
setwd("F:/Ashley_Projects/Metabolism_Immune/04.MRGsClusters/03_LassoRegression/")
rm(list=ls())
library("glmnet")
library("survival")
rt<-read.csv("uniSigExpID.csv")
rownames(rt)<-rt$Id
rt$Id<-NULL
rt$age<-NULL
rt$gender<-NULL
rt$stage<-NULL
rt$T<-NULL
rt$M<-NULL
rt$N<-NULL

multiCox<-coxph(Surv(futime,fustat) ~ ., data = rt)
multiCox<-step(multiCox,direction="both")
multiCoxSum<-summary(multiCox)
save(multiCox, file = "multiCox.RData")
save(multiCoxSum, file = "multiCoxSum.RData")
outTab<-data.frame()
outTab<-cbind(
  coef=multiCoxSum$coefficients[,"coef"],
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
outTab=gsub('"',"",outTab)
write.table(outTab,file="multiCox.xls",sep="\t",row.names = F,quote=F)

riskScore=predict(multiCox,type="risk",newdata=rt)
rt$riskScore<-riskScore
rt_saved<-rt
#rt_saved$riskCLass<-ifelse(rt$riskScore>1.144427,"high","low")
#write.table(rt_saved,"risk_class.txt",sep="\t",row.names = T, quote = F)
write.table(rt,"risk.txt",sep="\t",row.names = T, quote = F)
#计算最佳截点
rt=rt[,c(1,2,632)]

library("survminer")
res.cut<-surv_cutpoint(rt,
                       time = "futime",
                       event = "fustat",
                       variables = "riskScore")
summary(res.cut)
pdf("riskScore_distribution.pdf",height = 10, width =6)
plot(res.cut,"riskScore",palette="npg")
dev.off()

res.cat<-surv_categorize(res.cut)
head(res.cat)
table(res.cat$riskScore)
#high  low 
#  81  274
fit<-survfit(Surv(futime,fustat)~riskScore,data=res.cat)
pdf("survival_curve.pdf",height = 6, width = 6)
ggsurvplot(fit,
           data=res.cat,
           risk.table = TRUE,
           pval=T)
dev.off()

outTab<-as.data.frame(outTab)
dim(outTab)
#[1] 629   6
outTab_matrix<-outTab
outTab_matrix<-outTab[, -1]
new1=apply(outTab_matrix,2,function(x) as.numeric(as.character(x)))
new1<-as.data.frame(new1)
row.names(new1)<-row.names(outTab_matrix)
outTab_filt<-subset(new1, pvalue <= 0.05)
dim(outTab_filt)
#[1] 340   5
write.csv(outTab_matrix,"outTab_matrix.csv",quote = F)
write.csv(outTab_filt,"outTab_sig.csv",quote = F)

