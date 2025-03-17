setwd("F:/Ashley_Projects/Metabolism_Immune/06.MR_relatedRiskModel/01_RiskModel")
rm(list=ls())

library("glmnet")
library("survival")
rt<-read.table("F:/Ashley_Projects/Metabolism_Immune/04.MRGsClusters/03_LassoRegression/lassoSigExp.txt",
               header = T, row.names = 1)
head(rt)
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

pdf("riskScore_distribution.pdf",height = 10, width =6)
plot(res.cut,"riskScore",palette="npg")
dev.off()

res.cat<-surv_categorize(res.cut)
head(res.cat)
table(res.cat$riskScore)
#high  low 
#  73  282

sample_cluster<-cbind(rownames(res.cat), res.cat$fustat, res.cat$riskScore)
colnames(sample_cluster)<-c("SampleName", "fustat", "riskClass")
write.csv(sample_cluster, "Sample_Cluster.csv", row.names = F, quote = F)


fit<-survfit(Surv(futime,fustat)~riskScore,data=res.cat)
pdf("survival_curve.pdf",height = 6, width = 6)
ggsurvplot(fit,
           data=res.cat,
           risk.table = TRUE,
           pval=T)
dev.off()

outTab<-as.data.frame(outTab)
dim(outTab)
#[1] 7   6
outTab_matrix<-outTab
outTab_matrix<-outTab[, -1]
new1=apply(outTab_matrix,2,function(x) as.numeric(as.character(x)))
new1<-as.data.frame(new1)
row.names(new1)<-row.names(outTab_matrix)
outTab_filt<-subset(new1, pvalue <= 0.05)
dim(outTab_filt)
#[1] 6   5
write.csv(outTab_matrix,"outTab_matrix.csv",quote = F)
write.csv(outTab_filt,"outTab_sig.csv",quote = F)


############################visualization
library("ggplot2")
library("patchwork")
head(outTab_filt)
#                coef        HR    HR.95L    HR.95H       pvalue
#KIF20A   0.104559619 1.1102216 1.0322678 1.1940622 4.878348e-03
#PON1    -0.001765397 0.9982362 0.9965363 0.9999389 4.233813e-02
#GNG4     0.032699270 1.0332398 1.0090313 1.0580290 6.867043e-03
#HAVCR1   0.031005556 1.0314912 1.0130108 1.0503088 7.754526e-04
#FAM180A  0.159333903 1.1727295 1.0737712 1.2808077 3.964810e-04
#FAM163B  0.284075737 1.3285335 1.1768433 1.4997760 4.382561e-06

mydata<-outTab_filt
mydata$HR<-NULL
mydata$HR.95L<-NULL
mydata$HR.95H<-NULL
mydata$pvalue<-NULL
mydata$group<-ifelse(mydata$coef>0, "Pos", "Neg")
mydata$genes<-row.names(mydata)
mydata
#                coef group   genes
#KIF20A   0.104559619   Pos  KIF20A
#PON1    -0.001765397   Neg    PON1
#GNG4     0.032699270   Pos    GNG4
#HAVCR1   0.031005556   Pos  HAVCR1
#FAM180A  0.159333903   Pos FAM180A
#FAM163B  0.284075737   Pos FAM163B
ggplot(mydata, aes(genes, coef, fill = group)) + geom_bar(stat = "identity")

pdf("Coefficients.pdf")
ggplot(mydata, aes(genes, coef, fill = group)) + geom_bar(stat = "identity") + coord_flip() +
  scale_fill_manual(values = c("Pos" = "red", "Neg" = "green"), guide = FALSE) + xlab("Genes") +
  ylab("Coefficient") + theme_bw() + theme(panel.grid = element_blank()) + 
  theme(panel.border = element_rect(size = 0.6)) +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())
dev.off()


