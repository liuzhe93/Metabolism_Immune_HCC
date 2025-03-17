setwd("F:/Ashley_Projects/Metabolism_Immune/06.MR_relatedRiskModel/02_ExternalDataSet_GSE124535/")
rm(list=ls())

geneExp<-read.table("dataset/GSE124535_HCC.RNA-seq.35.samples.fpkm.txt", header = 1, row.names = 1)
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
#   Group.1     L001P    L003P     L004P
#1   DELEC1  0.000000  0.00000  0.000000
#2   MTARC1 22.003490 20.55836 25.453910
#3    SEPT1  3.194365  0.72740  1.313468
#4 MARCHF10  0.000000  0.00000  0.000000
#5 SEPTIN10 11.819900 21.15380 15.012700
#6 MARCHF11  0.000000  0.00000  0.000000

library("openxlsx")
clicalInfo<-read.xlsx("dataset/2017-04-04932E-Supplementary Table 1.xlsx", startRow = 2, colNames = TRUE)
dim(clicalInfo)
#[1] 111  22
clicalInfo[1:6,1:4]
#  Case.No. Gender Age HBV.a
#1     L001      M  67     P
#2     L003      M  45     P
#3     L004      M  42     P
#4     L005      M  51     P
#5     L006      M  57     P
#6     L007      M  56     N

clinicalInfo_sel<-clicalInfo[clicalInfo$`RNA-Seqb`=="1",]
dim(clinicalInfo_sel)
#[1] 36 22
dim(clicalInfo)
#[1] 111  22
head(clinicalInfo_sel)
clinicalInfo_sel<-subset(clinicalInfo_sel, select = c("Case.No.", "Gender", "Age", "Cancer.recurrenc.e.b",
                                                "Disease.free.survival.(m)"))
colnames(clinicalInfo_sel)<-c("Patient", "Gender", "Age", "CancerRecurrence", "DiseaseFreeSurvival_m")
head(clinicalInfo_sel)
#  Patient Gender Age CancerRecurrence DiseaseFreeSurvival_m
#1    L001      M  67                1                  19.5
#2    L003      M  45                0                  28.8
#3    L004      M  42                0                  28.7
#4    L005      M  51                0                  28.5
#7    L008      M  48                0                  27.4
#9    L010      M  25                0                  26.9
row.names(geneExp_uniq)<-geneExp_uniq$Group.1
geneExp_uniq$Group.1<-NULL
geneExp_t<-t(geneExp_uniq)
dim(geneExp_t)
#[1]    70 56615
geneExp_t[1:6,1:5]
#         DELEC1   MTARC1    SEPT1  MARCHF10 SEPTIN10
#L001P 0.0000000 22.00349 3.194365 0.0000000  11.8199
#L003P 0.0000000 20.55836 0.727400 0.0000000  21.1538
#L004P 0.0000000 25.45391 1.313468 0.0000000  15.0127
#L005P 0.0000000 19.99555 0.737425 0.1340000  14.8748
#L008P 0.0423502 21.41159 2.654530 0.0000000  20.0456
#L010P 0.0000000 20.57538 2.740135 0.0785101  19.8072
geneExp_t<-as.data.frame(geneExp_t)
geneExp_t$SampleID<-row.names(geneExp_t)
dim(geneExp_t)
#[1]    70 56616
geneExp_t$SampleID<-substr(geneExp_t$SampleID, 1, 4)
geneExp_t[1:6, 56610:56616]
#        ZYG11A ZYG11AP1  ZYG11B      ZYX   ZZEF1     ZZZ3 SampleID
#L001P 0.338463        0 3.53786 19.31170 4.04251  5.01862     L001
#L003P 0.161476        0 9.88813  8.37228 2.25808 11.12160     L003
#L004P 0.369763        0 7.93874 11.51890 2.69169  7.04392     L004
#L005P 0.470220        0 7.72440 15.40630 3.70013  9.13343     L005
#L008P 0.349502        0 9.11658 20.82670 5.45587  8.65731     L008
#L010P 0.273014        0 7.38196 28.80320 5.00410  7.17857     L010

merged_data<-merge(x = geneExp_t, y = clinicalInfo_sel, by.x = "SampleID", by.y = "Patient")
dim(merged_data)
#[1]    70 56620
merged_data[1:6, c(1:8,56614:56620)]
#  SampleID   DELEC1   MTARC1    SEPT1  MARCHF10 SEPTIN10 MARCHF11 SEPTIN11      ZYX   ZZEF1     ZZZ3 Gender
#1     L001 0.000000 22.00349 3.194365 0.0000000  11.8199        0  7.49662 19.31170 4.04251  5.01862      M
#2     L001 0.110243  5.89019 3.078259 0.0169441  17.1612        0 24.30710 48.29770 4.36698  9.99004      M
#3     L003 0.000000 23.53306 0.968030 0.0000000  21.9539        0 36.29560 29.28320 2.27734 17.21810      M
#4     L003 0.000000 20.55836 0.727400 0.0000000  21.1538        0 21.12940  8.37228 2.25808 11.12160      M
#5     L004 0.000000 10.77281 1.474450 0.0000000  14.1617        0 22.43910 50.57010 2.76070 15.40040      M
#6     L004 0.000000 25.45391 1.313468 0.0000000  15.0127        0 10.15950 11.51890 2.69169  7.04392      M
# Age CancerRecurrence DiseaseFreeSurvival_m
#1  67                1                  19.5
#2  67                1                  19.5
#3  45                0                  28.8
#4  45                0                  28.8
#5  42                0                  28.7
#6  42                0                  28.7
merged_data_sel<-merged_data[, c(1:8,56614:56620)]
write.csv(merged_data_sel, "clinical_infor.csv", quote = F, row.names = F)
for(i in 1:nrow(merged_data_sel)){
  merged_data_sel$id[i]<-paste(merged_data_sel$SampleID[i], i, sep = "_")
}
head(merged_data_sel)
rownames(merged_data_sel)<-merged_data_sel$id
merged_data_sel$id<-NULL
merged_data_sel$SampleID<-NULL
rt<-cbind(merged_data_sel$DiseaseFreeSurvival_m, merged_data_sel$CancerRecurrence, merged_data_sel[, 1:10])
colnames(rt)[1:2]<-c("futime", "fustat")
#       futime fustat   DELEC1   MTARC1    SEPT1  MARCHF10 SEPTIN10 MARCHF11 SEPTIN11      ZYX   ZZEF1
#L001_1   19.5      1 0.000000 22.00349 3.194365 0.0000000  11.8199        0  7.49662 19.31170 4.04251
#L001_2   19.5      1 0.110243  5.89019 3.078259 0.0169441  17.1612        0 24.30710 48.29770 4.36698
#L003_3   28.8      0 0.000000 23.53306 0.968030 0.0000000  21.9539        0 36.29560 29.28320 2.27734
#L003_4   28.8      0 0.000000 20.55836 0.727400 0.0000000  21.1538        0 21.12940  8.37228 2.25808
#L004_5   28.7      0 0.000000 10.77281 1.474450 0.0000000  14.1617        0 22.43910 50.57010 2.76070
#L004_6   28.7      0 0.000000 25.45391 1.313468 0.0000000  15.0127        0 10.15950 11.51890 2.69169
#           ZZZ3
#L001_1  5.01862
#L001_2  9.99004
#L003_3 17.21810
#L003_4 11.12160
#L004_5 15.40040
#L004_6  7.04392
library("glmnet")
library("survival")
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
rt=rt[,c(1,2,13)]
head(rt)
library("survminer")
res.cut<-surv_cutpoint(rt,
                       time = "futime",
                       event = "fustat",
                       variables = "riskScore")
summary(res.cut)
#          cutpoint statistic
#riskScore 6.636094   3.49593
pdf("riskScore_distribution.pdf",height = 10, width =6)
plot(res.cut,"riskScore",palette="npg")
dev.off()

res.cat<-surv_categorize(res.cut)
head(res.cat)
table(res.cat$riskScore)
#high  low 
#   7   63
fit<-survfit(Surv(futime,fustat)~riskScore,data=res.cat)
pdf("survival_curve.pdf",height = 6, width = 6)
ggsurvplot(fit,
           data=res.cat,
           risk.table = TRUE,
           pval=T)
dev.off()

