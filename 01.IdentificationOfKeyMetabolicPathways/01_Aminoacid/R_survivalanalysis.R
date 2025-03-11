setwd("F:/Ashley_Projects/Metabolism_Immune/01.IdentificationOfKeyMetabolicPathways/")
rm(list=ls())

dataExp <- read.csv("F:/Ashley_Projects/Metabolism_Immune/DataCollection/TCGA_LIHC_final.csv", header = T,check.names = FALSE,row.names = 1)


##################################gene expression profile#############################################################
#读取表型信息
barcode_cancer<-read.csv("F:/Ashley_Projects/Metabolism_Immune/DataCollection/Cancer_barcode.csv")
barcode_normal<-read.csv("F:/Ashley_Projects/Metabolism_Immune/DataCollection/Normal_barcode.csv")
barcode<-rbind(barcode_cancer,barcode_normal)
dim(barcode)
#[1] 421   1
dim(barcode_cancer)
#[1] 371   1
dim(barcode_normal)
#[1] 50  1
barcode$y<-c(rep("Cancer",371),rep("Normal",50))
colnames(barcode)<-c("SampleID","SampleType")
head(barcode)
tail(barcode)
#处理两个文件
logExp <- log2(dataExp+1)#把表达矩阵log
#logExp <- dataExp
#logExp<- logExp[,order(colnames(logExp))]#根据行名排序表达矩阵
AllSamples<-colnames(logExp)
length(AllSamples)
#[1] 390
allSamples_t<-as.data.frame(AllSamples)
dim(allSamples_t)
colnames(allSamples_t)<-"x"
library(dplyr)
intersect_cancer<-dplyr::intersect(allSamples_t,barcode_cancer)
intersect_normal<-dplyr::intersect(allSamples_t,barcode_normal)
dim(intersect_cancer)
#[1] 340   1
dim(intersect_normal)
#[1] 50  1
logExp_cancer <- subset(logExp, select=c(intersect_cancer$x))
logExp_normal <- subset(logExp, select=c(intersect_normal$x))
dim(logExp_cancer)
#[1] 12925   340
dim(logExp_normal)
#[1] 12925    50
logExp_process<-cbind(logExp_cancer,logExp_normal)
save(logExp_process,barcode,file = 'uni_matrix.Rdata')##把处理好的数据存好


##################################clinical information profile#############################################################
load("F:/Ashley_Projects/Metabolism_Immune/DataCollection/LIHC_clinical.Rdata")
surv_data<-read.csv("F:/Ashley_Projects/Metabolism_Immune/DataCollection/TCGAbiolinks-LIHC-clinical.csv",header=T)

#alive的样本采用days_to_last_follow_up;对于dead样本，overall survival采用day_to_death
surv_data$os<-ifelse(surv_data$vital_status=='Alive',surv_data$days_to_last_follow_up,surv_data$days_to_death)
library(dplyr)
surv_selected<-surv_data %>%
  select(submitter_id,os,vital_status,gender,ajcc_pathologic_stage,ajcc_pathologic_t,ajcc_pathologic_m,ajcc_pathologic_n,age_at_index)
#0=alive, 1=dead.
surv_selected$fustat<-ifelse(surv_selected$vital_status=='Alive',0,1)
surv_selected$futime<-surv_selected$os
result_time<-surv_selected %>%
  select(submitter_id,futime,fustat,age_at_index,gender,ajcc_pathologic_stage,ajcc_pathologic_t,ajcc_pathologic_m,ajcc_pathologic_n)		
colnames(result_time)<-c("Id","futime","fustat","age","gender","stage","T","M","N")
result_time$gender_num<-ifelse(result_time$gender=='female',0,1)
result_time$stage_num<-ifelse(result_time$stage=='Stage I',1,
                              ifelse(result_time$stage=='Stage II',2,
                                     ifelse(result_time$stage=='Stage III',3,
                                            ifelse(result_time$stage=='Stage IIIA',3,
                                                   ifelse(result_time$stage=='Stage IIIB',3,
                                                          ifelse(result_time$stage=='Stage IIIC',3,
                                                                 ifelse(result_time$stage=='Stage IV',4,
                                                                        ifelse(result_time$stage=='Stage IVA',4,
                                                                               ifelse(result_time$stage=='Stage IVB',4,
                                                                                      "NA")))))))))
result_time$T_num<-ifelse(result_time$T=='T1',1,
                          ifelse(result_time$T=='T2',2,
                                 ifelse(result_time$T=='T2a',2,
                                        ifelse(result_time$T=='T2b',2,
                                               ifelse(result_time$T=='T3',3,
                                                      ifelse(result_time$T=='T3a',3,
                                                             ifelse(result_time$T=='T3b',3,
                                                                    ifelse(result_time$T=='T4',4,
                                                                           ifelse(result_time$T=='TX',5,
                                                                                  "NA")))))))))
result_time$M_num<-ifelse(result_time$M=='M0',0,
                          ifelse(result_time$M=='M1',1,
                                 ifelse(result_time$M=='MX',2,
                                        "NA")))
result_time$N_num<-ifelse(result_time$N=='N0',0,
                          ifelse(result_time$N=='N1',1,
                                 ifelse(result_time$N=='NX',2,
                                        "NA")))
result_time<-result_time %>%
  select(Id,futime,fustat,age,gender_num,stage_num,T_num,M_num,N_num)		
colnames(result_time)<-c("Id","futime","fustat","age","gender","stage","T","M","N")
write.table(result_time,"F:/Ashley_Projects/Metabolism_Immune/DataCollection/clinical.txt",row.names=F,quote=F)
x<-read.table("F:/Ashley_Projects/Metabolism_Immune/DataCollection/clinical.txt",header=T)
x_rmNA<-na.omit(x)
head(x_rmNA)
colnames(x_rmNA)[2]<-"futime_days"
x_rmNA$futime_years<-x_rmNA$futime_days/365
x_rmNA_selected<-x_rmNA %>%
  select(Id,futime_years,fustat,age,gender,stage,T,M,N)
colnames(x_rmNA_selected)[2]<-"futime"
write.table(x_rmNA_selected,"F:/Ashley_Projects/Metabolism_Immune/DataCollection/clinical_rmNA.txt",row.names=F,quote=F)


##载入数据
load('uni_matrix.Rdata')
gene_set<-read.csv("F:/Ashley_Projects/Metabolism_Immune/01.IdentificationOfKeyMetabolicPathways/01_Aminoacid/TableS1_AminoAcid.csv")
colnames(gene_set)<-c("Metagene","Pathway")
genelist<- split(as.matrix(gene_set)[,1], gene_set[,2])

library("GSVA")
ssgsea_par <- ssgseaParam(as.matrix(logExp_process), genelist)  # all other values are default values
ssgsea_scores <- gsva(ssgsea_par)
write.csv(ssgsea_scores,"F:/Ashley_Projects/Metabolism_Immune/01.IdentificationOfKeyMetabolicPathways/01_Aminoacid/noNormalize.csv",quote=F)

normalization<-function(x){
  return((x-min(x))/(max(x)-min(x)))}#设定normalization函数
nor_gsva_matrix1 <- normalization(ssgsea_scores)
write.csv(nor_gsva_matrix1,"F:/Ashley_Projects/Metabolism_Immune/01.IdentificationOfKeyMetabolicPathways/01_Aminoacid/Normalize.csv",quote=F)

dim(nor_gsva_matrix1)
#[1]   1 390
nor_gsva_matrix1_df<-as.data.frame(nor_gsva_matrix1)
pathway_data<-t(nor_gsva_matrix1_df)
pathway_data<-as.data.frame(pathway_data)
colnames(pathway_data)<-"Amino_acid"

for(i in 1:nrow(pathway_data)){
  temp<-rownames(pathway_data)[i]
  temp<-substr(temp,1,12)
  #  temp<-gsub("\\.","-",temp)
  pathway_data$Sample[i]<-temp
}
dim(pathway_data)
#[1] 390   2
summary(pathway_data)
#Amino acid        Sample         
#Min.   :0.0000   Length:390        
#1st Qu.:0.6622   Class :character  
#Median :0.7554   Mode  :character  
#Mean   :0.7234                     
#3rd Qu.:0.8244                     
#Max.   :1.0000 

merged_cli_score<-merge(x_rmNA_selected, pathway_data, by.x = "Id", by.y = "Sample")
dim(merged_cli_score)
#[1] 355  10

#Calculate the optimal cutoff
library(survival)
library(survminer)
res.cut<-surv_cutpoint(merged_cli_score,time="futime",event="fustat",variables="Amino_acid")
summary(res.cut)
#             cutpoint statistic
#Amino.acid  0.6621026  2.951503

pdf("F:/Ashley_Projects/Metabolism_Immune/01.IdentificationOfKeyMetabolicPathways/01_Aminoacid/Amino_acid_cutoff.pdf")
plot(res.cut,"Amino_acid",palette="npg")
dev.off()

res.cat<-surv_categorize(res.cut)
fit=survfit(Surv(futime, fustat) ~Carbohydrate, data = res.cat)
head(res.cat)
diff<-survdiff(Surv(futime, fustat)~Amino_acid, data = res.cat)
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
                   legend.title="Amino acid",
                   xlab="Time (years)",
                   break.time.by = 1,
                   risk.table.title="",
                   palette=c("red", "blue"),
                   risk.table=T,
                   risk.table.height=.25)
outFile="F:/Ashley_Projects/Metabolism_Immune/01.IdentificationOfKeyMetabolicPathways/01_Aminoacid/Amino_acid_survival.pdf"
pdf(file=outFile,onefile = FALSE,width = 8,height =5)
print(surPlot)
dev.off()


