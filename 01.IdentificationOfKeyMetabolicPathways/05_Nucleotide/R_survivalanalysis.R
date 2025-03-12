setwd("F:/Ashley_Projects/Metabolism_Immune/01.IdentificationOfKeyMetabolicPathways/")
rm(list=ls())
x_rmNA_selected<-read.table("F:/Ashley_Projects/Metabolism_Immune/DataCollection/clinical_rmNA.txt", header = T)


##载入数据
load('uni_matrix.Rdata')
gene_set<-read.csv("F:/Ashley_Projects/Metabolism_Immune/01.IdentificationOfKeyMetabolicPathways/05_Nucleotide/TableS1_Nucleotide.csv")
colnames(gene_set)<-c("Metagene","Pathway")
genelist<- split(as.matrix(gene_set)[,1], gene_set[,2])

library("GSVA")
ssgsea_par <- ssgseaParam(as.matrix(logExp_process), genelist)  # all other values are default values
ssgsea_scores <- gsva(ssgsea_par)
write.csv(ssgsea_scores,"F:/Ashley_Projects/Metabolism_Immune/01.IdentificationOfKeyMetabolicPathways/05_Nucleotide/noNormalize.csv",quote=F)

normalization<-function(x){
  return((x-min(x))/(max(x)-min(x)))}#设定normalization函数
nor_gsva_matrix1 <- normalization(ssgsea_scores)
write.csv(nor_gsva_matrix1,"F:/Ashley_Projects/Metabolism_Immune/01.IdentificationOfKeyMetabolicPathways/05_Nucleotide/Normalize.csv",quote=F)

dim(nor_gsva_matrix1)
#[1]   1 390
nor_gsva_matrix1_df<-as.data.frame(nor_gsva_matrix1)
pathway_data<-t(nor_gsva_matrix1_df)
pathway_data<-as.data.frame(pathway_data)
colnames(pathway_data)<-"Nucleotide"

for(i in 1:nrow(pathway_data)){
  temp<-rownames(pathway_data)[i]
  temp<-substr(temp,1,12)
  #  temp<-gsub("\\.","-",temp)
  pathway_data$Sample[i]<-temp
}
dim(pathway_data)
#[1] 390   2
summary(pathway_data)
#  Nucleotide        Sample         
#Min.   :0.0000   Length:390        
#1st Qu.:0.4619   Class :character  
#Median :0.5854   Mode  :character  
#Mean   :0.5752                     
#3rd Qu.:0.6970                     
#Max.   :1.0000  

merged_cli_score<-merge(x_rmNA_selected, pathway_data, by.x = "Id", by.y = "Sample")
dim(merged_cli_score)
#[1] 355  10

#Calculate the optimal cutoff
library(survival)
library(survminer)
res.cut<-surv_cutpoint(merged_cli_score,time="futime",event="fustat",variables="Nucleotide")
summary(res.cut)
#            cutpoint statistic
#Nucleotide 0.4808013  3.603195

pdf("F:/Ashley_Projects/Metabolism_Immune/01.IdentificationOfKeyMetabolicPathways/05_Nucleotide/Nucleotide_cutoff.pdf")
plot(res.cut,"Nucleotide",palette="npg")
dev.off()

res.cat<-surv_categorize(res.cut)

fit=survfit(Surv(futime, fustat) ~Nucleotide, data = res.cat)

head(res.cat)
diff<-survdiff(Surv(futime, fustat)~Nucleotide, data = res.cat)
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
                   legend.title="Nucleotide",
                   xlab="Time (years)",
                   break.time.by = 1,
                   risk.table.title="",
                   palette=c("red", "blue"),
                   risk.table=T,
                   risk.table.height=.25)
outFile="F:/Ashley_Projects/Metabolism_Immune/01.IdentificationOfKeyMetabolicPathways/05_Nucleotide/Nucleotide_survival.pdf"
pdf(file=outFile,onefile = FALSE,width = 8,height =5)
print(surPlot)
dev.off()


