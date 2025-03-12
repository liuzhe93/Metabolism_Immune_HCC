setwd("F:/Ashley_Projects/Metabolism_Immune/02.PPInetworkHubgenes/01_Aminoacid/")
rm(list=ls())

#STRING: Version: 12.0
#Cytoscape: Version 3.10.2

PPI_all<-read.table("string_interactions_short.tsv", header = T, sep = "\t")
head(PPI_all)
dim(PPI_all)
#[1] 9139   13
library("dplyr")

PPI_fil_highestcon<-subset(PPI_all, combined_score>0.90)
dim(PPI_fil_highestcon)
#[1] 5258   13
PPI_fil_highcon<-subset(PPI_all, combined_score>0.70)
dim(PPI_fil_highcon)
#[1] 6142   13
PPI_fil_medcon<-subset(PPI_all, combined_score>0.40)
dim(PPI_fil_medcon)
#[1] 9111   13
PPI_fil_lowcon<-subset(PPI_all, combined_score>0.15)
dim(PPI_fil_lowcon)
#[1] 9139   13

write.csv(PPI_fil_highcon, "AA_input_highcon.csv", quote = F, row.names = F)

#####Cytoscape output
top10<-read.csv("AA_input_highcon.csv_MCC_top10 default node.csv")

top10_genes<-top10$name
#[1] "RPL38"   "RPL35A"  "RPL36AL" "RPL24"   "RPL14"   "RPL6"    "RPL21"   "RPL34"   "RPL31"   "RPL27"
PPI_sel<-data.frame()
for (i in 1:dim(PPI_fil_highcon)[1]) {
  n1 <- PPI_fil_highcon[i, 1]
  n2 <- PPI_fil_highcon[i, 2]
  if(n1 %in% top10$name){
    if(n2 %in% top10$name){
      temp <- PPI_fil_highcon[i, c(1,2,13)]
      PPI_sel<-rbind(PPI_sel, temp)
    }
  }
}
write.csv(PPI_sel, "top10_withscore.csv", quote = F, row.names = F)


#####Tools-->Analyze Network-->treat the network as undirected
##name + Degree ==> Excel
##"name" column: 数据-->删除重复值
##top10_genes_degrees.csv
top10_degree<-read.csv("top10_genes_degrees.csv", header = T)

PPI_deg_sel<-data.frame()
for (i in 1:dim(top10_degree)[1]) {
  n <- top10_degree[i, 1]
  if(n %in% top10$name){
      temp <- top10_degree[i, c(1,2)]
      PPI_deg_sel<-rbind(PPI_deg_sel, temp)
  }
}
write.csv(PPI_deg_sel, "top10_withdegree.csv", quote = F, row.names = F)



###################################################################################
mydata<-read.csv("../00_external/data4survival.csv", row.names = 1)
#Calculate the optimal cutoff
library(survival)
library(survminer)
genelist_aa<-read.csv("top10_withdegree.csv")
genelist_aa$Name
mydata_aa<-mydata[,c(genelist_aa$Name)]



library(corrplot)
M = cor(mydata_aa, use = "everything", method = "spearman")
testRes = cor.mtest(mydata_aa,  use = "everything", method = "spearman", conf.level = 0.95)

pdf("Cor_AA.pdf")
corrplot(M, p.mat = testRes$p, method = 'circle', type = 'full', insig='blank',
         addCoef.col ='black', number.cex = 0.8, order = 'original', diag=T)
dev.off()


###################################################################################
mydata<-read.csv("../00_external/data4survival.csv", row.names = 1)
genelist_aa<-read.csv("top10_withdegree.csv")
genelist_aa$Name
mydata_aa<-mydata[,c(genelist_aa$Name, "Time", "Status")]
library("survival")
library("forestplot")
rt<-mydata_aa

outTab=data.frame()
for(i in colnames(rt[,1:10])){
  cox<-coxph(Surv(Time,Status)~rt[,i],data=rt)
  coxSummary=summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  outTab=rbind(outTab,
               cbind(id=i,
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  )
}
write.table(outTab,file="uniCox.xls",sep="\t",row.names = F,quote=F)

rt=read.table("uniCox.xls",header=T,sep="\t",row.names = 1,check.names = F)
data=as.matrix(rt)
HR=log2(data[,1:3])
hr=sprintf("%.3f",HR[,"HR"])
hrLow=sprintf("%.3f",HR[,"HR.95L"])
hrHigh=sprintf("%.3f",HR[,"HR.95H"])
pVal=data[,"pvalue"]
pVal=ifelse(pVal<0.001,"<0.001",sprintf("%.3f",pVal))
tabletext<-
  list(c(NA,rownames(HR)),
       append("pvalue",pVal),
       append("Hazard ratio",paste0(hr,"(",hrLow,"-",hrHigh,")")))
clrs<-fpColors(box="red",line="darkblue",summary="royalblue")
pdf(file="forest.pdf",onefile = FALSE,
    width = 8,
    height = 4)
forestplot(tabletext,
           rbind(rep(NA,3),HR),
           col=clrs,
           graphwidth=unit(50,"mm"),
           xlog=F,
           lwd.ci=2,
           boxsize=0.3,
           xlab="log2(Harard ratio(95%CI))")
dev.off()

data[,1:3]


