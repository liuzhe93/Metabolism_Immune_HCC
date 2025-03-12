setwd("F:/Ashley_Projects/Metabolism_Immune/01.IdentificationOfKeyMetabolicPathways/08_EstimateScores")
rm(list=ls())

#install.packages("D:/ProgramFiles/Rpackages/estimate_1.0.13.tar.gz", repos = NULL, type = "source")
library("estimate")

x_rmNA_selected<-read.table("F:/Ashley_Projects/Metabolism_Immune/DataCollection/clinical_rmNA.txt", header = T)

dataExp <- read.csv("F:/Ashley_Projects/Metabolism_Immune/DataCollection/TCGA_LIHC_final.csv", header = T,check.names = FALSE,row.names = 1)
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
AllSamples<-colnames(dataExp)
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
dataExp_cancer <- subset(dataExp, select=c(intersect_cancer$x))
dataExp_normal <- subset(dataExp, select=c(intersect_normal$x))
dim(dataExp_cancer)
#[1] 12925   340
dim(dataExp_normal)
#[1] 12925    50
dataExp_process<-cbind(dataExp_cancer,dataExp_normal)
save(dataExp_process,barcode,file = 'fpkm_uni_matrix.Rdata')##把处理好的数据存好
write.table(dataExp_process, "dataexp.txt", sep = "\t", quote = F)

filterCommonGenes(input.f="dataexp.txt", 
                  output.f="LIHC_12925genes.gct", 
                  id="GeneSymbol")
#[1] "Merged dataset includes 8429 genes (1983 mismatched)."
estimateScore(input.ds = "LIHC_12925genes.gct",
              output.ds="LIHC_estimate_score.gct", 
              platform="illumina")
#[1] "1 gene set: StromalSignature  overlap= 122"
#[1] "2 gene set: ImmuneSignature  overlap= 132"
plotPurity(scores = "LIHC_12925genes.gct")
scores=read.table("LIHC_estimate_score.gct",skip = 2,header = T)
rownames(scores)=scores[,1]
scores=t(scores[,3:ncol(scores)])
dim(scores)
#[1] 390   3
scores<-as.data.frame(scores)
scores$sample_name<-row.names(scores)
scores$sample_name<-substr(scores$sample_name,1,12)
for(i in 1:390){
  scores$sample_name[i]<-gsub("\\.", "-", scores$sample_name[i])
}
head(scores)
#                             StromalScore ImmuneScore ESTIMATEScore  sample_name
#TCGA.FV.A3I0.01A.11R.A22L.07   -751.08441   -86.56871  -837.6531246 TCGA-FV-A3I0
#TCGA.BD.A3ER.01A.11R.A213.07   -124.01140   123.06654    -0.9448576 TCGA-BD-A3ER
#TCGA.CC.5261.01A.01R.A131.07    -87.94042   599.19969   511.2592656 TCGA-CC-5261
#TCGA.DD.AAVZ.01A.11R.A41C.07   -992.15123  -423.31182 -1415.4630509 TCGA-DD-AAVZ
#TCGA.DD.AADN.01A.11R.A41C.07   -928.81321  1369.00794   440.1947264 TCGA-DD-AADN
#TCGA.DD.A1EB.01A.11R.A131.07   -784.69262  -446.59174 -1231.2843646 TCGA-DD-A1EB
scores$TumourPurity<-cos(0.6049872018+0.0001467884*scores$ESTIMATEScore)
#the mathematical formula is from https://www.nature.com/articles/ncomms3612
dim(scores)
write.csv(scores, "scores_estimate.csv", quote = F)

scores_aa<-read.csv("F:/Ashley_Projects/Metabolism_Immune/01.IdentificationOfKeyMetabolicPathways/01_Aminoacid/Normalize.csv", row.names = 1)
scores_aa_df<-as.data.frame(t(scores_aa))
colnames(scores_aa_df)<-"amino_acid"
head(scores_aa_df)

scores_lp<-read.csv("F:/Ashley_Projects/Metabolism_Immune/01.IdentificationOfKeyMetabolicPathways/04_Lipid/Normalize.csv", row.names = 1)
scores_lp_df<-as.data.frame(t(scores_lp))
colnames(scores_lp_df)<-"lipid"
head(scores_lp_df)

scores_nt<-read.csv("F:/Ashley_Projects/Metabolism_Immune/01.IdentificationOfKeyMetabolicPathways/05_Nucleotide/Normalize.csv", row.names = 1)
scores_nt_df<-as.data.frame(t(scores_nt))
colnames(scores_nt_df)<-"nucleotide"
head(scores_nt_df)

dim(scores_aa_df)
#[1] 390   1
dim(scores_lp_df)
#[1] 390   1
dim(scores_nt_df)
#[1] 390   1

merged_scores_met<-cbind(scores_aa_df, scores_lp_df, scores_nt_df)
merged_scores_met$sample<-row.names(merged_scores_met)
head(merged_scores_met)

data_merged<-cbind(merged_scores_met, scores)
data_merged_sel<-subset(data_merged, select = c("amino_acid", "lipid", "nucleotide", "StromalScore", "ImmuneScore", "ESTIMATEScore", "TumourPurity"))
write.csv(data_merged_sel, "data_merged.csv", quote = F)


  
  
