setwd("F:/Ashley_Projects/Metabolism_Immune/06.MR_relatedRiskModel/05_EstimateScores/")
rm(list=ls())

geneExp<-read.table("F:/Ashley_Projects/Metabolism_Immune/06.MR_relatedRiskModel/02_ExternalDataSet_GSE124535/dataset/GSE124535_HCC.RNA-seq.35.samples.fpkm.txt", header = 1, row.names = 1)
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

row.names(geneExp_uniq)<-geneExp_uniq$Group.1
geneExp_uniq$Group.1<-NULL
geneExp_uniq[1:10,1:5]
#             L001P     L003P     L004P     L005P      L008P
#DELEC1    0.000000   0.00000  0.000000  0.000000  0.0423502
#MTARC1   22.003490  20.55836 25.453910 19.995550 21.4115950
#SEPT1     3.194365   0.72740  1.313468  0.737425  2.6545300
#MARCHF10  0.000000   0.00000  0.000000  0.134000  0.0000000
#SEPTIN10 11.819900  21.15380 15.012700 14.874800 20.0456000
#MARCHF11  0.000000   0.00000  0.000000  0.000000  0.0000000
#SEPTIN11  7.496620  21.12940 10.159500  9.047150 16.5686000
#SEPTIN12  0.102725   0.00000  0.000000  0.000000  0.0000000
#SEPTIN14  0.000000   0.00000  0.000000  0.000000  0.0000000
#SELENOF  57.303200 122.34700 74.037700 67.336600 80.4225000
dim(geneExp_uniq)
#[1] 56615    70
write.table(geneExp_uniq, "dataExp.txt", sep = "\t", quote = F)

library("estimate")
filterCommonGenes(input.f="dataexp.txt", 
                  output.f="LIHC_56615genes.gct", 
                  id="GeneSymbol")
#[1] "Merged dataset includes 10235 genes (177 mismatched)."
estimateScore(input.ds = "LIHC_56615genes.gct",
              output.ds="LIHC_estimate_score.gct", 
              platform="illumina")
#[1] "1 gene set: StromalSignature  overlap= 139"
#[1] "2 gene set: ImmuneSignature  overlap= 141"

plotPurity(scores = "LIHC_56615genes.gct")
scores=read.table("LIHC_estimate_score.gct",skip = 2,header = T)
rownames(scores)=scores[,1]
scores=t(scores[,3:ncol(scores)])
dim(scores)
#[1] 70   3
scores<-as.data.frame(scores)
scores$sample_name<-row.names(scores)
scores$sample_name<-substr(scores$sample_name,1,4)
head(scores)
#      StromalScore ImmuneScore ESTIMATEScore sample_name
#L001P   -457.52806    774.5618      317.0338        L001
#L003P   -196.62549   1265.7836     1069.1581        L003
#L004P   -215.87662   1101.6826      885.8060        L004
#L005P   -169.44134    901.3164      731.8750        L005
#L008P     21.71349   1597.9419     1619.6554        L008
#L010P    -11.58435   1295.7245     1284.1402        L010

scores$TumourPurity<-cos(0.6049872018+0.0001467884*scores$ESTIMATEScore)
#the mathematical formula is from https://www.nature.com/articles/ncomms3612
dim(scores)
write.csv(scores, "scores_estimate.csv", quote = F)



