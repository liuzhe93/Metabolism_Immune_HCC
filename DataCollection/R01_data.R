setwd("F:/Ashley_Projects/Metabolism_Immune/DataCollection")
rm(list=ls())
library(TCGAbiolinks)
library(DT)
library(dplyr)
library(SummarizedExperiment)
TCGAbiolinks::getGDCprojects()$project_id
getGDCprojects()$project_id
cancer_type="TCGA-LIHC"
#####################################临床信息##################################################
clinical<-GDCquery_clinic(project=cancer_type,type="clinical")
dim(clinical)
#[1] 377  69
#View(clinical)
save(clinical,file="LIHC_clinical.Rdata")
write.csv(clinical, file="TCGAbiolinks-LIHC-clinical.csv",row.names = F)
data_type <- "Gene Expression Quantification"
data_category <- "Transcriptome Profiling"
query <- GDCquery(project = cancer_type,
                  data.category = data_category,
                  data.type = data_type,
                  workflow.type = "STAR - Counts")
samplesDown <- getResults(query,cols=c("cases"))  
#getResults(query, rows, cols)根据指定行名或列名从query中获取结果,此处用来获得样本的barcode
# 从samplesDown中筛选出TP（实体肿瘤）样本的barcodes
# TCGAquery_SampleTypes(barcode, typesample)
# TP代表PRIMARY SOLID TUMOR；NT-代表Solid Tissue Normal（其他组织样本可参考学习文档）
dataSmTP <- TCGAquery_SampleTypes(barcode = samplesDown,
                                  typesample = "TP")
# 从samplesDown中筛选出NT(正常组织)样本的barcode
dataSmNT <- TCGAquery_SampleTypes(barcode = samplesDown,
                                  typesample = "NT")
length(samplesDown)
#[1] 424
length(dataSmNT)
#[1] 50
length(dataSmTP)
#[1] 371
###设置barcodes参数，筛选符合要求的371个肿瘤样本数据和50个正常组织数据
write.csv(dataSmTP,"Cancer_barcode.csv",quote=F,row.names=F)
write.csv(dataSmNT,"Normal_barcode.csv",quote=F,row.names=F)

######################表达谱数据下载########################################################
queryDown <- GDCquery(project = cancer_type,
                 data.category = "Transcriptome Profiling",
                 data.type = "Gene Expression Quantification",
                 workflow.type = "STAR - Counts",
                 barcode = c(dataSmTP, dataSmNT))
GDCdownload(queryDown, method = "client", directory = "GDCdata", files.per.chunk = 10)
dataPrep1 <- GDCprepare(query = queryDown, save = TRUE, save.filename = "LIHC_case.rda")
testData <- load("LIHC_case.rda")
testData
data
assays(data)
names(assays(data))[5]
dataPrep2 <- TCGAanalyze_Preprocessing(object = data,
                                      cor.cut = 0.6,
                                      datatype = "fpkm_unstrand")
write.csv(dataPrep2, file = "LIHC_dataPrep.csv", quote = F)

purityDATA <- TCGAtumor_purity(colnames(data), 0, 0, 0, 0, 0.6)
Purity.LIHC<-purityDATA$pure_barcodes
normal.LIHC<-purityDATA$filtered
puried_data <-dataPrep2[,c(Purity.LIHC,normal.LIHC)]
library("SummarizedExperiment")
rowData(data)
rownames(puried_data)<-rowData(data)$gene_name
write.csv(puried_data,file = "puried.LIHC.csv",quote = FALSE)

library("EDASeq")
dataNorm <- TCGAbiolinks::TCGAanalyze_Normalization(tabDF = puried_data,
                                                    geneInfo = geneInfo,
                                                    method = "gcContent")
#将标准化后的数据再过滤，去除掉表达量较低（count较低）的基因，得到最终的数据
dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  method = "quantile", 
                                  qnt.cut =  0.25)
dim(dataNorm)
#[1] 17317   390
dim(dataFilt)
#[1] 12925   390
str(dataFilt)
write.csv(dataFilt,file = "TCGA_LIHC_final.csv",quote = FALSE)
intersect(dataSmTP,colnames(dataFilt_LIHC_final))
intersect(dataSmNT,colnames(dataFilt_LIHC_final))
#保留的是390个样本（前340个样本来源于肿瘤，后50个样本来源于正常组织）

#读入表达矩阵文件
dataFilt_LIHC_final <- read.csv("TCGA_LIHC_final.csv", header = T,check.names = FALSE)
# 定义行名
rownames(dataFilt_LIHC_final) <- dataFilt_LIHC_final[,1]
dataFilt_LIHC_final <- dataFilt_LIHC_final[,-1]
View(dataFilt_LIHC_final)



