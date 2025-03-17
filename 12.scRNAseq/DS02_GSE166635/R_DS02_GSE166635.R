setwd("F:/Ashley_Projects/Metabolism_Immune/12.scRNAseq/DS02_GSE166635/")
rm(list=ls())

library("ggplot2")
dat<-read.table("LIHC_GSE166635_CellMetainfo_table.tsv", sep = "\t", header = T)
dim(dat)
#[1] 22631    10
head(dat)
#############################malignancy##################################################
pdf("Mali_immu_stro.pdf")
ggplot(data = dat, mapping = aes(x = UMAP_1, y = UMAP_2, colour= Celltype..malignancy.)) + 
  geom_point(size = 0.5) + theme(legend.position = "top") + 
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"))

dev.off()

#############################cell type####################################################
pdf("Celltype.pdf")
ggplot(data = dat, mapping = aes(x = UMAP_1, y = UMAP_2, colour= Celltype..major.lineage.)) + 
  geom_point(size = 0.5) + theme(legend.position = "top") + 
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"))
dev.off()

############################Risk score##################################################
library("Seurat")
library("hdf5r")
library("tidyverse")
library("GSVA")
library("msigdbr")

dat<-Read10X_h5("LIHC_GSE166635_expression.h5")
dim(dat)
#[1] 21394 22631
dat <- CreateSeuratObject(dat, project = "LIHC_DS2", min.cells = 3)
dim(dat)
#[1] 21394 22631
dat <- NormalizeData(dat) %>% FindVariableFeatures() %>% ScaleData()

#gene_set<-c("PON1", "KIF20A", "HAVCR1", "GNG4", "FAM180A", "FAM163B")
gene_set<-c("PON1")
gene_set<-as.data.frame(gene_set)
#gene_set$y<-rep("Risk Factor",6)
gene_set$y<-rep("Risk Factor",1)
colnames(gene_set)<-c("Metagene","Item")
list<- split(as.matrix(gene_set)[,1], gene_set[,2])

## 测试counts数据
expr <- as.matrix(dat@assays$RNA$counts)

ssgsea_par <- ssgseaParam(expr, list)  # all other values are default values
ssgsea_scores <- gsva(ssgsea_par)

#system.time({res.counts = gsva(expr, list, method="ssgsea", parallel.sz=10)}) 
# 用户   系统   流逝 
#219.42   0.71 220.47

res.counts<-ssgsea_scores
riskfactor_score<-t(res.counts)
RFS<-data.frame(riskfactor_score)
RFS$Cell<-row.names(RFS)
celltype_data<-read.table("LIHC_GSE166635_CellMetainfo_table.tsv", sep = "\t", header = T)
dim(RFS)
#[1] 22631     2
dim(celltype_data)
#[1] 22631    10
input_data<-merge(celltype_data,RFS,by="Cell")
dim(input_data)
#[1] 22631    11
input_data$RiskFactorScore<-input_data$Risk.Factor
library("ggplot2")
pdf("RiskFactors_UMAP_score.pdf")
ggplot(data = input_data, mapping = aes(x = UMAP_1, y = UMAP_2, colour= RiskFactorScore)) + 
  geom_point(size = 0.5) + theme(legend.position = "top") +
  scale_color_gradient(low = "white", high = "red") +
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"))
dev.off()
write.csv(input_data,"RiskFactors_data_scRNA.csv",quote=F,row.names = F)

#####################################distribution##########################################
RFS_data<-read.csv("RiskFactors_data_scRNA.csv")

library("ggplot2")
library("ggpubr")
library("ggstatsplot")
library("ggsci")

pdf("RiskFactors_score.pdf",height=4,width=10)
ggplot(data = RFS_data, mapping = aes(x = Celltype..major.lineage., y = RiskFactorScore, fill= Celltype..major.lineage.)) +
  geom_violin() + geom_boxplot(width=0.2)+
  #  scale_fill_jco()+geom_jitter(shape=16,size=2,position=position_jitter(0.2)) +
  stat_compare_means()+guides(fill=FALSE)+theme_classic()+
  geom_point(size = 1) + theme(legend.position = "top") +
  theme(legend.position = "top", axis.text.x = element_text(angle = 30, hjust =1, vjust = 1)) +
  scale_color_gradient(low = "white", high = "red")
dev.off()



