setwd("F:/Ashley_Projects/Metabolism_Immune/04.MRGsClusters/04_MRG_clusters/")
rm(list=ls())

load('F:/Ashley_Projects/Metabolism_Immune/01.IdentificationOfKeyMetabolicPathways/uni_matrix.Rdata')
logExp_process[1:5,1:5]
write.csv(logExp_process, "geneExp.csv", quote = F)
#row: gene
#column: sample

d<-logExp_process
genes_sig_prognosis<-read.csv("F:/Ashley_Projects/Metabolism_Immune/04.MRGsClusters/03_LassoRegression/outTab_sig.csv")
genelist_sig_prognosis<-genes_sig_prognosis$X
length(genelist_sig_prognosis)
write.csv(genelist_sig_prognosis, "genelist_sig_prognosis.csv", quote = F)

d_sel<-data.frame()
for(i in 1:12925){
  if(row.names(d)[i] %in% genelist_sig_prognosis){
    temp<-d[i,]
    d_sel<-rbind(d_sel, temp)
  }
}
dim(d_sel)
#[1] 340 390

data_tmp<-d_sel
mad(as.numeric(d_sel[1,]))
mads<-apply(d_sel, 1, mad)
d_sel=d_sel[rev(order(mads))[1:340],]
dim(d_sel)
d_sel = sweep(d_sel,1, apply(d_sel,1,median,na.rm=T))
dim(d_sel)
#[1] 340 390
d_sel[1:5,1:5]
library(ConsensusClusterPlus)
boxplot(d[,1:20])
boxplot(d_sel[,1:20])
ccres <- ConsensusClusterPlus(as.matrix(d_sel),
                              maxK=9,
                              reps=1000,
                              pItem=0.8,
                              pFeature=1,
                              tmyPal = c("white","blue"),
                              title='ConsensusCluster/',
                              clusterAlg="km",
                              distance="euclidean",
                              seed=123456,
                              plot="pdf"
)
## end fraction
iclres <- calcICL(ccres,title="ics of ssgsea res")

###determine the optimal number of clusters
###PAC=proportion of ambiguous clustering模糊聚类比例
Kvec = 2:9
x1 = 0.1; x2 = 0.9 
PAC = rep(NA,length(Kvec)) 
names(PAC) = paste("K=",Kvec,sep="") 
for(i in Kvec){
  M = ccres[[i]]$consensusMatrix
  Fn = ecdf(M[lower.tri(M)])
  PAC[i-1] = Fn(x2) - Fn(x1)
}
optK = Kvec[which.min(PAC)]
optK
## [1] 2
#根据PAC和上面一致性聚类给出的图来看，分成2个亚型是最合适的

###the data after clustering
sample_subtypes <- ccres[[3]][["consensusClass"]]
table(sample_subtypes)
#sample_subtypes
# 1   2   3 
#88 193 109 
###88个样本是第1型，193个样本是第2型，109个样本是第3型

sample_subtypes<-as.data.frame(sample_subtypes)
sample_subtypes$sample_name<-row.names(sample_subtypes)
row.names(sample_subtypes)<-NULL
head(sample_subtypes)
colnames(sample_subtypes)<-c("cluster", "sample_name")
write.csv(sample_subtypes, "Sample_Cluster.csv", quote = F, row.names = F)

