setwd("F:/Ashley_Projects/Metabolism_Immune/04.MRGsClusters/01_DEGs")
rm(list=ls())

#limma,edgeR和DESeq2三大包基本是做转录组差异分析的金标准
#大多数转录组的文章都是用这三个R包进行差异分析
#edgeR差异分析速度快，得到的基因数目比较多，假阳性高（实际不差异结果差异）
#DESeq2差异分析速度慢，得到的基因数目比较少，假阴性高（实际差异结果不差异）

#需要注意的是制作分组信息的因子向量时，因子水平的前后顺序
#在R的很多模型中，默认将因子向量的第一个水平看作对照组

#logFC值:
#limma使用对数化的表达矩阵作为输入，所以将零假设指定在平均log值（对数几何平均数）。
#edgeR和DESeq2使用原始的count矩阵作为输入，所以将原假设指定在count的平均值上（算术平均数）。
#logFC一般是在1.2到2之间，筛选到差异基因的数目在500-1000左右为宜•
#根据logFC的统计指标确定阈值的方法是计算logFC绝对值的平均数与2倍标准差之和。（正态分布）

#DEG表示差异表达基因，differential gene express

library("stringr")
library("ggplotify")
library("patchwork")
library("cowplot")
library("DESeq2")
library("edgeR")
library("limma")

library(TCGAbiolinks)
library(DT)
library(dplyr)
library(SummarizedExperiment)
testData <- load("F:/Ashley_Projects/Metabolism_Immune/DataCollection/LIHC_case.rda")
testData
data
assays(data)
names(assays(data))[5]
#[1] "fpkm_unstrand"
names(assays(data))[1]
#[1] "unstranded"
dataPrep2 <- TCGAanalyze_Preprocessing(object = data,
                                       cor.cut = 0.6,
                                       datatype = "unstranded")
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
#[1] 12987   390
str(dataFilt)

write.csv(dataFilt,file = "TCGA_LIHC_final.csv",quote = FALSE)
dataFilt_LIHC_final <- read.csv("TCGA_LIHC_final.csv", header = T,check.names = FALSE)
# 定义行名
rownames(dataFilt_LIHC_final) <- dataFilt_LIHC_final[,1]
dataFilt_LIHC_final <- dataFilt_LIHC_final[,-1]
View(dataFilt_LIHC_final)
dim(dataFilt_LIHC_final)
#[1] 12987   390

sample_cluster<-read.csv("F:/Ashley_Projects/Metabolism_Immune/03.ClusterTME/01_cluster/Sample_Cluster.csv", header = T) 
sample_cluster$sample_subtypes<-paste0("C",sample_cluster$sample_subtypes)
head(sample_cluster)
dim(sample_cluster)
#[1] 390   2
#参考资料：
#https://www.cnblogs.com/Ixiaozhu/p/16736353.html
#https://developer.aliyun.com/article/929229
row.names(sample_cluster)<-sample_cluster[,1]
sample_cluster<-sample_cluster[c(-1)]
group_list<-factor(sample_cluster$sample_subtypes, levels = c("C1", "C2"))
design <- model.matrix( ~ 0 + group_list)
colnames(design)<-c("C1", "C2")
rownames(design)<-row.names(sample_cluster)

dge <- DGEList(counts=dataFilt_LIHC_final) %>%  # 将count表达矩阵转换为edgeR的DGEList对象
  calcNormFactors()  # 进行标准化和表达量计算。calcNormFactors函数并不会标准化数据，只是计算标准化因子。
# 标准化、线性模型拟合
fit <- voom(dge, design, normalize = "quantile") %>%  # voom是limma中的一种标准化方法
  lmFit(design)  # 线性模型拟合，出图
# 设置需要进行对比的分组
comp = paste(rev(levels(group_list)), collapse = "-")
cont.matrix <- makeContrasts(contrasts = comp, levels = design)
# 比较每个基因
fit2 <- contrasts.fit(fit, cont.matrix) %>%  # 根据对比模型进行差值计算
  eBayes()  # 贝叶斯检验
# 使用plotSA绘制了log2残差标准差与log-CPM均值的关系。平均log2残差标准差由水平蓝线标出
plotSA(fit2, main = "Final model: Mean-variance trend")
DEG = topTable(fit2, coef = comp, n=Inf) %>%  # 提取差异分析结果
  na.omit()  # 去除差异分析结果中包含NA值的行
length(which(DEG$adj.P.Val < 0.05))  # p值<0.05的基因有多少个
# 添加change列标记基因上调下调
logFC_cutoff <- with(DEG,mean(abs(logFC)) + 2*sd(abs(logFC)))  # 设置logFC的阈值
#limma包的一种方法，精确权重法（voom）,然后把筛选得到的差异表达基因写入csv文件中。
dge <- DGEList(counts = express_rec)
dge <- calcNormFactors(dge)#表达矩阵进行标准化；
v <- voom(dge, design,plot=TRUE)#利用limma_voom方法进行差异分析；
fit <- lmFit(v, design)#对数据进行线性拟合；
fit <- eBayes(fit)#贝叶斯算法组建
all <- topTable(fit, coef=ncol(design),n=Inf)#从高到低排名；
sig.limma <- subset(all_diff,abs(all$logFC)>1.5&all$P.Value<0.05)#进行差异基因筛选；
write.csv(sig.limma,'C:/Users/FREEDOM/Desktop/TCGA_data/limm_diff.csv')#写入csv文件中；
k1 = (DEG$P.Value < 0.05)&(DEG$logFC < -logFC_cutoff)  # 判断下调基因
k2 = (DEG$P.Value < 0.05)&(DEG$logFC > logFC_cutoff)  # 判断上调基因
DEG$change = ifelse(k1,"DOWN",ifelse(k2,"UP","NOT"))  # 标记上、下调基因
table(DEG$change)
head(DEG)
limma_voom_DEG <- DEG
boxplot(as.integer(dataFilt_LIHC_final[rownames(limma_voom_DEG)[1],]) ~ group_list)
## 取exp的第一行基因表达量，以group_list为分组依据，做箱线图。
limma_voom_DEG$logFC[1]
# 得到DESeq2_DEG的第一个logFC值。如果箱线图中normal大tumor小，则logFC应为负值，反之为正值。
save(limma_voom_DEG, group_list, dataFilt_LIHC_final, file = paste0("LIHC","DEG.Rdata"))


load("LIHCDEG.Rdata")
library(ggplot2)
library(tinyarray)
dataFilt_LIHC_final[1:4,1:4] # 原始count表达矩阵数据
dat = log2(dataFilt_LIHC_final+1)
pca.plot = draw_pca(dat, group_list); pca.plot # 作PCA图
pdf("LIHC_pca.pdf")
pca.plot
dev.off()

save(pca.plot, file = paste0("LIHC", "_PCA.Plot.Rdata")) # 保存PCA图数据以便后续拼图

cg3 = rownames(limma_voom_DEG)[limma_voom_DEG$change != "NOT"]
h3 = draw_heatmap(dat[cg3,], group_list, scale_before = T)
# 构建“均数+2倍标准差”函数
m2d = function(x){
  mean(abs(x))+2*sd(abs(x))
}
v3 = draw_volcano(limma_voom_DEG, pkg = 3, logFC_cutoff = m2d(limma_voom_DEG$logFC))
v3
pdf("LIHC_heatmap.pdf")
v3
dev.off()



