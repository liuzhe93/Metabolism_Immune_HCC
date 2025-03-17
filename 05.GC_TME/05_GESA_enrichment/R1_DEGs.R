setwd("F:/Ashley_Projects/Metabolism_Immune/05.GC_TME/05_GESA_enrichment")
rm(list=ls())

dataFilt_LIHC_final <- read.csv("F:/Ashley_Projects/Metabolism_Immune/04.MRGsClusters/01_DEGs/TCGA_LIHC_final.csv", header = T,check.names = FALSE)
# 定义行名
rownames(dataFilt_LIHC_final) <- dataFilt_LIHC_final[,1]
dataFilt_LIHC_final <- dataFilt_LIHC_final[,-1]
View(dataFilt_LIHC_final)
dim(dataFilt_LIHC_final)
#[1] 12987   390

sample_cluster<-read.csv("F:/Ashley_Projects/Metabolism_Immune/04.MRGsClusters/04_MRG_clusters/Sample_Cluster.csv", header = T) 
sample_cluster$sample_subtypes<-paste0("C",sample_cluster$cluster)
sample_cluster$cluster_update<-sample_cluster$sample_subtypes
for(i in 1:nrow(sample_cluster)){
  if(sample_cluster$sample_subtypes[i] == "C1"){
    sample_cluster$cluster_update[i]="C1C2"
  }else if(sample_cluster$sample_subtypes[i] == "C2"){
    sample_cluster$cluster_update[i]="C1C2"
  }else{
    sample_cluster$cluster_update[i]="C3"
  }
}

head(sample_cluster)
dim(sample_cluster)
#[1] 390   2
#参考资料：
#https://www.cnblogs.com/Ixiaozhu/p/16736353.html
#https://developer.aliyun.com/article/929229
row.names(sample_cluster)<-sample_cluster$sample_name
sample_cluster<-sample_cluster[c(-1,-2,-3)]
head(sample_cluster)
colnames(sample_cluster)<-"sample_subtypes"
head(sample_cluster)
group_list<-factor(sample_cluster$sample_subtypes, levels = c("C1C2", "C3"))
table(sample_cluster$sample_subtypes)
#C1C2   C3 
# 281  109
design <- model.matrix( ~ 0 + group_list)
head(design)
colnames(design)<-c("C1C2", "C3")
rownames(design)<-row.names(sample_cluster)
head(design)
#                             C1C2 C3
#TCGA-FV-A3I0-01A-11R-A22L-07    1  0
#TCGA-BD-A3ER-01A-11R-A213-07    1  0
#TCGA-CC-5261-01A-01R-A131-07    1  0
#TCGA-DD-AAVZ-01A-11R-A41C-07    0  1
#TCGA-DD-AADN-01A-11R-A41C-07    1  0
#TCGA-DD-A1EB-01A-11R-A131-07    1  0

library("stringr")
library("ggplotify")
library("patchwork")
library("cowplot")
library("DESeq2")
library("edgeR")
library("limma")

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
#[1] 2611
# 添加change列标记基因上调下调
logFC_cutoff <- with(DEG,mean(abs(logFC)) + 2*sd(abs(logFC)))  # 设置logFC的阈值
#limma包的一种方法，精确权重法（voom）,然后把筛选得到的差异表达基因写入csv文件中。
dge <- DGEList(counts = express_rec)
dge <- calcNormFactors(dge)#表达矩阵进行标准化；
v <- voom(dge, design,plot=TRUE)#利用limma_voom方法进行差异分析；
fit <- lmFit(v, design)#对数据进行线性拟合；
fit <- eBayes(fit)#贝叶斯算法组建
all <- topTable(fit, coef=ncol(design),n=Inf)#从高到低排名；
sig.limma <- subset(all,abs(all$logFC)>1.5&all$P.Value<0.05)#进行差异基因筛选；
write.csv(sig.limma,'limm_diff.csv')#写入csv文件中；
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
dim(limma_voom_DEG)
#[1] 12987     7

head(limma_voom_DEG)
#            logFC   AveExpr        t      P.Value    adj.P.Val        B change
#FLVCR1  0.9084112  3.570566 7.656272 1.509244e-13 1.090885e-09 20.24737     UP
#UGT2B11 1.8810018  2.162512 7.580276 2.519947e-13 1.090885e-09 19.73964     UP
#INTS7   0.4904496  4.833714 7.581572 2.498087e-13 1.090885e-09 19.73041    NOT
#FSTL4   2.0159216 -1.235806 7.457451 5.726077e-13 1.859114e-09 18.09406     UP
#ACBD3   0.4154360  5.405486 7.228011 2.585033e-12 6.714363e-09 17.45649    NOT
#IQGAP3  1.3441181  3.068288 7.107915 5.612491e-12 9.308996e-09 16.75814     UP




