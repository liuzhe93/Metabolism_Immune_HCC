setwd("F:/Ashley_Projects/Metabolism_Immune/03.ClusterTME/01_cluster/")
rm(list=ls())

#https://cloud.tencent.com/developer/article/2019378
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ConsensusClusterPlus")


load('../../01.IdentificationOfKeyMetabolicPathways/uni_matrix.Rdata')
logExp_process[1:5,1:5]
write.csv(logExp_process, "geneExp.csv", quote = F)
#row: gene
#column: sample
d<-logExp_process
genes_aa<-read.csv("../../02.PPInetworkHubgenes/01_Aminoacid/gene_set_01.csv")
genes_lp<-read.csv("../../02.PPInetworkHubgenes/04_Lipid/gene_set_04.csv")
genes_nt<-read.csv("../../02.PPInetworkHubgenes/05_Nucleotide/gene_set_05.csv")
geneslist_aa<-genes_aa$Metagene
geneslist_lp<-genes_lp$Metagene
geneslist_nt<-genes_nt$Metagene
genes_all<-c(geneslist_aa, geneslist_lp, geneslist_nt)
genes_all_uniq<-unique(genes_all)
length(genes_all_uniq)
#[1] 1194
write.csv(genes_all_uniq, "allgenes.csv", quote = F)
d_sel<-data.frame()
for(i in 1:12925){
  if(row.names(d)[i] %in% genes_all_uniq){
    temp<-d[i,]
    d_sel<-rbind(d_sel, temp)
  }
}
dim(d_sel)
#[1] 1054  390
genes_all_sel<-rownames(d_sel)
write.csv(genes_all_sel, "selgenes.csv", quote = F)

checkedgenes<-setdiff(genes_all_uniq, genes_all_sel)
length(checkedgenes)
#[1] 140
write.csv(checkedgenes, "checkedgenes.csv", quote = F)

data_tmp<-d_sel

mad(as.numeric(d_sel[1, ]))
mads<-apply(d_sel, 1, mad)
d_sel=d_sel[rev(order(mads))[1:1054],]
dim(d_sel)
d_sel = sweep(d_sel,1, apply(d_sel,1,median,na.rm=T))
dim(d_sel)
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
sample_subtypes <- ccres[[2]][["consensusClass"]]
table(sample_subtypes)
## sample_subtypes
##   1   2 
## 126 264
###126个样本是第1型，264个样本是第2型

###boxplot
# https://stackoverflow.com/questions/47651868/split-violin-plot-with-ggplot2-with-quantiles
library(ggplot2)
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin,
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             # Original function by Jan Gleixner (@jan-glx)
                             # Adjustments by Wouter van der Bijl (@Axeman)
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 1))
                               quantiles <- create_quantile_segment_frame(data, draw_quantiles, split = TRUE, grp = grp)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           }
)

create_quantile_segment_frame <- function(data, draw_quantiles, split = FALSE, grp = NULL) {
  dens <- cumsum(data$density) / sum(data$density)
  ecdf <- stats::approxfun(dens, data$y)
  ys <- ecdf(draw_quantiles)
  violin.xminvs <- (stats::approxfun(data$y, data$xminv))(ys)
  violin.xmaxvs <- (stats::approxfun(data$y, data$xmaxv))(ys)
  violin.xs <- (stats::approxfun(data$y, data$x))(ys)
  if (grp %% 2 == 0) {
    data.frame(
      x = ggplot2:::interleave(violin.xs, violin.xmaxvs),
      y = rep(ys, each = 2), group = rep(ys, each = 2)
    )
  } else {
    data.frame(
      x = ggplot2:::interleave(violin.xminvs, violin.xs),
      y = rep(ys, each = 2), group = rep(ys, each = 2)
    )
  }
}

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, position = position, 
        show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

data_tmp_t<-t(data_tmp)
data_tmp_t<-as.data.frame(data_tmp_t)
data_tmp_t$ID<-row.names(data_tmp_t)


genes_top10_aa<-read.csv("F:/Ashley_Projects/Metabolism_Immune/02.PPInetworkHubgenes/01_Aminoacid/top10_withdegree.csv")
genelist_top10_aa<-genes_top10_aa$Name
data_aa<-data_tmp_t[,genelist_top10_aa]
dim(data_aa)
#[1] 390  10

genes_top10_lp<-read.csv("F:/Ashley_Projects/Metabolism_Immune/02.PPInetworkHubgenes/04_Lipid/top10_withdegree.csv")
lp_genes<-c("MED14", "MED1", "CCNC", "MED31", "MED28", "MED10", "MED4", "MED15")
data_lp<-data_tmp_t[,lp_genes]
dim(data_lp)
#[1] 390   8

genes_top10_nt<-read.csv("F:/Ashley_Projects/Metabolism_Immune/02.PPInetworkHubgenes/05_Nucleotide/top10_withdegree.csv")
nt_genes<-c("NT5E", "HPRT1", "ITPA", "AMPD2", "AMPD3", "ATIC", "GMPR", "IMPDH2")
data_nt<-data_tmp_t[,nt_genes]
dim(data_nt)
#[1] 390   8

##########################metabolic pathways--aa#######################################
library("tidyr")
library("dplyr")
data_aa$ID<-row.names(data_aa)
tmp_aa <- data_aa %>% mutate(sample_subtypes = factor(sample_subtypes))
pdf("AA_diff_meta_path.pdf")
data_aa %>% 
  mutate(sample_subtypes = factor(sample_subtypes)) %>% 
  pivot_longer(-c(ID,sample_subtypes), names_to = "cell_type",values_to = "value") %>% 
  ggplot(aes(cell_type,value,fill=sample_subtypes))+
  geom_split_violin(draw_quantiles = c(0.25, 0.5, 0.75))+
  theme_bw()+
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45,hjust = 1)
  )
dev.off()

cluster1<-filter(tmp_aa, sample_subtypes == "1")
cluster2<-filter(tmp_aa, sample_subtypes == "2")

res_RPL14<-t.test(cluster1$RPL14,cluster2$RPL14)
res_RPL21<-t.test(cluster1$RPL21,cluster2$RPL21)
res_RPL24<-t.test(cluster1$RPL24,cluster2$RPL24)
res_RPL27<-t.test(cluster1$RPL27,cluster2$RPL27)
res_RPL31<-t.test(cluster1$RPL31,cluster2$RPL31)
res_RPL34<-t.test(cluster1$RPL34,cluster2$RPL34)
res_RPL35A<-t.test(cluster1$RPL35A,cluster2$RPL35A)
res_RPL36AL<-t.test(cluster1$RPL36AL,cluster2$RPL36AL)
res_RPL38<-t.test(cluster1$RPL38,cluster2$RPL38)
res_RPL6<-t.test(cluster1$RPL6,cluster2$RPL6)
res_RPL14$p.value
#[1] 6.805867e-17
res_RPL21$p.value
#[1] 1.839841e-13
res_RPL24$p.value
#[1] 1.295145e-28
res_RPL27$p.value
#[1] 2.07578e-33
res_RPL31$p.value
#[1] 3.153239e-24
res_RPL34$p.value
#[1] 2.38682e-14
res_RPL35A$p.value
#[1] 3.553479e-34
res_RPL36AL$p.value
#[1] 0.0001552919
res_RPL38$p.value
#[1] 1.405278e-23
res_RPL6$p.value
#[1] 4.55874e-33



##########################metabolic pathways--lp#######################################
data_lp$ID<-row.names(data_lp)
tmp_lp <- data_lp %>% mutate(sample_subtypes = factor(sample_subtypes))
pdf("LP_diff_meta_path.pdf")
data_lp %>% 
  mutate(sample_subtypes = factor(sample_subtypes)) %>% 
  pivot_longer(-c(ID,sample_subtypes), names_to = "cell_type",values_to = "value") %>% 
  ggplot(aes(cell_type,value,fill=sample_subtypes))+
  geom_split_violin(draw_quantiles = c(0.25, 0.5, 0.75))+
  theme_bw()+
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45,hjust = 1)
  )
dev.off()

cluster1<-filter(tmp_lp, sample_subtypes == "1")
cluster2<-filter(tmp_lp, sample_subtypes == "2")
res_MED14<-t.test(cluster1$MED14,cluster2$MED14)
res_MED1<-t.test(cluster1$MED1,cluster2$MED1)
res_CCNC<-t.test(cluster1$CCNC,cluster2$CCNC)
res_MED31<-t.test(cluster1$MED31,cluster2$MED31)
res_MED28<-t.test(cluster1$MED28,cluster2$MED28)
res_MED10<-t.test(cluster1$MED10,cluster2$MED10)
res_MED4<-t.test(cluster1$MED4,cluster2$MED4)
res_MED15<-t.test(cluster1$MED15,cluster2$MED15)
res_MED14$p.value
#[1] 0.8013701
res_MED1$p.value
#[1] 0.00749798
res_CCNC$p.value
#[1] 0.1471964
res_MED31$p.value
#[1] 0.4127845
res_MED28$p.value
#[1] 3.06814e-15
res_MED10$p.value
#[1] 3.078573e-18
res_MED4$p.value
#[1] 0.07612336
res_MED15$p.value
#[1] 4.407032e-11


##########################metabolic pathways--nt#######################################
data_nt$ID<-row.names(data_nt)
tmp_nt <- data_nt %>% mutate(sample_subtypes = factor(sample_subtypes))
pdf("NT_diff_meta_path.pdf")
data_nt %>% 
  mutate(sample_subtypes = factor(sample_subtypes)) %>% 
  pivot_longer(-c(ID,sample_subtypes), names_to = "cell_type",values_to = "value") %>% 
  ggplot(aes(cell_type,value,fill=sample_subtypes))+
  geom_split_violin(draw_quantiles = c(0.25, 0.5, 0.75))+
  theme_bw()+
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45,hjust = 1)
  )
dev.off()

cluster1<-filter(tmp_nt, sample_subtypes == "1")
cluster2<-filter(tmp_nt, sample_subtypes == "2")
res_NT5E<-t.test(cluster1$NT5E,cluster2$NT5E)
res_HPRT1<-t.test(cluster1$HPRT1,cluster2$HPRT1)
res_ITPA<-t.test(cluster1$ITPA,cluster2$ITPA)
res_AMPD2<-t.test(cluster1$AMPD2,cluster2$AMPD2)
res_AMPD3<-t.test(cluster1$AMPD3,cluster2$AMPD3)
res_AATIC<-t.test(cluster1$ATIC,cluster2$ATIC)
res_GMPR<-t.test(cluster1$GMPR,cluster2$GMPR)
res_IMPDH2<-t.test(cluster1$GMPR,cluster2$IMPDH2)
res_NT5E$p.value
#[1] 0.5030945
res_HPRT1$p.value
#[1] 0.471418
res_ITPA$p.value
#[1] 4.353597e-23
res_AMPD2$p.value
#[1] 2.336914e-05
res_AMPD3$p.value
#[1] 0.04505571
res_AATIC$p.value
#[1] 9.39401e-30
res_GMPR$p.value
#[1] 4.803444e-05
res_IMPDH2$p.value
#[1] 1.235971e-86

sample_cluster<-tmp_aa[, c("ID", "sample_subtypes")]
write.csv(sample_cluster, "Sample_Cluster.csv", quote = F, row.names = F)


