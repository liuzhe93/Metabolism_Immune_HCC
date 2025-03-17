setwd("F:/Ashley_Projects/Metabolism_Immune/05.GC_TME/05_GESA_enrichment/")
rm(list=ls())

#参考资料：
#http://www.biocloudservice.com/wordpress/?p=54527

library(clusterProfiler)
library(ReactomePA)    #加载所需R包
library(tidyverse)
library(data.table)
library(org.Hs.eg.db)
library(clusterProfiler)
library(biomaRt)
library(enrichplot)

load("F:/Ashley_Projects/Metabolism_Immune/05.GC_TME/05_GESA_enrichment/LIHCDEG.Rdata")
dim(limma_voom_DEG)
#[1] 12987     7
head(limma_voom_DEG)
genelist_input<-limma_voom_DEG
genelist_input$genesymbol<-rownames(genelist_input)
head(genelist_input)
#            logFC   AveExpr        t      P.Value    adj.P.Val        B change genesymbol
#FLVCR1  0.9084112  3.570566 7.656272 1.509244e-13 1.090885e-09 20.24737     UP     FLVCR1
#UGT2B11 1.8810018  2.162512 7.580276 2.519947e-13 1.090885e-09 19.73964     UP    UGT2B11
#INTS7   0.4904496  4.833714 7.581572 2.498087e-13 1.090885e-09 19.73041    NOT      INTS7
#FSTL4   2.0159216 -1.235806 7.457451 5.726077e-13 1.859114e-09 18.09406     UP      FSTL4
#ACBD3   0.4154360  5.405486 7.228011 2.585033e-12 6.714363e-09 17.45649    NOT      ACBD3
#IQGAP3  1.3441181  3.068288 7.107915 5.612491e-12 9.308996e-09 16.75814     UP     IQGAP3
dim(genelist_input)
#[1] 12987     8
library("clusterProfiler")
gene_entrezid <- bitr(geneID = genelist_input$genesymbol
                      , fromType = "SYMBOL" # 从symbol
                      , toType = "ENTREZID" # 转成ENTREZID
                      , OrgDb = "org.Hs.eg.db"
)
head(gene_entrezid)
#   SYMBOL ENTREZID
#1  FLVCR1    28982
#2 UGT2B11    10720
#3   INTS7    25896
#4   FSTL4    23105
#5   ACBD3    64746
#6  IQGAP3   128239

#做GSEA分析对数据格式有要求，之前也说过，需要是一个有序的数值型向量，其名字是基因的ID
gene_entrezid <- merge(gene_entrezid,genelist_input,by.x = "SYMBOL", by.y = "genesymbol")
genelist <- gene_entrezid$logFC
names(genelist) <- gene_entrezid$ENTREZID
genelist <- sort(genelist,decreasing = T)
head(genelist)
#    4622     2719   124872    10814    57016    10642 
#3.442706 2.694464 2.613488 2.551697 2.458270 2.236153
library(msigdbr)
m_t2g <- msigdbr(species = "Homo sapiens", category = "C2") %>% 
  dplyr::select(gs_name, entrez_gene)
head(m_t2g)
## A tibble: 6 × 2
#gs_name               entrez_gene
#<chr>                       <int>
#1 HALLMARK_ADIPOGENESIS          19
#2 HALLMARK_ADIPOGENESIS       11194
#3 HALLMARK_ADIPOGENESIS       10449
#4 HALLMARK_ADIPOGENESIS          33
#5 HALLMARK_ADIPOGENESIS          34
#6 HALLMARK_ADIPOGENESIS          35
gsea_res <- GSEA(genelist, 
                 TERM2GENE = m_t2g,
                 minGSSize = 10,
                 maxGSSize = 500,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",
                 seed = 456
)
# 第一个条目的所有基因
gsea_res[[gsea_res$ID[[1]]]]
library(enrichplot)
library(ggplot2)
p <- gseaplot2(gsea_res, geneSetID = c(182, 215, 230, 385,  400), pvalue_table = TRUE)
write.csv(gsea_res, "GSEA_c2_result.csv")
pdf("GSEA_c2.pdf")
gseaplot2(gsea_res, geneSetID = 182, pvalue_table = TRUE)
dev.off()


