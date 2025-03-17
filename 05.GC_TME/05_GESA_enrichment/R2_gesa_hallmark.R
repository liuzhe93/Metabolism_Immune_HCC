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
m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
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
p <- gseaplot2(gsea_res, geneSetID = 4:5, pvalue_table = TRUE)
write.csv(gsea_res, "GSEA_result.csv")
pdf("GSEA.pdf")
p
dev.off()




# 
# 
# 
# genename <- as.character(rownames(genelist_input))   #提取第一列基因名
# gene_map <- select(org.Hs.eg.db, keys=genename, keytype="SYMBOL", columns=c("ENTREZID"))  #将SYMBOL格式的ID换成ENTREZ格式的ID。
# dim(gene_map)
# #[1] 12988     2
# length(genename)
# #[1] 12987
# 
# non_duplicates_idx <- which(duplicated(gene_map$SYMBOL) == FALSE)  
# gene_map <- gene_map[non_duplicates_idx, ]   #去除重复值
# colnames(gene_map)[1]<-"Gene"  #输出mapping结果
# dim(gene_map)
# #[1] 12987     2
# 
# genelist_input$Gene<-row.names(genelist_input)
# 
# temp<-inner_join(gene_map,genelist_input,by = "Gene")
# temp<-temp[,-1]
# temp<-na.omit(temp)
# temp$logFC<-sort(temp$logFC,decreasing = T)
# 
# geneList = temp[,2]
# names(geneList) = as.character(temp[,1])
# head(geneList)
# 
# library("msigdbr")
# 
# m_t2g<-msigdbr(species = "Homo sapiens", category = "C5") %>%
#   dplyr::select(gs_name, entrez_gene)
# head(m_t2g)
# 
# gsea_res <- GSEA(geneList, 
#                  TERM2GENE = m_t2g,
#                  minGSSize = 10,
#                  maxGSSize = 500,
#                  pvalueCutoff = 0.05,
#                  pAdjustMethod = "BH",
#                  seed = 456
# )
# # 第一个条目的所有基因
# gsea_res[[gsea_res$ID[[1]]]]
# 
# library(enrichplot)
# library(ggplot2)
# ridgeplot(gsea_res,
#           showCategory = 20,
#           fill = "p.adjust", #填充色 "pvalue", "p.adjust", "qvalue" 
#           core_enrichment = TRUE,#是否只使用 core_enriched gene
#           label_format = 30,
#           orderBy = "NES",
#           decreasing = FALSE
# )+
#   theme(axis.text.y = element_text(size=8))
# ## Picking joint bandwidth of 0.212
# 
# #展示基因集的logFC分布
# ids <- gsea_res@result$ID[10:15]
# gseadist(gsea_res,
#          IDs = ids,
#          type="density" # boxplot
# )+
#   theme(legend.direction = "vertical")
# 
# #展示基因的排序以及富集分数的变化。
# gsearank(gsea_res,
#          geneSetID = 1 # 要展示的基因集
# )
# 
# #这个函数还可以直接返回画图数据:
# aa <- gsearank(gsea_res, 1, title = gsea_res[1, "Description"],output = "table")
# head(aa)
# 
# #gseaplot函数可以画两个图：ES或者ranked-gene-list，通过参数by设置，默认是两个图都画出来
# #如果by="runningScore"，则是画出ES的图，如果是by = "preranked",则是画出ranked gene list的图，
# p <- gseaplot(gsea_res, geneSetID = 1, by = "preranked", 
#               title = gsea_res$Description[1])
# p
# p <- gseaplot(gsea_res, geneSetID = 1, by = "runningScore", 
#               title = gsea_res$Description[1])
# p
# 
# 
# Go_gseresult <- gseGO(geneList, 'org.Hs.eg.db', keyType = "ENTREZID", ont="BP", pvalueCutoff=0.05)   #使用GSEA进行GO富集分析（'org.Hs.eg.db'：对应物种的数据库；ont：选择输出条目，可选“BP,MF,CC或者ALL”，pvalueCutoff：设置P的阈值）
# KEGG_gseresult <- gseKEGG(geneList, pvalueCutoff=0.05)  #使用GSEA进行KEGG富集分析    
# go_results<-as.data.frame(Go_gseresult)
# kegg_results<-as.data.frame(KEGG_gseresult)
# write.csv (go_results, file ="GO_gseresult.csv")
# write.csv (kegg_results, file ="KEGG_gseresult.csv")
# 
# 
# gseaplot2(Go_gseresult, 1:3, title = "Specific GO Biological Process in G123 group", pvalue_table = T)  #1:3：这表示在图上显示前3条富集结果，也可以根据自己分析需要指定输出某一条结果；Go_gseresult：GO富集分析结果；title：加上标题；pvalue_table：是否在图上显示P值列表。
# 
# 
# fgseaRes1 <- fgsea(examplePathways[[175]], exampleRanks)
# plotEnrichment(kegg_results)
# 
# 
# gseaplot2(kegg_results, 1:3, title = "Specific KEGG signaling pathways in G123 group", pvalue_table = T)
# 
# #######################################################
# #reference: https://cran.r-project.org/web/packages/msigdbr/vignettes/msigdbr-intro.html
# library("msigdb")
# library("msigdbr")
# library("fgsea")
# all_gene_sets = msigdbr(species = "Homo sapiens")
# head(all_gene_sets)
# msigdbr_species()
# 
# h_gene_sets = msigdbr(species = "human", category = "H")
# head(h_gene_sets)
# msigdbr_collections()
# all_gene_sets %>%
#   dplyr::filter(gs_cat == "H") %>%
#   head()
# msigdbr_df <- msigdbr(species = "human", category = "H")
# head(msigdbr_df)
# msigdbr_t2g = msigdbr_df %>% dplyr::distinct(gs_name, entrez_gene) %>% as.data.frame()
# enricher(gene = geneList, TERM2GENE = msigdbr_t2g)
# 
# 
# 


