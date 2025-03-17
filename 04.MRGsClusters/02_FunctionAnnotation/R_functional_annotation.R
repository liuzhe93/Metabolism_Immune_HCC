setwd("F:/Ashley_Projects/Metabolism_Immune/04.MRGsClusters/02_FunctionAnnotation")
rm(list=ls())

library("Rgraphviz")
library("ggplot2")
library("stringr")
#install.packages("D:/ProgramFiles/Rpackages/AnnotationDbi_1.69.0.tar.gz", repos = NULL)
packageVersion("AnnotationDbi")
#install.packages("D:/ProgramFiles/Rpackages/org.Hs.eg.db_3.20.0.tar.gz", repos = NULL)
library("org.Hs.eg.db")
library("clusterProfiler")
library("topGO")
library("pathview")
library("openxlsx")

load("../01_DEGs/LIHCDEG.Rdata")
head(limma_voom_DEG)
cg3 = rownames(limma_voom_DEG)[limma_voom_DEG$change != "NOT"]
data_sig<-limma_voom_DEG[cg3,]
head(data_sig)
table(data_sig$change)
#DOWN   UP 
# 240  420 
write.csv(data_sig, "limma_res_sig.csv", quote = F)
table(limma_voom_DEG$change)
#DOWN   NOT    UP 
# 240 12327   420 
data_sig$gene_name<-row.names(data_sig)

keytypes(org.Hs.eg.db)
entrezid_all<-mapIds(x = org.Hs.eg.db,
                     keys = data_sig$gene_name,
                     keytype = "SYMBOL",
                     column = "ENTREZID")
entrezid_all<-na.omit(entrezid_all)
entrezid_all<-data.frame(entrezid_all)
head(entrezid_all)

#####################GO富集分析######################################
GO_enrich<-enrichGO(gene = entrezid_all[,1],
                    OrgDb = org.Hs.eg.db,
                    keyType = "ENTREZID",
                    ont = "ALL",
                    pvalueCutoff = 1, qvalueCutoff = 1,
                    readable = T)
GO_enrich<-data.frame(GO_enrich)
write.csv(GO_enrich, "GO_enrich.csv")

GO_enrich$term<-paste(GO_enrich$ID, GO_enrich$Description, sep = ": ")
GO_enrich$term<-factor(GO_enrich$term, levels = GO_enrich$term, ordered = T)

GO_enrich_BP<-GO_enrich[GO_enrich$ONTOLOGY=="BP",]
GO_enrich_CC<-GO_enrich[GO_enrich$ONTOLOGY=="CC",]
GO_enrich_MF<-GO_enrich[GO_enrich$ONTOLOGY=="MF",]
GO_enrich_top8<-rbind(GO_enrich_BP[1:8,], GO_enrich_CC[1:8,], GO_enrich_MF[1:8,])


#纵向柱状图
p3<-ggplot(GO_enrich_top8, aes(y = term, x = Count, fill = p.adjust)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(title = "GO Term Enrich",
       x = "Gene number",
       y = "GO Terms") +
  theme(axis.title.x = element_text(face = "bold", size = 4),
        axis.title.y = element_text(face = "bold", size = 4),
        legend.title = element_text(face = "bold", size = 16)) +
  theme_bw()
p3
pdf("GO_MFCCBP_v.pdf", height = 5, width = 20)
p3 + facet_grid(ONTOLOGY~., scale = "free_y", space = "free_y")
dev.off()

p4<-ggplot(GO_enrich_top8,
           aes(x = term, y = Count, fill = p.adjust)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_gradient(low = "blue", high = "red") +
  xlab("GO term") +
  ylab("Gene number") +
  labs(title = "GO Terms Enrich") +
  theme_bw() +
  theme(axis.text.x = element_text(family = "sans", face = "bold", color = "grey50", angle = 70, 
                                   vjust = 1, hjust = 1, size = 4))
p4
pdf("GO_MFCCBP_h.pdf")
p4 + facet_grid(.~ONTOLOGY, scales = "free_x", space = "free_x")
dev.off()

#########################metascape
#https://metascape.org/gp/index.html#/main/step1
#input genes #660 and submit
#Express Analysis
#Analysis Report Page


#####################KEGG富集分析######################################
#ref: https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/
ids <- bitr(data_sig$gene_name, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
dedup_ids <- ids[!duplicated(ids[c("SYMBOL")]),]
tmp=data.frame(data_sig$gene_name, logfc=data_sig$logFC)
tmp=merge(tmp, dedup_ids, by.x = "data_sig.gene_name", by.y = "SYMBOL")
colnames(tmp)[1]<-"SYMBOL"
CB_pathway=tmp$logfc
names(CB_pathway)=dedup_ids$ENTREZID
CB_pathway=sort(CB_pathway, decreasing = T)
kegg_gene_list<-CB_pathway
#https://www.genome.jp/kegg/tables/br08606.html
kegg_organism = "hsa"
kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = "hsa",
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")
pdf("KEGG_res.pdf")
dotplot(kk2, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
dev.off()
write.csv(kk2@result, "KEGG_res.csv", quote = F)
