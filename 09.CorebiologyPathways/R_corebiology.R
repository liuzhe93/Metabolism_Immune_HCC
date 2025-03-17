setwd("F:/Ashley_Projects/Metabolism_Immune/09.CorebiologyPathways/")
rm(list=ls())


genelist<-read.csv("F:/Ashley_Projects/Metabolism_Immune/01.IdentificationOfKeyMetabolicPathways/11_ssGSEA_corebiology/corebio_genes.csv")
head(genelist)
library("dplyr")

genelist_sel<-genelist %>% select(Gene.Signature, Genes) 
dim(genelist_sel)
#[1] 20  2

library("tidyverse")
dataExp<-read.table("F:/Ashley_Projects/Metabolism_Immune/06.MR_relatedRiskModel/05_EstimateScores/dataexp.txt", header = T, row.names = 1)

#01.CD8 T effector
genelist_01<-genelist_sel[1, ]
genelist_01
df<-genelist_01
split_strings <- strsplit(df$Genes, ", ")
df_new <- data.frame(matrix(unlist(split_strings), nrow = length(split_strings), byrow = TRUE))
path_name <- data.frame(matrix(rep("CD8 T effector", dim(df_new)[2]), byrow = TRUE))
df_new_t <- t(df_new)
merged_01 <- cbind(df_new_t, path_name)
colnames(merged_01) <- c("Metagene","Pathway")
merged_01

#03.Pan-F-TBRS
genelist_03<-genelist_sel[3, ]
genelist_03
df<-genelist_03
split_strings <- strsplit(df$Genes, ", ")
df_new <- data.frame(matrix(unlist(split_strings), nrow = length(split_strings), byrow = TRUE))
path_name <- data.frame(matrix(rep("Pan-F-TBRS", dim(df_new)[2]), byrow = TRUE))
df_new_t <- t(df_new)
merged_03 <- cbind(df_new_t, path_name)
colnames(merged_03) <- c("Metagene","Pathway")
merged_03

#04.Antigen processing machinery
genelist_04<-genelist_sel[4, ]
genelist_04
df<-genelist_04
split_strings <- strsplit(df$Genes, ", ")
df_new <- data.frame(matrix(unlist(split_strings), nrow = length(split_strings), byrow = TRUE))
path_name <- data.frame(matrix(rep("Antigen processing machinery", dim(df_new)[2]), byrow = TRUE))
df_new_t <- t(df_new)
merged_04 <- cbind(df_new_t, path_name)
colnames(merged_04) <- c("Metagene","Pathway")
merged_04

#05.Immune checkpoint
genelist_05<-genelist_sel[5, ]
genelist_05
df<-genelist_05
split_strings <- strsplit(df$Genes, ", ")
df_new <- data.frame(matrix(unlist(split_strings), nrow = length(split_strings), byrow = TRUE))
path_name <- data.frame(matrix(rep("Immune checkpoint", dim(df_new)[2]), byrow = TRUE))
df_new_t <- t(df_new)
merged_05 <- cbind(df_new_t, path_name)
colnames(merged_05) <- c("Metagene","Pathway")
merged_05

#06.EMT (1)
genelist_06<-genelist_sel[6, ]
genelist_06
df<-genelist_06
split_strings <- strsplit(df$Genes, ", ")
df_new <- data.frame(matrix(unlist(split_strings), nrow = length(split_strings), byrow = TRUE))
path_name <- data.frame(matrix(rep("EMT_1", dim(df_new)[2]), byrow = TRUE))
df_new_t <- t(df_new)
merged_06 <- cbind(df_new_t, path_name)
colnames(merged_06) <- c("Metagene","Pathway")
merged_06

#07.FGFR3-related genes
genelist_07<-genelist_sel[7, ]
genelist_07
df<-genelist_07
split_strings <- strsplit(df$Genes, ", ")
df_new <- data.frame(matrix(unlist(split_strings), nrow = length(split_strings), byrow = TRUE))
path_name <- data.frame(matrix(rep("FGFR3-related genes", dim(df_new)[2]), byrow = TRUE))
df_new_t <- t(df_new)
merged_07 <- cbind(df_new_t, path_name)
colnames(merged_07) <- c("Metagene","Pathway")
merged_07

#09.Angiogenesis
genelist_09<-genelist_sel[9, ]
genelist_09
df<-genelist_09
split_strings <- strsplit(df$Genes, ", ")
df_new <- data.frame(matrix(unlist(split_strings), nrow = length(split_strings), byrow = TRUE))
path_name <- data.frame(matrix(rep("Angiogenesis", dim(df_new)[2]), byrow = TRUE))
df_new_t <- t(df_new)
merged_09 <- cbind(df_new_t, path_name)
colnames(merged_09) <- c("Metagene","Pathway")
merged_09

#10.Fanconi anemia
genelist_10<-genelist_sel[10, ]
genelist_10
df<-genelist_10
split_strings <- strsplit(df$Genes, ", ")
df_new <- data.frame(matrix(unlist(split_strings), nrow = length(split_strings), byrow = TRUE))
path_name <- data.frame(matrix(rep("Fanconi anemia", dim(df_new)[2]), byrow = TRUE))
df_new_t <- t(df_new)
merged_10 <- cbind(df_new_t, path_name)
colnames(merged_10) <- c("Metagene","Pathway")
merged_10

#11.Cell cycle
genelist_11<-genelist_sel[11, ]
genelist_11
df<-genelist_11
split_strings <- strsplit(df$Genes, ", ")
df_new <- data.frame(matrix(unlist(split_strings), nrow = length(split_strings), byrow = TRUE))
path_name <- data.frame(matrix(rep("Cell cycle", dim(df_new)[2]), byrow = TRUE))
df_new_t <- t(df_new)
merged_11 <- cbind(df_new_t, path_name)
colnames(merged_11) <- c("Metagene","Pathway")
merged_11

#12.DNA replication
genelist_12<-genelist_sel[12, ]
genelist_12
df<-genelist_12
split_strings <- strsplit(df$Genes, ", ")
df_new <- data.frame(matrix(unlist(split_strings), nrow = length(split_strings), byrow = TRUE))
path_name <- data.frame(matrix(rep("DNA replication", dim(df_new)[2]), byrow = TRUE))
df_new_t <- t(df_new)
merged_12 <- cbind(df_new_t, path_name)
colnames(merged_12) <- c("Metagene","Pathway")
merged_12

#13.Nucleotide excision repair
genelist_13<-genelist_sel[13, ]
genelist_13
df<-genelist_13
split_strings <- strsplit(df$Genes, ", ")
df_new <- data.frame(matrix(unlist(split_strings), nrow = length(split_strings), byrow = TRUE))
path_name <- data.frame(matrix(rep("Nucleotide excision repair", dim(df_new)[2]), byrow = TRUE))
df_new_t <- t(df_new)
merged_13 <- cbind(df_new_t, path_name)
colnames(merged_13) <- c("Metagene","Pathway")
merged_13

#14.Homologous recombination
genelist_14<-genelist_sel[14, ]
genelist_14
df<-genelist_14
split_strings <- strsplit(df$Genes, ", ")
df_new <- data.frame(matrix(unlist(split_strings), nrow = length(split_strings), byrow = TRUE))
path_name <- data.frame(matrix(rep("Homologous recombination", dim(df_new)[2]), byrow = TRUE))
df_new_t <- t(df_new)
merged_14 <- cbind(df_new_t, path_name)
colnames(merged_14) <- c("Metagene","Pathway")
merged_14

#15.Mismatch repair
genelist_15<-genelist_sel[15, ]
genelist_15
df<-genelist_15
split_strings <- strsplit(df$Genes, ", ")
df_new <- data.frame(matrix(unlist(split_strings), nrow = length(split_strings), byrow = TRUE))
path_name <- data.frame(matrix(rep("Mismatch repair", dim(df_new)[2]), byrow = TRUE))
df_new_t <- t(df_new)
merged_15 <- cbind(df_new_t, path_name)
colnames(merged_15) <- c("Metagene","Pathway")
merged_15

#16.EMT2
genelist_16<-genelist_sel[16, ]
genelist_16
df<-genelist_16
split_strings <- strsplit(df$Genes, ", ")
df_new <- data.frame(matrix(unlist(split_strings), nrow = length(split_strings), byrow = TRUE))
path_name <- data.frame(matrix(rep("EMT_2", dim(df_new)[2]), byrow = TRUE))
df_new_t <- t(df_new)
merged_16 <- cbind(df_new_t, path_name)
colnames(merged_16) <- c("Metagene","Pathway")
merged_16

#17.EMT3
genelist_17<-genelist_sel[17, ]
genelist_17
df<-genelist_17
split_strings <- strsplit(df$Genes, ", ")
df_new <- data.frame(matrix(unlist(split_strings), nrow = length(split_strings), byrow = TRUE))
path_name <- data.frame(matrix(rep("EMT_3", dim(df_new)[2]), byrow = TRUE))
df_new_t <- t(df_new)
merged_17 <- cbind(df_new_t, path_name)
colnames(merged_17) <- c("Metagene","Pathway")
merged_17

#18.WNT target
genelist_18<-genelist_sel[18, ]
genelist_18
df<-genelist_18
split_strings <- strsplit(df$Genes, ", ")
df_new <- data.frame(matrix(unlist(split_strings), nrow = length(split_strings), byrow = TRUE))
path_name <- data.frame(matrix(rep("WNT target", dim(df_new)[2]), byrow = TRUE))
df_new_t <- t(df_new)
merged_18 <- cbind(df_new_t, path_name)
colnames(merged_18) <- c("Metagene","Pathway")
merged_18

#19.Cell cycle regulators
genelist_19<-genelist_sel[19, ]
genelist_19
df<-genelist_19
split_strings <- strsplit(df$Genes, ", ")
df_new <- data.frame(matrix(unlist(split_strings), nrow = length(split_strings), byrow = TRUE))
path_name <- data.frame(matrix(rep("Cell cycle regulators", dim(df_new)[2]), byrow = TRUE))
df_new_t <- t(df_new)
merged_19 <- cbind(df_new_t, path_name)
colnames(merged_19) <- c("Metagene","Pathway")
merged_19

#20.DNA damage repair
genelist_20<-genelist_sel[20, ]
genelist_20
df<-genelist_20
split_strings <- strsplit(df$Genes, ", ")
df_new <- data.frame(matrix(unlist(split_strings), nrow = length(split_strings), byrow = TRUE))
path_name <- data.frame(matrix(rep("DNA damage repair", dim(df_new)[2]), byrow = TRUE))
df_new_t <- t(df_new)
merged_20 <- cbind(df_new_t, path_name)
colnames(merged_20) <- c("Metagene","Pathway")
merged_20
#02.DNA damage repair
genelist_02<-genelist_sel[2, ]
genelist_02
df<-genelist_02
split_strings <- strsplit(df$Genes, ", ")
df_new <- data.frame(matrix(unlist(split_strings), nrow = length(split_strings), byrow = TRUE))
path_name <- data.frame(matrix(rep("DNA damage repair", dim(df_new)[2]), byrow = TRUE))
df_new_t <- t(df_new)
merged_02 <- cbind(df_new_t, path_name)
colnames(merged_02) <- c("Metagene","Pathway")
merged_02

merged_all <- rbind(merged_01, merged_03, merged_04, merged_05, merged_06, merged_07, merged_09, merged_10, merged_11, merged_12, merged_13, merged_14, merged_15, merged_16, 
                    merged_17, merged_18, merged_19, merged_20, merged_02)
table(merged_all$Pathway)
#Angiogenesis Antigen processing machinery               CD8 T effector                   Cell cycle        Cell cycle regulators            DNA damage repair 
#           4                            6                            8                          123                           10                          164 
#DNA replication                        EMT_1                        EMT_2                        EMT_3               Fanconi anemia          FGFR3-related genes 
#             36                            8                            7                            6                           51                            3 
#Homologous recombination            Immune checkpoint              Mismatch repair   Nucleotide excision repair                   Pan-F-TBRS                   WNT target 
#                      28                            7                           23                           44                           19                            4
list<-split(as.matrix(merged_all)[,1], merged_all[,2])
list

library("GSVA")
ssgsea_par <- ssgseaParam(as.matrix(dataExp), list)
gsva_matrix <- gsva(ssgsea_par)
dim(gsva_matrix)
#[1]  18 70
gsva_matrix<-as.data.frame(gsva_matrix)
gsva_matrix_t<-t(gsva_matrix)
dim(gsva_matrix_t)
#[1] 70  18
gsva_matrix_t<-as.data.frame(gsva_matrix_t)
gsva_matrix_t$sampleName<-row.names(gsva_matrix_t)

riskLabel<-read.table("F:/Ashley_Projects/Metabolism_Immune/06.MR_relatedRiskModel/05_EstimateScores/risk_class.txt")
riskLabel$sampleName<-c("L001P", "L001T", "L003T", "L003P", "L004T", 
                         "L004P", "L005T", "L005P", "L008P", "L008T", 
                         "L010T", "L010P", "L012T", "L012P", "L014T", 
                         "L014P", "L016P", "L016T", "L017T", "L017P", 
                         "L019T", "L019P", "L021T", "L021P", "L023T", 
                         "L023P", "L024P", "L024T", "L025T", "L025P", 
                         "L026T", "L026P", "L028T", "L028P", "L030P", 
                         "L030T", "L032T", "L032P", "L033T", "L033P", 
                         "L035T", "L035P", "L037P", "L037T", "L038T", 
                         "L038P", "L039T", "L039P", "L040T", "L040P", 
                         "L041T", "L041P", "L042P", "L042T", "L044T", 
                         "L044P", "L045T", "L045P", "L046P", "L046T", 
                         "L047T", "L047P", "L048T", "L048P", "L052T", 
                         "L052P", "L056P", "L056T", "L076T", "L076P")

merged_data<-merge(riskLabel, gsva_matrix_t, by = "sampleName")

merged_data_sel<-subset(merged_data, select = c("sampleName", "riskScore", "Angiogenesis", 
                                                "Antigen processing machinery", "CD8 T effector",
                                                "Cell cycle", "Cell cycle regulators", "DNA damage repair",
                                                "DNA replication", "EMT_1", "EMT_2", "EMT_3",
                                                "Fanconi anemia", "FGFR3-related genes", 
                                                "Homologous recombination", "Immune checkpoint", 
                                                "Mismatch repair", "Nucleotide excision repair", 
                                                "Pan-F-TBRS", "WNT target"))
dim(merged_data_sel)
#[1] 70 20
write.csv(merged_data_sel, "data_merged.csv", quote = F, row.names = F)

row.names(merged_data_sel)<-merged_data_sel$sampleName
merged_data_sel$sampleName<-NULL

M = cor(merged_data_sel, use = "everything", method = "spearman")
library(corrplot)
testRes = cor.mtest(merged_data_sel,  use = "everything", method = "spearman", conf.level = 0.95)

## add significant level stars
pdf("Cor_18corebio.pdf", height = 8, width = 8)
corrplot(M, p.mat = testRes$p, method = 'circle', type = 'lower', diag = T, insig='blank', pch.col = 'grey20',
         order = 'original')
dev.off()




