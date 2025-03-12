setwd("F:/Ashley_Projects/Metabolism_Immune/02.PPInetworkHubgenes/")
rm(list=ls())

gene_set_01<-read.csv("F:/Ashley_Projects/Metabolism_Immune/01.IdentificationOfKeyMetabolicPathways/01_Aminoacid/TableS1_AminoAcid.csv")
colnames(gene_set_01)<-c("Metagene","Pathway")

gene_set_04<-read.csv("F:/Ashley_Projects/Metabolism_Immune/01.IdentificationOfKeyMetabolicPathways/04_Lipid/TableS1_Lipid.csv")
colnames(gene_set_04)<-c("Metagene","Pathway")

gene_set_05<-read.csv("F:/Ashley_Projects/Metabolism_Immune/01.IdentificationOfKeyMetabolicPathways/05_Nucleotide/TableS1_Nucleotide.csv")
colnames(gene_set_05)<-c("Metagene","Pathway")

dim(gene_set_01)
#[1] 348   2
dim(gene_set_04)
#[1] 766   2
dim(gene_set_05)
#[1] 90  2

merged_genes<-c(gene_set_01$Metagene, gene_set_04$Metagene, gene_set_05$Metagene)
length(merged_genes)
#[1] 1204
length(unique(merged_genes))
#[1] 1194

allgenes_uniq<-unique(merged_genes)
write.csv(allgenes_uniq, "allgenes.csv", quote = F, row.names = F)

write.csv(gene_set_01, "gene_set_01.csv", quote = F, row.names = F)
write.csv(gene_set_04, "gene_set_04.csv", quote = F, row.names = F)
write.csv(gene_set_05, "gene_set_05.csv", quote = F, row.names = F)

#https://string-db.org/cgi/input?sessionId=bE6w2HnPhFJA&input_page_active_form=multiple_identifiers
#organism； Homo sapiens
