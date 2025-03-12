setwd("F:/Ashley_Projects/Metabolism_Immune/02.PPInetworkHubgenes/00_external")
rm(list=ls())

#GEO: GSE124535
#clinical characteristics: https://www.nature.com/articles/s41586-019-0987-8#Sec47
#Supplementary Table 1

library("dplyr")
mydata<-read.table("GSE124535_HCC.RNA-seq.35.samples.fpkm.txt/GSE124535_HCC.RNA-seq.35.samples.fpkm.txt", 
                   header = T)
dim(mydata)
#[1] 63652    74
length(unique(mydata$gene_id))
#[1] 63652
length(unique(mydata$gene_symbol))
#[1] 56615
mydata_sel<-mydata[, -which(colnames(mydata) %in% c("gene_id", "Chr", "Biotype"))]
dim(mydata_sel)
#[1] 63652    71
mydata_sel_uniq<-aggregate(mydata_sel[,2:71], list(mydata_sel[,"gene_symbol"]), mean)
dim(mydata_sel_uniq)
#[1] 56615    71
row.names(mydata_sel_uniq)<-mydata_sel_uniq$Group.1
mydata_sel_uniq<-mydata_sel_uniq[,-1]

genelist_aa<-read.csv("../01_Aminoacid/top10_withdegree.csv")
genelist_lp<-read.csv("../04_Lipid/top10_withdegree.csv")
genelist_nt<-read.csv("../05_Nucleotide/top10_withdegree.csv")
gene_alas<-c("ADSSL1", "ADSS")

data_sel<-mydata_sel_uniq[c(genelist_aa$Name, genelist_lp$Name, genelist_nt$Name, gene_alas),]
dim(data_sel)
#[1] 32 70
data_sel_t<-t(data_sel)
data_sel_t<-as.data.frame(data_sel_t)
data_sel_t<-data_sel_t[, -which(colnames(data_sel_t) %in% c("NA", "NA.1"))]
data_sel_t$Sample<-substr(row.names(data_sel_t),1,4)
data_sel_t$Name<-row.names(data_sel_t)

cli_inf<-read.csv("2017-04-04932E-Supplementary Table 1.csv", header = T)
class(cli_inf)
table(cli_inf$RNA.Seqb)
# 0  1 
#75 35
table(cli_inf$Cancer)
# 0  1 
#70 40
cli_inf_sel<-subset(cli_inf, RNA.Seqb==1)
dim(cli_inf)
#[1] 111  23
dim(cli_inf_sel)
#[1] 35 23
cli_inf_sel<-cli_inf_sel[, c("Case.No.", "Disease.free.survival..m.", "Cancer.recurrenc.e.b")]

merged_data<-merge(data_sel_t, cli_inf_sel, by.x = "Sample", by.y = "Case.No.")
dim(merged_data)
#[1] 70 34
colnames(merged_data)[33]<-"Time"
colnames(merged_data)[34]<-"Status"
row.names(merged_data)<-merged_data$Name
merged_data<-merged_data[,-1]
write.csv(merged_data, "data4survival.csv", quote = F)

#col_sel<-cli_inf_sel$Case.No.
#cols<-c("gene_id", "gene_symbol", "Biotype")
#for(i in 1:35){
#  temp1<-paste0(col_sel[i],"P")
#  temp2<-paste0(col_sel[i],"T")
#  cols<-append(cols,temp1)
#  cols<-append(cols,temp2)
#}
#mydata_sel<-mydata[, cols]
#dim(mydata)
#[1] 63652    74
#dim(mydata_sel)
#[1] 63652    73
#mydata_sel_t<-t(mydata_sel)
#mydata_sel_t<-as.data.frame(mydata_sel_t)
#mydata_sel_t$Sample<-row.names(mydata_sel_t)


