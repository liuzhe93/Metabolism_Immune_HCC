setwd("F:/Ashley_Projects/Metabolism_Immune/03.ClusterTME/06_heatmap_28immunecells")
rm(list=ls())


mydata_imm_sco<-read.csv("F:/Ashley_Projects/Metabolism_Immune/01.IdentificationOfKeyMetabolicPathways/09_ssGSEA_immunecells/data_merged.csv")
head(mydata_imm_sco)
dim(mydata_imm_sco)
#[1] 390  32


mydata_cluster<-read.csv("F:/Ashley_Projects/Metabolism_Immune/03.ClusterTME/01_cluster/Sample_Cluster.csv")
head(mydata_cluster)
dim(mydata_cluster)
#[1] 390   2


sur_time<-read.table("F:/Ashley_Projects/Metabolism_Immune/DataCollection/clinical_rmNA.txt",header=T)
head(sur_time)
dim(sur_time)
#[1] 349   9

merged_clu_imm<-cbind(mydata_cluster, mydata_imm_sco)
merged_clu_imm$Id<-substring(merged_clu_imm$ID, 1, 12)
head(merged_clu_imm)

merged_data<-merge(sur_time, merged_clu_imm, by = "Id")
dim(merged_data)
#[1] 355  43
merged_data_sel<-subset(merged_data, select = c("ID", "sample_subtypes", "futime", "fustat", 
                                                "age", "gender", "stage", "T", "M", "N",
                                                "Activated.B.cell", "Activated.CD4.T.cell", "Activated.CD8.T.cell",
                                                "Activated.dendritic.cell", "CD56bright.natural.killer.cell", 
                                                "CD56dim.natural.killer.cell", "Central.memory.CD4.T.cell","Central.memory.CD8.T.cell",
                                                "Effector.memeory.CD4.T.cell", "Effector.memeory.CD8.T.cell", "Eosinophil", 
                                                "Gamma.delta.T.cell", "Immature..B.cell", "Immature.dendritic.cell", "Macrophage",
                                                "Mast.cell", "MDSC", "Memory.B.cell", "Monocyte", "Natural.killer.cell",
                                                "Natural.killer.T.cell", "Neutrophil","Plasmacytoid.dendritic.cell", "Regulatory.T.cell",            
                                                "T.follicular.helper.cell", "Type.1.T.helper.cell", "Type.17.T.helper.cell", 
                                                "Type.2.T.helper.cell"))

merged_srt<-merged_data_sel[order(merged_data_sel[,"sample_subtypes"]),]
rownames(merged_srt)<-merged_srt$ID
merged_srt<-merged_srt[,-1]
merged_srt_t<-t(merged_srt)
dim(merged_srt_t)
#[1]  37 355


genelist<-read.csv("F:/Ashley_Projects/Metabolism_Immune/01.IdentificationOfKeyMetabolicPathways/09_ssGSEA_immunecells/ImmunecellMetagenes.csv")
table(genelist$Immunity)
#Adaptive   Innate 
#     431      351
loc_adap<-which(genelist[,3]=="Adaptive")
loc_adap
genelist_adap<-genelist[loc_adap,]
head(genelist_adap)
loc_inna<-which(genelist[,3]=="Innate")
loc_inna
genelist_inna<-genelist[loc_inna,]
head(genelist_inna)
celltypes_adap<-unique(genelist_adap$Cell.type)
#[1] "Activated B cell"            "Activated CD4 T cell"        "Activated CD8 T cell"       
#[4] "Central memory CD4 T cell"   "Central memory CD8 T cell"   "Effector memeory CD4 T cell"
#[7] "Effector memeory CD8 T cell" "Gamma delta T cell"          "Immature  B cell"           
#[10] "Memory B cell"               "Regulatory T cell"           "T follicular helper cell"   
#[13] "Type 1 T helper cell"        "Type 17 T helper cell"       "Type 2 T helper cell"       
celltypes_inna<-unique(genelist_inna$Cell.type)
#[1] "Activated dendritic cell"       "CD56bright natural killer cell" "CD56dim natural killer cell"   
#[4] "Eosinophil"                     "Immature dendritic cell"        "Macrophage"                    
#[7] "Mast cell"                      "MDSC"                           "Monocyte"                      
#[10] "Natural killer cell"            "Natural killer T cell"          "Neutrophil"                    
#[13] "Plasmacytoid dendritic cell" 


row.names(merged_srt_t)<-gsub("\\.", " ", rownames(merged_srt_t))
data_adap<-merged_srt_t[celltypes_adap, ]
data_inna<-merged_srt_t[celltypes_inna, ]
dim(data_adap)
#[1]  15 355
dim(data_inna)
#[1]  13 355
data_allimmune<-rbind(data_adap, data_inna)
dim(data_allimmune)
#[1]  28 35
data_label<-merged_srt_t[c("sample_subtypes", "futime", "fustat", "age", "gender", "stage", "T", "M", "N"),]
data_label[1:9,1:3]
data_sel<-rbind(data_label, data_allimmune)

library(pheatmap)

annotation_col_info<-as.matrix(data_sel[1:9,])
annotation_col<-t(annotation_col_info)
annotation_col<-as.data.frame(annotation_col)
annotation_col$sample_subtypes<-paste0("C", annotation_col$sample_subtypes)
annotation_col$fustat<-ifelse(annotation_col$fustat==0, "alive", "dead")
annotation_col$gender<-ifelse(annotation_col$gender==0, "female", "male")

annotation_col$stage<-ifelse(annotation_col$stage==1, "I",
                             ifelse(annotation_col$stage==2, "II",
                                    ifelse(annotation_col$stage==3, "III",
                                           ifelse(annotation_col$stage==4, "IV",
                                                  "NA"))))
annotation_col$T<-ifelse(annotation_col$T==1, "T1",
                         ifelse(annotation_col$T==2, "T2",
                                ifelse(annotation_col$T==3, "T3",
                                       ifelse(annotation_col$T==4, "T4",
                                              "NA"))))
annotation_col$M<-ifelse(annotation_col$M==0, "M0",
                         ifelse(annotation_col$M==1, "M1",
                                ifelse(annotation_col$M==2, "MX",
                                       "NA")))
annotation_col$N<-ifelse(annotation_col$N==0, "N0",
                         ifelse(annotation_col$N==1, "N1",
                                ifelse(annotation_col$N==2, "NX",
                                       "NA")))
colnames(annotation_col)[1]<-"Cluster"
head(annotation_col)

annotation_row_adap<-cbind(as.data.frame(celltypes_adap), as.data.frame(rep("Adaptive immunity",15)))
colnames(annotation_row_adap)<-c("ImmuneCells", "Immunity")
annotation_row_inna<-cbind(as.data.frame(celltypes_inna), as.data.frame(rep("Innate immunity",13)))
colnames(annotation_row_inna)<-c("ImmuneCells", "Immunity")
annotation_row<-rbind(annotation_row_adap, annotation_row_inna)
row.names(annotation_row)<-annotation_row$ImmuneCells
annotation_row$ImmuneCells<-NULL
head(annotation_row)

input_data<-data_sel[10:37,1:355]
write.csv(data_sel, "Cluster_immunecells.csv", quote = F)


pdf("heatmap_28immune.pdf")
pheatmap(input_data, cluster_rows = F, cluster_cols = F, show_colnames = F, show_rownames = T, scale = "row",
         annotation_col = annotation_col, annotation_row = annotation_row, border_color = NA, gaps_row = 15)
dev.off()
