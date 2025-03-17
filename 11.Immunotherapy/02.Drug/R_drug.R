###R4.1.3
###Server 53

setwd("/mnt/zhangzheng_group/liuz-53/other/metabolism")
rm(list=ls())

#install R packages: https://blog.csdn.net/sk8_jarv/article/details/124384850
library(pRRophetic)

cluster<-read.csv("Sample_Cluster.csv")
LIHC_Expr<-read.csv("geneExp.csv", row.names = 1)
dim(LIHC_Expr)
#[1] 12925   390


expMatrix<-as.matrix(LIHC_Expr)
library(ggplot2)
data(cgp2016ExprRma) 
dim(cgp2016ExprRma)
#[1] 17419  1018

data(PANCANCER_IC_Tue_Aug_9_15_28_57_2016)
possibleDrugs2016 <- unique(drugData2016$Drug.name)
remove_drug<-c("PHA-665752","MG-132","VX-680","TAE684","Crizotinib","Saracatinib","S-Trityl-L-cysteine",
               "Z-LLNle-CHO","GNF-2","CGP-60474","CGP-082996","A-770041","WH-4-023","WZ-1-84","BI-2536",
               "BMS-536924","BMS-509744","JW-7-52-1","A-443654","MS-275","KIN001-135","TGX221","XMD8-85",
               "Mitomycin C","NSC-87877","CP466722","CHIR-99021","AP-24534","JNK-9L","PF-562271","HG-6-64-1",
               "JQ1","JQ12","FTI-277","OSU-03012","AKT inhibitor VIII","PAC-1","IPA-3","GSK-650394",
               "BAY 61-3606","5-Fluorouracil","Obatoclax Mesylate","BMS-754807","Lisitinib","LFM-A13",
               "GW-2580","Phenformin","Bryostatin 1","LAQ824","Epothilone B","GSK1904529A","BMS345541",
               "BMS-708163","Ruxolitinib","Ispinesib Mesylate","TL-2-105","AT-7519","TAK-715","BX-912",
               "ZSTK474","AS605240","Genentech Cpd 10","GSK1070916","KIN001-102","LY317615","GSK429286A","FMK",
               "QL-XII-47","CAL-101","UNC0638","XL-184","WZ3105","XMD14-99","AC220","CP724714","JW-7-24-1",
               "NPK76-II-72-1","STF-62247","NG-25","TL-1-85","VX-11e","FR-180204","Tubastatin A","Zibotentan",
               "YM155","NSC-207895","VNLG/124","AR-42","CUDC-101","Belinostat","I-BET-762","CAY10603",
               "Linifanib ","BIX02189","CH5424802","EKB-569","GSK2126458","KIN001-236","KIN001-244",
               "KIN001-055","KIN001-260","KIN001-266","Masitinib","MP470","MPS-1-IN-1","BHG712","OSI-930",
               "OSI-027","CX-5461","PHA-793887","PI-103","PIK-93","SB52334","TPCA-1","TG101348","Foretinib",
               "Y-39983","YM201636","Tivozanib","GSK690693","SNX-2112","QL-XI-92","XMD13-2","QL-X-138",
               "XMD15-27","T0901317","EX-527","THZ-2-49","KIN001-270","THZ-2-102-1","Navitoclax","CI-1040",
               "Olaparib","Veliparib","VX-702","AMG-706","KU-55933","Afatinib","GDC0449","BX-795","NU-7441",
               "SL 0101-1","BIRB 0796","JNK Inhibitor VIII","681640","Nutlin-3a (-)","PD-173074","ZM-447439",
               "RO-3306","MK-2206","PD-0332991","BEZ235","PD-0325901","selumetinib","EHT 1864","Cetuximab",
               "PF-4708671","JNJ-26854165","HG-5-113-01","HG-5-88-01","TW 37","XMD11-85h","ZG-10","XMD8-92",
               "QL-VIII-58","AG-014699","SB 505124","Tamoxifen","QL-XII-61","PFI-1","IOX2","YK 4-279",
               "(5Z)-7-Oxozeaenol","piperlongumine","FK866","Talazoparib","rTRAIL","UNC1215","SGC0946",
               "XAV939","Trametinib","Dabrafenib","Temozolomide","Bleomycin (50 uM)","SN-38","MLN4924",
               "AZD7762","GW 441756","CEP-701","SB 216763","17-AAG")
drug_list<-setdiff(possibleDrugs2016,remove_drug)
varnames<-paste("predictedPtype_",drug_list,sep = "")
for(drug_each in drug_list){
  print(drug_each)
  assign(paste0("predictedPtype_",drug_each),pRRopheticPredict(testMatrix=expMatrix, 
                                                               drug=drug_each,
                                                               tissueType = "all", 
                                                               batchCorrect = "eb",
                                                               selection=1,
                                                               dataset = "cgp2014"))
}
results<-get0(varnames[1])
for(i in 2:length(drug_list)){
  results<-cbind(results,get0(varnames[i]))
}
colnames(results)<-drug_list
write.csv(results,"drugSensitivity.csv",quote=F)
###R4.4
###Server 53

setwd("/mnt/zhangzheng_group/liuz-53/other/metabolism")
rm(list=ls())

#install R packages: https://blog.csdn.net/sk8_jarv/article/details/124384850
library(pRRophetic)



cluster<-read.csv("Sample_Cluster.csv")
LIHC_Expr<-read.csv("geneExp.csv", row.names = 1)
dim(LIHC_Expr)
#[1] 12925   390


expMatrix<-as.matrix(LIHC_Expr)
library(ggplot2)
data(cgp2016ExprRma) 
dim(cgp2016ExprRma)
#[1] 17419  1018

data(PANCANCER_IC_Tue_Aug_9_15_28_57_2016)
possibleDrugs2016 <- unique(drugData2016$Drug.name)
remove_drug<-c("PHA-665752","MG-132","VX-680","TAE684","Crizotinib","Saracatinib","S-Trityl-L-cysteine",
               "Z-LLNle-CHO","GNF-2","CGP-60474","CGP-082996","A-770041","WH-4-023","WZ-1-84","BI-2536",
               "BMS-536924","BMS-509744","JW-7-52-1","A-443654","MS-275","KIN001-135","TGX221","XMD8-85",
               "Mitomycin C","NSC-87877","CP466722","CHIR-99021","AP-24534","JNK-9L","PF-562271","HG-6-64-1",
               "JQ1","JQ12","FTI-277","OSU-03012","AKT inhibitor VIII","PAC-1","IPA-3","GSK-650394",
               "BAY 61-3606","5-Fluorouracil","Obatoclax Mesylate","BMS-754807","Lisitinib","LFM-A13",
               "GW-2580","Phenformin","Bryostatin 1","LAQ824","Epothilone B","GSK1904529A","BMS345541",
               "BMS-708163","Ruxolitinib","Ispinesib Mesylate","TL-2-105","AT-7519","TAK-715","BX-912",
               "ZSTK474","AS605240","Genentech Cpd 10","GSK1070916","KIN001-102","LY317615","GSK429286A","FMK",
               "QL-XII-47","CAL-101","UNC0638","XL-184","WZ3105","XMD14-99","AC220","CP724714","JW-7-24-1",
               "NPK76-II-72-1","STF-62247","NG-25","TL-1-85","VX-11e","FR-180204","Tubastatin A","Zibotentan",
               "YM155","NSC-207895","VNLG/124","AR-42","CUDC-101","Belinostat","I-BET-762","CAY10603",
               "Linifanib ","BIX02189","CH5424802","EKB-569","GSK2126458","KIN001-236","KIN001-244",
               "KIN001-055","KIN001-260","KIN001-266","Masitinib","MP470","MPS-1-IN-1","BHG712","OSI-930",
               "OSI-027","CX-5461","PHA-793887","PI-103","PIK-93","SB52334","TPCA-1","TG101348","Foretinib",
               "Y-39983","YM201636","Tivozanib","GSK690693","SNX-2112","QL-XI-92","XMD13-2","QL-X-138",
               "XMD15-27","T0901317","EX-527","THZ-2-49","KIN001-270","THZ-2-102-1","Navitoclax","CI-1040",
               "Olaparib","Veliparib","VX-702","AMG-706","KU-55933","Afatinib","GDC0449","BX-795","NU-7441",
               "SL 0101-1","BIRB 0796","JNK Inhibitor VIII","681640","Nutlin-3a (-)","PD-173074","ZM-447439",
               "RO-3306","MK-2206","PD-0332991","BEZ235","PD-0325901","selumetinib","EHT 1864","Cetuximab",
               "PF-4708671","JNJ-26854165","HG-5-113-01","HG-5-88-01","TW 37","XMD11-85h","ZG-10","XMD8-92",
               "QL-VIII-58","AG-014699","SB 505124","Tamoxifen","QL-XII-61","PFI-1","IOX2","YK 4-279",
               "(5Z)-7-Oxozeaenol","piperlongumine","FK866","Talazoparib","rTRAIL","UNC1215","SGC0946",
               "XAV939","Trametinib","Dabrafenib","Temozolomide","Bleomycin (50 uM)","SN-38","MLN4924",
               "AZD7762","GW 441756","CEP-701","SB 216763","17-AAG")
drug_list<-setdiff(possibleDrugs2016,remove_drug)
varnames<-paste("predictedPtype_",drug_list,sep = "")
for(drug_each in drug_list){
  print(drug_each)
  assign(paste0("predictedPtype_",drug_each),pRRopheticPredict(testMatrix=expMatrix, 
                                                               drug=drug_each,
                                                               tissueType = "all", 
                                                               batchCorrect = "eb",
                                                               selection=1,
                                                               dataset = "cgp2014"))
}
results<-get0(varnames[1])
for(i in 2:length(drug_list)){
  results<-cbind(results,get0(varnames[i]))
}
colnames(results)<-drug_list
write.csv(results,"drugSensitivity.csv",quote=F)

sessionInfo()
#R version 4.1.3 (2022-03-10)
#Platform: x86_64-conda-linux-gnu (64-bit)
#Running under: CentOS Linux 8

#Matrix products: default
#BLAS/LAPACK: /mnt/zhangzheng_group/liuz-53/miniconda3/envs/R4.1/lib/libopenblasp-r0.3.28.so

#locale:
# [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
# [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C            
#[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

#attached base packages:
#[1] stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#[1] ggplot2_3.4.2  pRRophetic_0.5

#loaded via a namespace (and not attached):
# [1] Rcpp_1.0.13-1          locfit_1.5-9.10        lattice_0.21-8        
# [4] png_0.1-8              Biostrings_2.62.0      utf8_1.2.3            
# [7] ridge_3.3              R6_2.5.1               GenomeInfoDb_1.30.1   
#[10] stats4_4.1.3           RSQLite_2.3.1          sva_3.42.0            
#[13] pillar_1.9.0           httr_1.4.6             zlibbioc_1.40.0       
#[16] rlang_1.1.1            annotate_1.72.0        car_3.1-3             
#[19] blob_1.2.4             S4Vectors_0.32.4       Matrix_1.5-4.1        
#[22] preprocessCore_1.56.0  splines_4.1.3          BiocParallel_1.28.3   
#[25] RCurl_1.98-1.12        bit_4.0.5              munsell_0.5.0         
#[28] compiler_4.1.3         pkgconfig_2.0.3        BiocGenerics_0.40.0   
#[31] mgcv_1.8-42            tidyselect_1.2.0       KEGGREST_1.34.0       
#[34] tibble_3.2.1           GenomeInfoDbData_1.2.7 edgeR_3.36.0          
#[37] IRanges_2.28.0         matrixStats_1.4.1      XML_3.99-0.14         
#[40] fansi_1.0.4            withr_2.5.0            dplyr_1.1.2           
#[43] crayon_1.5.2           bitops_1.0-7           grid_4.1.3            
#[46] nlme_3.1-162           xtable_1.8-4           gtable_0.3.3          
#[49] lifecycle_1.0.3        DBI_1.1.3              magrittr_2.0.3        
#[52] scales_1.2.1           cli_3.6.1              cachem_1.1.0          
#[55] carData_3.0-5          XVector_0.34.0         genefilter_1.76.0     
#[58] limma_3.50.3           generics_0.1.3         vctrs_0.6.5           
#[61] Formula_1.2-5          tools_4.1.3            bit64_4.0.5           
#[64] Biobase_2.54.0         glue_1.6.2             abind_1.4-8           
#[67] parallel_4.1.3         fastmap_1.2.0          survival_3.5-5        
#[70] AnnotationDbi_1.56.2   colorspace_2.1-0       memoise_2.0.1         

