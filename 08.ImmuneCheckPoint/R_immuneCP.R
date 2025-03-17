setwd("F:/Ashley_Projects/Metabolism_Immune/08.ImmuneCheckPoint")
rm(list=ls())

dataExp<-read.table("F:/Ashley_Projects/Metabolism_Immune/06.MR_relatedRiskModel/05_EstimateScores/dataExp.txt", 
                    header = T, row.names = 1)
icpg<-c("IDO1", "CD274", "HAVCR2", "PDCD1", "CTLA4", "LAG3", "PDCD1LG2")
dataExp_sel<-dataExp[icpg,]
dim(dataExp_sel)
#[1]   7 70
dataExp_sel_t<-t(dataExp_sel)
write.csv(dataExp_sel_t, "ImmuneCheckPointExp.csv", quote = F)


library(corrplot)
M = cor(dataExp_sel_t, use = "everything", method = "spearman")
testRes = cor.mtest(dataExp_sel_t,  use = "everything", method = "spearman", conf.level = 0.95)

pdf("Cor_ICPGs.pdf")
corrplot(M, p.mat = testRes$p, method = 'circle', type = 'full', insig='blank',
         addCoef.col ='black', number.cex = 0.8, order = 'original', diag=T)
dev.off()
#0.05-->*
#0.01-->**
#0.001-->***
write.csv(M, "cor.csv", quote = F)
write.csv(testRes$p, "pvalues.csv", quote = F)

#######################visualization############################
library("tidyr")
library("dplyr")
dataExp_sel_t<-as.data.frame(dataExp_sel_t)
dataExp_sel_t$ID<-row.names(dataExp_sel_t)
dataExp_sel_t<-as.data.frame(dataExp_sel_t)
head(dataExp_sel_t)
#          IDO1    CD274   HAVCR2     PDCD1     CTLA4    LAG3 PDCD1LG2    ID
#L001P 0.121534 0.196102 0.955999 0.4202530 0.3338460 1.85722 0.293138 L001P
#L003P 1.551850 1.281990 4.518670 0.2109940 0.0852604 1.61951 3.249600 L003P
#L004P 2.502590 1.320930 2.216040 0.7473150 0.5183470 2.10494 1.874700 L004P
#L005P 0.476383 1.178510 3.314640 0.0915344 0.0735700 1.60196 1.205040 L005P
#L008P 3.611470 2.413620 4.457170 0.8611140 0.8161360 2.21783 2.451010 L008P
#L010P 2.110280 1.292900 2.576740 0.6100300 0.4080240 1.54968 1.626120 L010P
risk_label<-read.table("F:/Ashley_Projects/Metabolism_Immune/06.MR_relatedRiskModel/05_EstimateScores/risk_class.txt", 
                       header = T)
table(risk_label$riskCLass)
#high  low 
#   7   63 
head(risk_label)
#       futime fustat   DELEC1   MTARC1    SEPT1  MARCHF10 SEPTIN10 MARCHF11 SEPTIN11      ZYX   ZZEF1
#L001_1   19.5      1 0.000000 22.00349 3.194365 0.0000000  11.8199        0  7.49662 19.31170 4.04251
#L001_2   19.5      1 0.110243  5.89019 3.078259 0.0169441  17.1612        0 24.30710 48.29770 4.36698
#L003_3   28.8      0 0.000000 23.53307 0.968030 0.0000000  21.9539        0 36.29560 29.28320 2.27734
#L003_4   28.8      0 0.000000 20.55836 0.727400 0.0000000  21.1538        0 21.12940  8.37228 2.25808
#L004_5   28.7      0 0.000000 10.77281 1.474450 0.0000000  14.1617        0 22.43910 50.57010 2.76070
#L004_6   28.7      0 0.000000 25.45391 1.313468 0.0000000  15.0127        0 10.15950 11.51890 2.69169
#           ZZZ3 riskScore riskCLass
#L001_1  5.01862 0.6973069       low
#L001_2  9.99004 0.4608506       low
#L003_3 17.21810 6.6360938       low
#L003_4 11.12160 6.8012530      high
#L004_5 15.40040 3.5807258       low
#L004_6  7.04392 3.9104384       low

risk_label$sampleName<-c("L001P", "L001T", "L003T", "L003P", "L004T", 
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
merged_data<-merge(dataExp_sel_t, risk_label, by.x = "ID", by.y = "sampleName")
dim(merged_data)
#[1] 70 22
head(merged_data)
write.csv(merged_data, "Merged_data.csv", quote = F, row.names = F)

library("tidyr")
library("dplyr")
mydata<-subset(merged_data, select = c("ID", "IDO1", "CD274", "HAVCR2", "PDCD1", "CTLA4", "LAG3", "PDCD1LG2", "riskCLass"))
head(mydata)
#     ID     IDO1    CD274   HAVCR2     PDCD1     CTLA4      LAG3 PDCD1LG2 riskCLass
#1 L001P 0.121534 0.196102 0.955999 0.4202530 0.3338460 1.8572200 0.293138       low
#2 L001T 1.080120 0.845056 3.807380 0.6332460 0.9786790 0.8067680 0.878279       low
#3 L003P 1.551850 1.281990 4.518670 0.2109940 0.0852604 1.6195100 3.249600      high
#4 L003T 0.531639 0.821252 2.296870 0.0681826 0.1685560 0.0490559 0.603099       low
#5 L004P 2.502590 1.320930 2.216040 0.7473150 0.5183470 2.1049400 1.874700       low
#6 L004T 1.060800 0.748879 2.685160 0.4224910 0.5352890 0.4844910 1.202480       low

tmp_mydata <- mydata %>% mutate(riskCLass = factor(riskCLass))

library("ggplot2")
library("ggunchained")
library("reshape2")
pdf("ImmuneCheckPoints_boxplot.pdf")
tmp_mydata %>% 
  mutate(riskCLass = factor(riskCLass)) %>% 
  pivot_longer(-c(ID,riskCLass), names_to = "Immune_check_point_genes",values_to = "value") %>% 
  ggplot(aes(Immune_check_point_genes,value,fill=riskCLass))+
  geom_split_violin(draw_quantiles = c(0.25, 0.5, 0.75))+
  theme_bw()+
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45,hjust = 1)
  )
dev.off()

cluster1<-filter(tmp_mydata, riskCLass == "high")
cluster2<-filter(tmp_mydata, riskCLass == "low")

res_IDO1<-t.test(cluster1$IDO1,cluster2$IDO1)
res_IDO1$p.value
#[1] 0.9898867

res_CD274<-t.test(cluster1$CD274,cluster2$CD274)
res_CD274$p.value
#[1] 0.8841081

res_HAVCR2<-t.test(cluster1$HAVCR2,cluster2$HAVCR2)
res_HAVCR2$p.value
#[1] 0.9381403

res_PDCD1<-t.test(cluster1$PDCD1,cluster2$PDCD1)
res_PDCD1$p.value
#[1] 0.6808476

res_CTLA4<-t.test(cluster1$CTLA4,cluster2$CTLA4)
res_CTLA4$p.value
#[1] 0.4732059

res_LAG3<-t.test(cluster1$LAG3,cluster2$LAG3)
res_LAG3$p.value
#[1] 0.01629391

res_PDCD1LG2<-t.test(cluster1$PDCD1LG2,cluster2$PDCD1LG2)
res_PDCD1LG2$p.value
#[1] 0.9848008





