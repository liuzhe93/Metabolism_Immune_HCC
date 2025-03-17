setwd("F:/Ashley_Projects/Metabolism_Immune/06.MR_relatedRiskModel/05_EstimateScores")
rm(list=ls())

library("estimate")

mydata<-read.csv("scores_estimate.csv")
head(mydata)
#      X StromalScore ImmuneScore ESTIMATEScore sample_name TumourPurity
#1 L001P   -457.52806    774.5618      317.0338        L001    0.7951605
#2 L003P   -196.62549   1265.7836     1069.1581        L003    0.7235070
#3 L004P   -215.87662   1101.6826      885.8060        L004    0.7418219
#4 L005P   -169.44134    901.3164      731.8750        L005    0.7567835
#5 L008P     21.71349   1597.9419     1619.6554        L008    0.6654246
#6 L010P    -11.58435   1295.7245     1284.1402        L010    0.7013661

risk<-read.table("F:/Ashley_Projects/Metabolism_Immune/06.MR_relatedRiskModel/02_ExternalDataSet_GSE124535/risk.txt")
head(risk)
#       futime fustat   DELEC1   MTARC1    SEPT1  MARCHF10 SEPTIN10 MARCHF11 SEPTIN11      ZYX   ZZEF1
#L001_1   19.5      1 0.000000 22.00349 3.194365 0.0000000  11.8199        0  7.49662 19.31170 4.04251
#L001_2   19.5      1 0.110243  5.89019 3.078259 0.0169441  17.1612        0 24.30710 48.29770 4.36698
#L003_3   28.8      0 0.000000 23.53307 0.968030 0.0000000  21.9539        0 36.29560 29.28320 2.27734
#L003_4   28.8      0 0.000000 20.55836 0.727400 0.0000000  21.1538        0 21.12940  8.37228 2.25808
#L004_5   28.7      0 0.000000 10.77281 1.474450 0.0000000  14.1617        0 22.43910 50.57010 2.76070
#L004_6   28.7      0 0.000000 25.45391 1.313468 0.0000000  15.0127        0 10.15950 11.51890 2.69169
#           ZZZ3 riskScore
#L001_1  5.01862 0.6973069
#L001_2  9.99004 0.4608506
#L003_3 17.21810 6.6360938
#L003_4 11.12160 6.8012530
#L004_5 15.40040 3.5807258
#L004_6  7.04392 3.9104384

risk$riskCLass<-ifelse(risk$riskScore>6.636094,"high","low")
write.table(risk,"risk_class.txt",sep="\t",row.names = T, quote = F)
table(risk$riskCLass)
#high  low 
#   7   63 
head(risk)
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
risk$sampleName<-c("L001P", "L001T", "L003T", "L003P", "L004T", 
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
merged_data<-merge(mydata, risk, by.x = "X", by.y = "sampleName")
dim(merged_data)
#[1] 70 20
head(merged_data)
#X       StromalScore ImmuneScore ESTIMATEScore sample_name TumourPurity futime fustat   DELEC1   MTARC1
#1 L001P    -457.5281    774.5618      317.0338        L001    0.7951605   19.5      1 0.000000 22.00349
#2 L001T     561.7740   1178.8698     1740.6438        L001    0.6520634   19.5      1 0.110243  5.89019
#3 L003P    -196.6255   1265.7836     1069.1581        L003    0.7235070   28.8      0 0.000000 20.55836
#4 L003T    -517.4538    285.2957     -232.1581        L003    0.8414101   28.8      0 0.000000 23.53307
#5 L004P    -215.8766   1101.6826      885.8060        L004    0.7418219   28.7      0 0.000000 25.45391
#6 L004T    -240.3811    744.3294      503.9484        L004    0.7782256   28.7      0 0.000000 10.77281
#     SEPT1  MARCHF10 SEPTIN10 MARCHF11 SEPTIN11      ZYX   ZZEF1     ZZZ3 riskScore riskCLass
#1 3.194365 0.0000000  11.8199        0  7.49662 19.31170 4.04251  5.01862 0.6973069       low
#2 3.078259 0.0169441  17.1612        0 24.30710 48.29770 4.36698  9.99004 0.4608506       low
#3 0.727400 0.0000000  21.1538        0 21.12940  8.37228 2.25808 11.12160 6.8012530      high
#4 0.968030 0.0000000  21.9539        0 36.29560 29.28320 2.27734 17.21810 6.6360938       low
#5 1.313468 0.0000000  15.0127        0 10.15950 11.51890 2.69169  7.04392 3.9104384       low
write.csv(merged_data, "Merged_data.csv", quote = F, row.names = F)

library("tidyverse")
library("gapminder")
library("ggsci")
library("ggprism")
library("rstatix")
library("ggpubr")
library("reshape2")
library("knitr")
library("rstatix")


###############################################################################################
StromalScore<-merged_data[, c("X", "riskCLass", "StromalScore")]
StromalScore$cluster<-as.factor(StromalScore$riskCLass)
df<-StromalScore
df_p_val1 <- df %>% 
  wilcox_test(StromalScore  ~ cluster) %>%
  adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  add_significance(p.col = "p.adj") %>% 
  add_xy_position(x = "Risk Level", dodge = 0.8) 
pdf("StromalScore_boxplot.pdf")
p <- ggboxplot(df, x = "cluster", y = "StromalScore",fill = "cluster", 
               color = 'cluster', palette = c("blue", "red","green"), add = "boxplot", 
               add.params = list(fill = "white")) + 
  labs(x = "Risk Level", y = "Stromal Score", fill = "cluster") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p + stat_compare_means(aes(group = cluster), label = "p.signif", paired = FALSE)
dev.off()


###############################################################################################
ImmuneScore<-merged_data[, c("X", "riskCLass", "ImmuneScore")]
ImmuneScore$cluster<-as.factor(ImmuneScore$riskCLass)
df<-ImmuneScore
df_p_val1 <- df %>% 
  wilcox_test(ImmuneScore  ~ cluster) %>%
  adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  add_significance(p.col = "p.adj") %>% 
  add_xy_position(x = "Risk Level", dodge = 0.8) 
pdf("ImmuneScore_boxplot.pdf")
p <- ggboxplot(df, x = "cluster", y = "ImmuneScore",fill = "cluster", 
               color = 'cluster', palette = c("blue", "red","green"), add = "boxplot", 
               add.params = list(fill = "white")) + 
  labs(x = "Risk Level", y = "Immune Score", fill = "cluster") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p + stat_compare_means(aes(group = cluster), label = "p.signif", paired = FALSE)
dev.off()


###############################################################################################
TumourPurity<-merged_data[, c("X", "riskCLass", "TumourPurity")]
TumourPurity$cluster<-as.factor(TumourPurity$riskCLass)
df<-TumourPurity
df_p_val1 <- df %>% 
  wilcox_test(TumourPurity  ~ cluster) %>%
  adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  add_significance(p.col = "p.adj") %>% 
  add_xy_position(x = "Risk Level", dodge = 0.8) 
pdf("TumourPurity_boxplot.pdf")
p <- ggboxplot(df, x = "cluster", y = "TumourPurity",fill = "cluster", 
               color = 'cluster', palette = c("blue", "red","green"), add = "boxplot", 
               add.params = list(fill = "white")) + 
  labs(x = "Risk Level", y = "Tumour Purity", fill = "cluster") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p + stat_compare_means(aes(group = cluster), label = "p.signif", paired = FALSE)
dev.off()



