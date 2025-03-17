setwd("F:/Ashley_Projects/Metabolism_Immune/06.MR_relatedRiskModel/04_tROC_external/")
rm(list=ls())

library(timeROC)
library(survival)
library(survivalROC)

risk<-read.table("F:/Ashley_Projects/Metabolism_Immune/06.MR_relatedRiskModel/02_ExternalDataSet_GSE124535/risk.txt",header = T, row.names = 1)
time_roc_res <- timeROC(
  T = risk$futime,
  delta = risk$fustat,
  marker = risk$riskScore,
  cause = 1,
  weighting="marginal",
  #  times = c(7, 8, 9),
  times = c(10, 15, 25),
  ROC = TRUE,
  iid = TRUE
)

time_roc_res$AUC
#     t=10      t=15      t=25 
#0.9191176 0.9687500 0.7530120

confint(time_roc_res, level = 0.95)$CI_AUC
#      2.5%  97.5%
#t=10 79.79 104.03
#t=15 91.16 102.59
#t=25 50.86  99.75

time_ROC_df <- data.frame(
  TP_10month = time_roc_res$TP[, 1],
  FP_10month = time_roc_res$FP[, 1],
  TP_15month = time_roc_res$TP[, 2],
  FP_15month = time_roc_res$FP[, 2],
  TP_25month = time_roc_res$TP[, 3],
  FP_25month = time_roc_res$FP[, 3]
)
library(ggplot2)
pdf("ROC_curve_month101525_TCGA.pdf")
ggplot(data = time_ROC_df) +
  geom_line(aes(x = FP_10month, y = TP_10month), size = 1, color = "#BC3C29FF") +
  geom_line(aes(x = FP_15month, y = TP_15month), size = 1, color = "#0072B5FF") +
  geom_line(aes(x = FP_25month, y = TP_25month), size = 1, color = "#E18727FF") +
  geom_abline(slope = 1, intercept = 0, color = "grey", size = 1, linetype = 2) +
  theme_bw() +
  annotate("text",
           x = 0.75, y = 0.25, size = 4.5,
           label = paste0("AUC at 10 months = ", sprintf("%.3f", time_roc_res$AUC[[1]])), color = "#BC3C29FF"
  ) +
  annotate("text",
           x = 0.75, y = 0.15, size = 4.5,
           label = paste0("AUC at 15 months = ", sprintf("%.3f", time_roc_res$AUC[[2]])), color = "#0072B5FF"
  ) +
  annotate("text",
           x = 0.75, y = 0.05, size = 4.5,
           label = paste0("AUC at 25 months = ", sprintf("%.3f", time_roc_res$AUC[[3]])), color = "#E18727FF"
  ) +
  labs(x = "False positive rate", y = "True positive rate") +
  theme(
    axis.text = element_text(face = "bold", size = 11, color = "black"),
    axis.title.x = element_text(face = "bold", size = 14, color = "black", margin = margin(c(15, 0, 0, 0))),
    axis.title.y = element_text(face = "bold", size = 14, color = "black", margin = margin(c(0, 15, 0, 0)))
  )
dev.off()


