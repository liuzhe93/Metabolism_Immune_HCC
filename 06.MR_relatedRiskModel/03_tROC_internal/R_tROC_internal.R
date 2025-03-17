setwd("F:/Ashley_Projects/Metabolism_Immune/06.MR_relatedRiskModel/03_tROC_internal")
rm(list=ls())

library(timeROC)
library(survival)
library(survivalROC)

risk<-read.table("F:/Ashley_Projects/Metabolism_Immune/06.MR_relatedRiskModel/01_RiskModel/risk.txt",header = T, row.names = 1)
time_roc_res <- timeROC(
  T = risk$futime,
  delta = risk$fustat,
  marker = risk$riskScore,
  cause = 1,
  weighting="marginal",
#  times = c(7, 8, 9),
  times = c(1,3,5),
  ROC = TRUE,
  iid = TRUE
)

time_roc_res$AUC
#      t=1       t=3       t=5 
#0.7580786 0.7504173 0.6461774

confint(time_roc_res, level = 0.95)$CI_AUC
#     2.5% 97.5%
#t=1 68.97 82.64
#t=3 68.09 81.99
#t=5 55.51 73.73

time_ROC_df <- data.frame(
  TP_1year = time_roc_res$TP[, 1],
  FP_1year = time_roc_res$FP[, 1],
  TP_3year = time_roc_res$TP[, 2],
  FP_3year = time_roc_res$FP[, 2],
  TP_5year = time_roc_res$TP[, 3],
  FP_5year = time_roc_res$FP[, 3]
)
library(ggplot2)
pdf("ROC_curve_year135_TCGA.pdf")
ggplot(data = time_ROC_df) +
  geom_line(aes(x = FP_1year, y = TP_1year), size = 1, color = "#BC3C29FF") +
  geom_line(aes(x = FP_3year, y = TP_3year), size = 1, color = "#0072B5FF") +
  geom_line(aes(x = FP_5year, y = TP_5year), size = 1, color = "#E18727FF") +
  geom_abline(slope = 1, intercept = 0, color = "grey", size = 1, linetype = 2) +
  theme_bw() +
  annotate("text",
           x = 0.75, y = 0.25, size = 4.5,
           label = paste0("AUC at 1 years = ", sprintf("%.3f", time_roc_res$AUC[[1]])), color = "#BC3C29FF"
  ) +
  annotate("text",
           x = 0.75, y = 0.15, size = 4.5,
           label = paste0("AUC at 3 years = ", sprintf("%.3f", time_roc_res$AUC[[2]])), color = "#0072B5FF"
  ) +
  annotate("text",
           x = 0.75, y = 0.05, size = 4.5,
           label = paste0("AUC at 5 years = ", sprintf("%.3f", time_roc_res$AUC[[3]])), color = "#E18727FF"
  ) +
  labs(x = "False positive rate", y = "True positive rate") +
  theme(
    axis.text = element_text(face = "bold", size = 11, color = "black"),
    axis.title.x = element_text(face = "bold", size = 14, color = "black", margin = margin(c(15, 0, 0, 0))),
    axis.title.y = element_text(face = "bold", size = 14, color = "black", margin = margin(c(0, 15, 0, 0)))
  )
dev.off()
