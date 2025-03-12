setwd("F:/Ashley_Projects/Metabolism_Immune/01.IdentificationOfKeyMetabolicPathways/08_EstimateScores")
rm(list=ls())

data_merged_sel<-read.csv("data_merged.csv", row.names = 1)
dim(data_merged_sel)
#[1] 390   7

#calculate the correlation matrix
cor(x = data_merged_sel, use = "everything", method = "spearman")
library("psych")
res<-corr.test(data_merged_sel, use = "pairwise", method = "spearman", adjust = "holm", alpha=.05, minlength=5, ci=TRUE, normal = F)
#https://www.math.pku.edu.cn/teachers/lidf/course/mvr/mvrnotes/html/_mvrnotes/multcomp.html
print(res,digits=6, short=F)
#输出相关性和P值
res$r
res$ci

#visualization
#https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html
library(corrplot)
M = cor(data_merged_sel, use = "everything", method = "spearman")
testRes = cor.mtest(data_merged_sel,  use = "everything", method = "spearman", conf.level = 0.95)

pdf("EstimateScore_cor.pdf")
corrplot(M, p.mat = testRes$p, method = 'circle', type = 'full', insig='blank',
         addCoef.col ='black', number.cex = 0.8, order = 'original', diag=T)
dev.off()
#0.05-->*
#0.01-->**
#0.001-->***
write.csv(M, "cor.csv", quote = F)
write.csv(testRes$p, "pvalues.csv", quote = F)



