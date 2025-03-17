setwd("F:/Ashley_Projects/Metabolism_Immune/11.Immunotherapy/01.TIDE/Responder")
rm(list=ls())



#####################################MC#############################################################
MC_res<-read.csv("F:/Ashley_Projects/Metabolism_Immune/11.Immunotherapy/01.TIDE/MC/TIDE_res.csv")
head(MC_res)
library(stringr)
MC_res$Risk<-ifelse(str_sub(MC_res$Patient, 1, 2) == "C1", "C1", "C2")
MC_res<-subset(MC_res, select = c("Risk", "Responder"))
MC_SP<-table(MC_res$Risk, MC_res$Responder)
MC_NO<-addmargins(MC_SP)
MC_NO_matrix<-as.matrix(MC_NO)[1:2,1:2]
MC_NO_df<-as.data.frame.array(MC_NO_matrix)
MC_NO_df
MC_NO_df$MC<-row.names(MC_NO_df)
library(tidyverse)
MC_count <- MC_NO_df %>% pivot_longer(cols=c(False:True),
                                      names_to = 'Responder',
                                      values_to = 'freq')
MC_count

library(plyr)
MC_count_transfer = ddply(MC_count,'MC',transform,percent_con=freq/sum(freq)*100)
MC_count_transfer

pdf("MC_percentage.pdf")
p3 <- ggplot(MC_count_transfer,aes(x=MC,y=percent_con,fill=Responder))+
  geom_bar(stat = 'identity',width = 0.5,colour='black')
p4 <- p3+labs(x='MC',y='Percentage')+
  theme(axis.title = element_text(size=12),
        axis.text = element_text(size=11))+
  scale_y_continuous(breaks=seq(0,100,25),
                     labels=c('0','25%','50%','75%','100%'))

p4
dev.off()





#####################################GC#############################################################
GC_res<-read.csv("F:/Ashley_Projects/Metabolism_Immune/11.Immunotherapy/01.TIDE/GC/TIDE_res.csv")
head(GC_res)
library(stringr)
GC_res$Risk<-GC_res$Patient
for(i in 1:nrow(GC_res)){
  GC_res$Risk[i]<-substr(GC_res$Patient[i], 1, 2)
}

GC_res<-subset(GC_res, select = c("Risk", "Responder"))
GC_SP<-table(GC_res$Risk, GC_res$Responder)
GC_NO<-addmargins(GC_SP)
GC_NO_matrix<-as.matrix(GC_NO)[1:3,1:2]
GC_NO_df<-as.data.frame.array(GC_NO_matrix)
GC_NO_df
GC_NO_df$GC<-row.names(GC_NO_df)
library(tidyverse)
GC_count <- GC_NO_df %>% pivot_longer(cols=c(False:True),
                               names_to = 'Responder',
                               values_to = 'freq')
GC_count

library(plyr)
GC_count_transfer = ddply(GC_count,'GC',transform,percent_con=freq/sum(freq)*100)
GC_count_transfer

pdf("GC_percentage.pdf")
p3 <- ggplot(GC_count_transfer,aes(x=GC,y=percent_con,fill=Responder))+
  geom_bar(stat = 'identity',width = 0.5,colour='black')
p4 <- p3+labs(x='GC',y='Percentage')+
  theme(axis.title = element_text(size=12),
        axis.text = element_text(size=11))+
  scale_y_continuous(breaks=seq(0,100,25),
                     labels=c('0','25%','50%','75%','100%'))

p4
dev.off()



