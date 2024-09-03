library(ggplot2) 
library(dplyr)
library(ggsignif) 
library(ggdist) 
library(tidyr)
library("ggpubr")
library("ggthemes")
library("reshape2")

setwd("D:/ITD/1J")

#1J: the number of typical and atypical in sITD and mITD 
df<-read.table("typical_atypical.txt",head=TRUE,sep="\t")
data<-data.frame(df)
Custom.color <- c("#F8766D", "#009DC1")
pdf("1J_typical_atypical.pdf",width = 4,height = 5)
df %>%
  mutate(Group = factor(Group,levels=c("sITD","mITD")))%>%  
  ggplot(aes(x = Group, y = Patients, fill=Type)) +
  geom_bar(position="stack",stat="identity", width=0.5)+ 
  scale_fill_manual(values = Custom.color)+   
  scale_color_manual(values = Custom.color)+  
  ylim(0,150)+ 
  ylab("Number of patients") +   #设置Y轴标题
  theme_classic() 
dev.off()

X-squared = 19.334, df = 1, p-value = 1.098e-05

# chisq.test
df<-read.table("typical_atypical_stat.txt",head=TRUE,sep="\t")
tb <- table(df[,2:3])
chisq.test(tb, correct = FALSE)

          Atypical
ITD_group   0  1
  Multiple 30 55
  Single   91 48
  
X-squared = 19.334, df = 1, p-value = 1.098e-05