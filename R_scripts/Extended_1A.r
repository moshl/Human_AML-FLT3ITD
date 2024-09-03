library(ggplot2) 
library(ggsignif) 
library(ggdist) 
library(tidyr)
library("ggpubr")
library("ggthemes")
library("reshape2")

# Extended_1A: the proportion of MM indel in each cancer type and MM gene
setwd("D:/ITD/Extended_1A")
info<-read.table("indel_MM_percent.txt",head=TRUE,sep="\t")
data<-data.frame(info)  

Custom.color <- c("#009DC1", "#F8766D")
pdf("Extended_1A_MM_Indel_percent.pdf",width = 4,height = 5)
data %>%
  mutate(Group = factor(Group,levels=c("low","high")))%>%  
  ggboxplot(x = "Group",y = "percent",
            color = "Group",   
            ylab = "MM Indel Percent",    
            add="jitter") + 
  scale_fill_manual(values = Custom.color)+   
  scale_color_manual(values = Custom.color)+  
  theme_classic()
dev.off()

low <- data[which(data$Group=="low"),"percent"]
high <- data[which(data$Group=="high"),"percent"]
wilcox.test(low,high,paired=F,alternative = "less") 
p-value = 0.02836
