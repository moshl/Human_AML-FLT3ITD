library(pheatmap)
library(ggplot2)
library(reshape2)
library(ggforce)
library("ggsci")

setwd("D:/ITD/3CD")

# 3C: fraction genome alterled in sITD and mITD
info<-read.table("3C_D0_CNV_36_genome.txt",head=TRUE,sep="\t")
data<-data.frame(info)  

Custom.color <- c("#71c9ce", "#ffaaa6")
pdf("3C_D0_36_boxplot_dot_genome.pdf",width = 4,height = 5)
data %>%
  mutate(Group = factor(Group,levels=c("sITD","mITD")))%>%  
  ggboxplot(x = "Group",y = "genome_fra",
            color = "Group",   
            ylab = "Fraction Genome Altered (%)",    
            add="jitter") + 
  scale_fill_manual(values = Custom.color)+   
  scale_color_manual(values = Custom.color)+  
  theme_classic()
dev.off()


# 3D: fraction genome alterled in sITD_CR, sITD_NR, mITD_CR and mITD_NR
info<-read.table("3D_D0_CNV_36_CR_NR.txt",head=TRUE,sep="\t")
data<-data.frame(info)  

Custom.color <- c("#71c9ce", "#ffaaa6", "#ffde7d", "#aa96da")
pdf("3D_D0_CNV_36_CR_NR_boxplot_dot.pdf",width = 6,height = 5)
data %>%
  mutate(Group = factor(Group,levels=c("sITD_CR","sITD_NR","mITD_CR","mITD_NR")))%>%  
  ggboxplot(x = "Group",y = "genome_fra",
            color = "Group",   
            ylab = "Fraction Genome Altered (%)",    
            add="jitter") +  
  scale_fill_manual(values = Custom.color)+   
  scale_color_manual(values = Custom.color)+  
  theme_classic()
dev.off()