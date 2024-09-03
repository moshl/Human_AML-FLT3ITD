library(ggplot2) 
library(ggsignif) 
library(ggdist) 
library(tidyr)
library("ggpubr")
library("ggthemes")
library("reshape2")

setwd("D:/ITD/2B_2H")

# 2B: clone number of sITD and mITD
info<-read.table("2D_clone_number_0507.txt",head=TRUE,sep="\t")
data<-data.frame(info)  

Custom.color <- c("#71c9ce", "#ffaaa6")
pdf("2B_clone_number.pdf",width = 3,height = 5)
data %>%
  mutate(Group = factor(Group,levels=c("sITD","mITD")))%>%  
  ggboxplot(x = "Group",y = "Clusters",
            color = "Group",   
            ylab = "Number of clones",    
            add="jitter") + 	
  scale_fill_manual(values = Custom.color)+   
  scale_color_manual(values = Custom.color)+  
  theme_classic()+
  geom_signif(comparisons = list(c("sITD", "mITD")),step_increase = .1,
              map_signif_level = TRUE,vjust = 0.5,hjust= 0)
dev.off()

sITD <- data[which(data$Group=="sITD"),"Clusters"]
mITD <- data[which(data$Group=="mITD"),"Clusters"]
wilcox.test(sITD,mITD,paired=F,alternative = "two.sided")  
p-value = 0.01282


# 2H: TPM of key genes of SSA and alt-NHEJ
info<-read.table("2H_sITD_VS_mITD_TPM.txt",head=TRUE,sep="\t")
data<-data.frame(info)  

data <- melt(data, id.vars=c("ID", "Group"), 
             variable.name="Gene", 
             value.name="exp")

Custom.color <- c("#F8766D","#009DC1")
pdf("2H_sITD_VS_mITD.pdf",width = 10,height = 4)
data %>%
  mutate(Group = factor(Group,levels=c("sITD","mITD")))%>%  
  ggboxplot(x = "Gene",y = "exp",
            color = "Group",  
            xlab = "Gene",ylab = "TPM", 
            add="jitter") +
  theme_classic() 
dev.off()

p <- NA
for (i in 3:6) {
  sITD <- info[which(info$Group=="sITD"),i]
  mITD <- info[which(info$Group=="mITD"),i]
  test <- wilcox.test(sITD,mITD,paired=F,alternative = "two.sided")  #c("two.sided", "less", "greater")
  p[i-2] <- test$p.value
}
> p
[1] 0.032984920 0.003323213 0.024371120 0.017125176