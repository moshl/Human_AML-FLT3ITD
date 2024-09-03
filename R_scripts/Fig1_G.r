library(ggplot2)
library(dplyr)

#1G: the distribution of VAF of sITD
setwd("F:/ITD/1G")
ITD_sum_mat <- read.csv("ITD_VAF_PCR_130.csv")

ITD_sum_mat$Group <- "Major"
sITD <- ITD_sum_mat$ID

> max(ITD_sum_mat$VAF)
[1] 95.52
> min(ITD_sum_mat$VAF)
[1] 0.83
> median(ITD_sum_mat$VAF)
[1] 37.255
> mean(ITD_sum_mat$VAF)
[1] 35.74715

list_order_s <- ITD_sum_mat[order(ITD_sum_mat$VAF,decreasing = T),]  

ITD_sum_mat$Order <- NA
for (i in 1:nrow(ITD_sum_mat)){
  sampleID_i <- ITD_sum_mat[i,"ID"]
  ITD_sum_mat[i,"Order"] = which(list_order_s == sampleID_i)  
}

p <- ITD_sum_mat %>%
  ggplot(aes(x=Order,y=VAF,fill=Group)) +
  geom_point(aes(fill=Group), shape = 21, color="white",size = 2) + 
  scale_fill_manual(values = "#009DC1") +   
  geom_segment(aes(x=Order,xend=Order,y=0,yend=VAF),linetype=1,col="gray30",alpha=0.1,linewidth=0.1) +   
  ylab("FLT3-ITD VAF(%)") +
  scale_x_continuous(breaks=c(1:length(sITD)),labels=1:130) +  
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, size = 5, vjust = .5))
p

ggsave("1G_sITD_VAF_PCR.pdf",
       units = "in", width = 10,height = 3)  