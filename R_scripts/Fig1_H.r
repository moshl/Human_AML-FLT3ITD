#library(forcats)
library(ggplot2)
library(dplyr)

setwd("D:/ITD/1H")
ITD_sum_mat <- read.csv("ITD_VAF_PCR_81.csv")

# 多ITD中VAF最高的是major，其余都是minor
ITD_sum_mat$Group <- "Minor"
mITD <- unique(ITD_sum_mat$ID)

for (i in 1:length(mITD)){
  mITD_seq <- ITD_sum_mat[which(ITD_sum_mat$ID == mITD[i]),]   
  nrow <- which(ITD_sum_mat$VAF == max(mITD_seq$VAF))     #每个mITD最大的VAF对应的行，可能有重复
  ITD_sum_mat[nrow, "Group"] <- "Major"
}

major_ITD <- ITD_sum_mat[which(ITD_sum_mat$Group == "Major"), ]   
nrow(major_ITD)   
dup <- major_ITD[duplicated(major_ITD$ID),]     
ITD_sum_mat[c(67,174),"Group"] <- "Minor"

table(ITD_sum_mat$Group)
Major Minor 
85   118

# diff of max_VAF and minor_VAF
diff <- as.data.frame(mITD)
diff$max_min <- NA
for (i in 1:length(mITD)){
  mITD_seq <- ITD_sum_mat[which(ITD_sum_mat$ID == mITD[i]),]   
  diff[i,2] <- max(mITD_seq$VAF)-min(mITD_seq$VAF)  
}

#PCR
max(diff[,2])
[1] 53.87
> min(diff[which(diff$max_min >0),2])
[1] 0.04
> median(diff[which(diff$max_min >0),2])
[1] 20.06
> mean(diff[which(diff$max_min >0),2])
[1] 22.81167

major_ITD <- ITD_sum_mat[which(ITD_sum_mat$Group == "Major"), ]  
list_order_m <- major_ITD[order(major_ITD$VAF,decreasing = T),]  

ITD_sum_mat$Order <- NA
for (i in 1:nrow(ITD_sum_mat)){
  sampleID_i <- ITD_sum_mat[i,"ID"]
  ITD_sum_mat[i,"Order"] = which(list_order_m == sampleID_i)    
}

Custom.color <- c("#F8766D", "#009DC1")
p <- ITD_sum_mat %>%
  ggplot(aes(x=Order,y=VAF,fill=Group)) +
  geom_point(aes(fill=Group), shape = 21, color="white",size = 3) + 
  scale_fill_manual(values = Custom.color) + 
  geom_segment(aes(x=Order,xend=Order,y=0,yend=VAF),linetype=1,col="gray30",alpha=0.1,linewidth=0.1) +   
  ylab("FLT3-ITD VAF(%)") +
  scale_x_continuous(breaks=c(1:length(mITD)),labels=1:85) +  
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, size = 5, vjust = .5))
p

ggsave("1H_MMITD_VAF_PCR.pdf",
       units = "in", width = 10,height = 3)  
