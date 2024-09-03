setwd("F:/ITD/1D")

# the diversity of FLT3-ITD sequences
ITD <- read.table("ITD_seq.txt", header=T, sep = "\t") 
length(unique(ITD$ITD))
[1] 820

unique_ITD <- as.data.frame(unique(ITD$ITD))   
unique_ITD$count <- 0     
for (i in 1:nrow(ITD)){
  n <- which(unique_ITD == ITD[i,])   
  unique_ITD[n,2] <- unique_ITD[n,2]+1   
}

sum(unique_ITD$count)
[1] 1116
unique_ITD$frequency <- round(unique_ITD$count/1116,4)
unique_ITD <- unique_ITD[order(unique_ITD$count, decreasing = T),]
colnames(unique_ITD)[1] <- "sequence"
unique_ITD$length <- nchar(unique_ITD$sequence)

unique_ITD$ID <- NA
for (i in 1:nrow(unique_ITD)){
  unique_ITD[i,"ID"] <- paste0("ITD",i)
}

write.csv(unique_ITD, "occurrence of each ITD sequence.csv", quote = F, row.names = F)

# 绘制多样性
library(ggplot2)

df <- read.table("input.txt",header=TRUE,sep="\t") 

Custom.color <- c("#F8766D",  "#F8B5B0", "#F8E3E2")
p1 <- ggplot(df,aes(x=factor(ID,levels=as.character(ID)),
                    y=Frequency,fill=Number))+
  geom_bar(position="stack",stat="identity")+ # , width=3
  scale_fill_manual(values = Custom.color)+
  xlab(" ") +   
  ylab("Frequency") +   
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45))   # , size = 10

pdf("Fig1D_occurrence of each ITD sequence.pdf",width = 10,height = 5)
print(p1)
dev.off()


ratio <- c(6.54, 33.48, 59.98)
group <- c(">10", "2-10", "1")
Custom.color <- c("#F8766D",  "#F8B5B0", "#F8E3E2")
pdf("Fig1D_pie.pdf",width = 4,height = 4)
pie(ratio, labels=group,radius = 1.0,clockwise=T,col=Custom.color,main = "Frequency")
dev.off()

