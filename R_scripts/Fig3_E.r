library(pheatmap)
library(ggplot2)
library(reshape2)
library(ggforce)
library("ggsci")


setwd("D:/ITD/3E")

# left: cluster of segment number and score of CNV signature
df <- read.table("seg_type_length_CN_D0_36_normal.txt",header=TRUE,sep="\t") 
rownames(df)<-df[,1]
df<-df[,-1]
cohort <- colnames(df)

group<-read.table("group.txt",head=TRUE,sep="\t")
Segment_number <- c()
ITD_group <- c()
Karyotype <- c()
Resistant_group <- c()
CN1 <- c()
CN5 <- c()
CN9 <- c()
CN21 <- c()

for (i in 1:length(cohort)) {  
  ITD_group[i] <- group[which(group[,1]==cohort[i]), 4]
  Resistant_group[i] <- group[which(group[,1]==cohort[i]), 6]
  Karyotype[i] <- group[which(group[,1]==cohort[i]), 5]
  Segment_number[i] <- group[which(group[,1]==cohort[i]), 3]
  CN1[i] <- group[which(group[,1]==cohort[i]), 7]
  CN5[i] <- group[which(group[,1]==cohort[i]), 8]
  CN9[i] <- group[which(group[,1]==cohort[i]), 9]
  CN21[i] <- group[which(group[,1]==cohort[i]), 10]  
}

Custom.color <- c("#71c9ce",  "#bbded6", "#ffaaa6", "#ffde7d", "#aa96da", "#fcbad3")
ann_colors = list(
  Segment_number=c(colorRampPalette(colors = c("white","firebrick3"))(100)),
  ITD_group = c("sITD" = "#71c9ce", "mITD" = "#ffaaa6"),
  Karyotype=c("NK" = "#ffde7d", "noNK" = "#aa96da"),
  Resistant_group=c("CR_CR" = "#bbded6", "PR_CR" = "#fae3d9", "NR_NR" = "#ffb6b9"),
  CN1=c(colorRampPalette(colors = c("white","firebrick3"))(100)),
  CN5=c(colorRampPalette(colors = c("white","firebrick3"))(100))
) 

################################ pheatmap ################################
annotation_col = data.frame( 
  ITD_group = factor(ITD_group), Resistant_group = factor(Resistant_group),
  Karyotype= factor(Karyotype), Segment_number = Segment_number, CN1=CN1, CN5=CN5) 
rownames(annotation_col) = colnames(df)

bk = unique(c(seq(-1,1, length=100))) 
callback = function(hc, mat){      
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}

pdf("Fig3E_pheatmap_36_D0_seg_length_type_CN_normal.pdf", width = 7, height = 4) 
p<-pheatmap(df, breaks=bk,scale = "row", # kmeans_k = 2，把列聚成两类
            annotation_col = annotation_col, annotation_names_col=TRUE,
            #annotation_row = annotation_row, annotation_names_row=TRUE,
            cluster_rows = T, cluster_cols =T,
            clustering_distance_rows="euclidean",
            clustering_distance_cols = "euclidean",
            clustering_callback = callback,
            annotation_colors = ann_colors,
            show_rownames = T, show_colnames=T, #angle_col=45,
            fontsize=5, fontsize_row=4, fontsize_col=4, #border_color="white",
            #color= c(colorRampPalette(c("navy","#7986CB"))(70),
            #         colorRampPalette(c("#7986CB",'#E5E5E5','#EF5350'))(50),
            #         colorRampPalette(c("#EF5350","#B71C1C"))(30)),
            width=14, height=15, clustering_method='ward.D')
print(p)
dev.off()



# right: SNV segment
data<-read.csv("all_rawsegs_合并后_data_D0_36.csv",head=TRUE,sep=",")
new_data <- data
chr_length<-read.table("chr_length.txt",head=TRUE,sep="\t")
group<-read.table("group.txt",head=TRUE,sep="\t")
samples <- unique(new_data[,"sampleID"])   
n_sample <- length(samples)     

df <- new_data   
chr_fraction <- matrix(NA,nrow(chr_length)+1,2)
chr_sum <- sum(chr_length[,2])
i=1
for (i in 1:nrow(chr_length)){
  chr_fraction[i+1,1] <- chr_length[i,2]/chr_sum
  chr_fraction[i+1,2] <- sum(chr_length[1:i,2])/chr_sum
}
chr_fraction[1,1:2] <- 0

samples <- group[,"sampleID"]   

x_y <- matrix(NA,nrow(df),4)
colnames(x_y) <- c("x0","x1","y0","y1")
for (i in 1:nrow(df)){
  chr <- as.numeric(df[i,"chromosome"])
  chr_len <- df[i,"chr_len"]
  start_pos <- df[i,"start.pos"]
  end_pos <- df[i,"end.pos"]
  x_y[i,"x0"] <- chr_fraction[chr,2]+start_pos/chr_len*chr_fraction[chr+1,1]  # "start.pos"
  x_y[i,"x1"] <- chr_fraction[chr,2]+end_pos/chr_len*chr_fraction[chr+1,1]   # "end.pos"
  sampleID_i <- df[i,"sampleID2"]
  if (sampleID_i %in% samples == TRUE){
    x_y[i,3:4] = which(samples== sampleID_i)   
  } 
}

df<-cbind(df,x_y)   

# LOSS，GAIN，UPD
new_type <- matrix(NA,nrow(df),1)
df <- cbind(df, new_type)
df[which(df$type == "NLOH"),"new_type"] <- "LOSS"    #0,1
df[which(df$type == "DLOH"),"new_type"] <- "GAIN"    #0,>2 
df[which(df$type == "HD"),"new_type"] <- "LOSS"     #0,0
df[which(df$type == "GAIN"),"new_type"] <- "GAIN"   
df[which(df$type == "UPD"),"new_type"] <- "UPD"  #UPD
df[which(df$type == "normal"),"new_type"] <- "normal"

pdf("3E_36_D0_seg.pdf",width=6,height=5)   
ggplot(df, aes(x=x0, y=y0, xend=x1, yend=y1,
               col=factor(new_type),
               main = "D0_36"))+
  geom_segment(aes(linewidth=1),linewidth = 3)+   
  scale_color_manual(values = c("#0074b3", "#db6968", "white",  "#459943"))+   
  scale_x_continuous(breaks=chr_fraction[2:23,2],labels=c(1:22))+  
  scale_y_continuous(breaks=c(1:36),labels=samples)+  #  
  labs(x = "chr", y = "Samples")+
  theme_bw()+   
  theme_classic()+   
  theme(axis.text.y = element_text(size = 5)     
        ,axis.text.x = element_text(size = 5))+    
  theme(strip.background = element_rect(fill = c("white"), 
                                        linetype = "solid"), 
        strip.text = element_text(color = "black", 
                                  face = "bold", 
                                  size = 8)) 
dev.off()
