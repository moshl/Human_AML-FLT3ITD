######
library(dplyr)

setwd("F:/ITD/1BC")
dat_dir="F:/ITD/1BC/"
Mut<-read.delim(paste0(dat_dir,"Pan.TotalGene.Freq.txt"),header=F,sep = "\t",stringsAsFactors = F,skip=1)  
 
#1 the proportion of indels in 35 cancer types and 49 MM genes
colnames(Mut)<-c('Cancer','Gene','Indel','SNV','Pts_with_Indel','Pts_with_SNV','Total_Pts','Indel_len')
Mut=Mut[-which(Mut$Gene=="Gene"),] 
Mut[,3:7]=apply(Mut[,3:7],2,as.numeric)
Mut<-Mut[which(Mut$Total_Pts>0),] 
Mut<-Mut[which(Mut$Gene!="."),]  

Nature_MM_merge <- c("PIK3CA","EGFR","CTNNB1","ERBB2","MTOR","ERBB3","NRAS","KRAS","KIT","PDGFRA","FGFR2","FGFR3","PPP2R1B","SPOP","FLT3","JAK2","CARD11","JAK3","NOTCH1","CHD4","APC","TP53","PTEN","MAP3K1","PIK3R1","ARID1B","CDK12","FBXW7","CIC","NF1","CDKN1B","BAP1","PTPRT","ARID2","SOX9","INPPL1","ERBB4","TSC2","CDKN2A","SMAD4","TBX3","NKX2-1","B2M","RUNX1","EPHA7","FOXA1","RB1","DICER1","HLA-B")

head(Mut)[1:2,]
  Cancer Gene Indel SNV Pts_with_Indel Pts_with_SNV Total_Pts Indel_len
1    ACC DDC8     0   1              0            1         1         ,
3    ACC AACS     0   2              0            1         1         ,
unique(Mut$Cancer)  
 [1] "ACC"      "BLCA"     "BRCA"     "CESC"     "CHOL"     "COADREAD" "DLBC"    
 [8] "ESCA"     "GBM"      "HMALL"    "HMAML"    "HMCLL2"   "HMCLL"    "HMDLBC"  
[15] "HMMM"     "HMMPN"    "HNSC"     "KICH"     "KIRC"     "KIRP"     "LAML"    
[22] "LIHC"     "LUAD"     "LUSC"     "MESO"     "OV"       "PAAD"     "PCPG"    
[29] "PRAD"     "SARC"     "SKCM"     "STAD"     "TGCT"     "THCA"     "THYM"    
[36] "UCEC"     "UCS"      "UVM" 
Mut$Cancer=gsub("HM","",Mut$Cancer)
Mut$Cancer=gsub("LAML","AML",Mut$Cancer)   #merge TCGA-AML and HMAML
Mut$Cancer=gsub("CLL2","CLL",Mut$Cancer)   #merge HMCLL and HMCLL2
unique(Mut$Cancer)  
 [1] "ACC"      "BLCA"     "BRCA"     "CESC"     "CHOL"     "COADREAD" "DLBC"    
 [8] "ESCA"     "GBM"      "ALL"      "AML"      "CLL"      "MM"       "MPN"     
[15] "HNSC"     "KICH"     "KIRC"     "KIRP"     "LIHC"     "LUAD"     "LUSC"    
[22] "MESO"     "OV"       "PAAD"     "PCPG"     "PRAD"     "SARC"     "SKCM"    
[29] "STAD"     "TGCT"     "THCA"     "THYM"     "UCEC"     "UCS"      "UVM"   

#1.1 the proportion of indels in each cancer type and each MM gene
tmp1=aggregate(Mut[,3:7],by=list(Mut$Cancer,Mut$Gene),FUN=sum)
colnames(tmp1)<-colnames(Mut)[1:7]
mut_dat <- tmp1
mut_dat$Percent_of_Indel=round(mut_dat$Indel/(mut_dat$Indel+mut_dat$SNV),3) #number of mutations
mut_dat$Percent_of_Pts_with_Indel=round(mut_dat$Pts_with_Indel/mut_dat$Total_Pts,3) #number of samples
mut_dat_49 <- mut_dat[which(mut_dat[,2] %in% Nature_MM_merge),]

#1.2 the proportion of indels in each MM gene
mut_stat=aggregate(mut_dat[,c(3:7)],by=list(mut_dat$Gene),FUN=sum)
colnames(mut_stat)<-colnames(mut_dat)[2:7]
mut_stat$Percent_of_Indel=round(mut_stat$Indel/(mut_stat$Indel+mut_stat$SNV),3) #number of mutations
mut_stat$Percent_of_Pts_with_Indel=round(mut_stat$Pts_with_Indel/mut_stat$Total_Pts,3) #number of samples
mut_stat_49 <- mut_stat[which(mut_stat$Gene %in% Nature_MM_merge),]

#1.3 the proportion of indels in each cancer type
mut_stat_cancer=aggregate(mut_dat[,c(3:7)],by=list(mut_dat$Cancer),FUN=sum)
colnames(mut_stat_cancer)<-colnames(mut_dat)[c(1,3:7)]
mut_stat_cancer$Percent_of_Indel=round(mut_stat_cancer$Indel/(mut_stat_cancer$Indel+mut_stat_cancer$SNV),3)  #all genes
mut_stat_cancer_49=aggregate(mut_dat_49[,c(3:7)],by=list(mut_dat_49$Cancer),FUN=sum)
colnames(mut_stat_cancer_49)<-colnames(mut_dat_49)[c(1,3:7)]
mut_stat_cancer_49$Percent_of_Indel=round(mut_stat_cancer_49$Indel/(mut_stat_cancer_49$Indel+mut_stat_cancer_49$SNV),3)  #49 MM genes

#1.4 the proportion of indels in all cancer types
round(sum(mut_stat$Indel)/(sum(mut_stat$Indel)+sum(mut_stat$SNV)),3)   
[1] 0.088

#1.5 plot
library(ggrepel)
library(ggplot2)

plot_1 <- mut_dat_49[which(mut_dat_49$Total_Pts >= 30),]  #number of samples with mutations ≥ 30
colnames(plot_1)[7] <- "Mut_Pts"

### 1C_percent_of_pts_with_indel
plot_1C <- plot_1[order(plot_1$Percent_of_Pts_with_Indel),]
plot_1C$change<-as.factor(ifelse(plot_1C$Percent_of_Pts_with_Indel > 0.5,'High','Low'))
plot_1C$id=1:nrow(plot_1)
plot_1C$label<-ifelse(plot_1C$Percent_of_Pts_with_Indel > 0.5, paste(plot_1C$Gene,"(",plot_1C$Cancer,")"),"")

p <- ggplot(data = plot_1C, aes(x = id, y = Percent_of_Pts_with_Indel, color = change)) +
  xlim(0,nrow(plot_1C)+10) +
  geom_point(alpha=0.8, size = plot_1C$Percent_of_Pts_with_Indel*2) +
  xlab("Gene ranked by fraction")+ylab("Fraction of patients with indel")+ 
  theme_bw() +theme(panel.grid = element_blank(),
                    panel.border=element_blank(),
                    plot.title = element_text(hjust = 0))+
  guides(color = guide_legend(title="")) +
  theme(legend.position = c(0.35,0.8),
        legend.text = element_text(size=8),
        legend.key.size =unit(0.35,'cm') ) +
  theme(axis.title = element_text(size=10,colour = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=10,colour = "black"),
        axis.ticks.y = element_line(size=0.8),
        axis.line = element_line(colour = "black",size = 0.8))
p1<-p+geom_text_repel(data=plot_1C,aes(x = id, y = Percent_of_Pts_with_Indel,label=label),colour = "black",cex=2)
print(p1)
pdf("Fig1C_percent_of_pts_with_indel.pdf",height = 4,width = 4)
print(p1)
dev.off()


# 1B_1_percent of pts with indels
cancer_sample_number <- read.table("cancer_sample_number.txt", header=T, sep="")
plot_1B <- full_join(plot_1,cancer_sample_number,by = colnames(plot_1)[1])[1:nrow(plot_1),]
plot_1B$Percent_of_Pts_with_Mut <- round(plot_1B$Mut_Pts/plot_1B$Total_Pts,3)

mut_stat_cancer_49 <- mut_stat_cancer_49[order(mut_stat_cancer_49$Percent_of_Indel),]
mut_stat_49 <- mut_stat_49[order(mut_stat_49$Percent_of_Indel),]

p3 <- plot_1B %>% 
  ggplot()+geom_point( aes(x=factor(Cancer,levels = mut_stat_cancer_49$Cancer), 
                           y=factor(Gene,levels =mut_stat_49$Gene),   
                           size=Percent_of_Pts_with_Mut,
                           fill=Percent_of_Pts_with_Indel),shape=21)+
  theme_bw()+
  scale_fill_gradientn(values = seq(0,1,0.25),colours = c('#3d84a8','#a6e3e9','white','#ffb6b9','#e23e57'))+
  scale_radius(range=c(3,8), name="SIZE")+
  labs(x = "",y="",title="")+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0),               
        legend.position="top",
        legend.box = "vertical",
        legend.direction = "horizontal")+
  theme(axis.ticks = element_line(size=0.8),
        axis.text = element_text(colour="black",size=10),
        axis.text.x = element_text(angle=90,size=10,vjust = 0.5,hjust=1))+
  guides(fill = guide_colorbar(title="Patients with Indel"),
         size=guide_legend(title="Patients with mutations")) 
print(p3)
pdf("Fig1B_1_percent of pts with indels.pdf",width=7,height = 10)
print(p3)
dev.off()

# 1B_2_percent of indels in each cancer type
pdf("Fig1B_2_percent of indels in each cancer type.pdf",width=7,height=4)     
ggplot(mut_stat_cancer_49, aes(x=factor(Cancer,levels =unique(Cancer)),  
               y = Percent_of_Indel)) +
  geom_bar(fill='#3d84a8',stat="identity",width=0.7) +
  labs(x=" ", y="Fraction of indel") + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
dev.off()

# 1B_3_percent of indels in each MM gene
pdf("Fig1B_3_percent of indels in each MM gene.pdf",width=7,height=4)     
ggplot(mut_stat_49, aes(x=factor(Gene,levels =unique(Gene)),  
                               y = Percent_of_Indel)) +
  geom_bar(fill='#3d84a8',stat="identity",width=0.7) +
  labs(x=" ", y="Fraction of indel") + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
dev.off()


