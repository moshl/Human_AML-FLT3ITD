library(ggplot2)
library(reshape2)
library("ggsci")  

setwd('D:/ITD/Extended_3AB')

### Extended_3A: number of SNV
info<-read.table("WES_WGS_31_SNV_number.txt",head=TRUE,sep="\t")
df<-data.frame(info)

df<-melt(
  df,
  id.vars=c("Samples", "Group"), 
  variable.name = "SNV", 
  value.name = "number" 
)

pdf("Extended_3A_WES_WGS_31_exon_SNV_number.pdf",width=8,height=6)
ggplot(df, aes(x=factor(Samples,levels =unique(Samples)),  
               y = number, fill=factor(SNV,levels = unique(SNV)), ))+
  labs(x=" ", y="Number of SNV", fill="SNV") + 
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5)) +  
  theme(legend.key = element_blank(), panel.background = element_blank()) + 
  geom_bar(position="stack",stat="identity") + 
  facet_grid(.~factor(Group,levels =unique(Group)), scales = "free_x") +  # , ncol=3
  scale_fill_locuszoom() +
  theme_classic()
dev.off()


### Extended_3B: proportion of SNV Signature
info<-read.table("WES_WGS_31_SNV_Signature.txt",head=TRUE,sep="\t")
df<-data.frame(info)

df<-melt(
  df,
  id.vars=c("Samples", "Group"),  
  variable.name = "Signature",#
  value.name = "Activity" 
)

pdf("Extended_3B_WES_WGS_31_signature_exon_v3.2_proportion.pdf",width=8,height=6)
ggplot(df, aes(x=factor(Samples,levels =unique(Samples)),  
               y = Activity, fill=factor(Signature,levels = unique(Signature)), ))+
  labs(x=" ", y="Proportion of  each signatrure", fill="COSMIC_SBS96_Score" ) + 
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5)) +  
  theme(legend.key = element_blank(), panel.background = element_blank()) + 
  geom_bar(position="fill",stat="identity") + 
  facet_grid(.~factor(Group,levels =unique(Group)), scales = "free_x") +  # , ncol=3
  scale_fill_manual(values=c(
    "SBS1"="#f69c9f",
    "SBS5"="#fedcbd",
    "SBS4"="#9b95c9",
    "SBS15"="#76becc"))+
  theme_classic()

dev.off()

