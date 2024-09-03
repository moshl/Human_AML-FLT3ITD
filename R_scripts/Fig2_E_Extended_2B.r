library(ggplot2)

setwd("F:/ITD/2E") 

info<-read.table("TCGA_ITDvsWT_up_top15.txt",head=TRUE,sep="\t")
data<-data.frame(info)
factor=data[,"NES"]

pdf("2E_ITDvsWT_up_top15.pdf",width=6.5,height=6)

ggplot(data,aes(x = NES, 
                y = reorder(NAME,NES), 
                size = SIZE,
                colour=FDR.q.val)) +
  geom_point(shape = 16) +                    
  labs(x = "Normalized Enrichment Score", y = "Pathway")+           
  scale_colour_continuous(                    #
    name="p.adjust",                       
    low="red",                              
    high="#FFBEB2")+
  scale_radius(                               
    range=c(3,8),                             
    name="SIZE")+                             
  guides(   
    color = guide_colorbar(order = 1),        
    size = guide_legend(order = 2)
  )+
  theme_bw()
dev.off()


info<-read.table("Extended_2B_TCGA_ITDvsWT_down_top15.txt",head=TRUE,sep="\t")
data<-data.frame(info)

pdf("ITDvsWT_down_top15.pdf",width=8,height=6)
ggplot(data,aes(x = NES, 
                y = reorder(NAME,NES), 
                size = SIZE,
                colour=FDR.q.val)) +
  geom_point(shape = 16) +                    
  labs(x = "Normalized Enrichment Score", y = "Pathway")+           
  scale_colour_continuous(                    
    name="p.adjust",                        
    low="red",                              
    high="#FFBEB2")+
  scale_radius(                               
    range=c(3,8),                             
    name="SIZE")+                             
  guides(   
    color = guide_colorbar(order = 1),        
    size = guide_legend(order = 2)
  )+
  theme_bw()
dev.off()



