library(ggplot2)

setwd("D:/ITD/Extended_2CD") 

info<-read.table("mITDvssITD_up_top15.txt",head=TRUE,sep="\t")
data<-data.frame(info)
factor=data[,"NES"]

pdf("Extended_2C_mITDvssITD_up_top15.pdf",width=8,height=6)

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
  scale_radius(range=c(3,8),                
    name="SIZE")+                             
  guides(   
    color = guide_colorbar(order = 1),        
    size = guide_legend(order = 2))+
  theme_bw()
dev.off()


info<-read.table("mITDvssITD_down_top15.txt",head=TRUE,sep="\t")
data<-data.frame(info)

pdf("Extended_2D_mITDvssITD_down_top15.pdf",width=8,height=6)
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