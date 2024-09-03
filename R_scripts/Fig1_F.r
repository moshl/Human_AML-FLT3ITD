setwd("F:/ITD/1F")

# 1F: Visualize the distribution of ITD sequences
info<-read.csv("ITD_input.csv",header=TRUE,sep=",")
data <- info
data<-data[-which(is.na(data$Protein_dupstart)=="TRUE"),]   
data<-data[-which(data$insertion_site_protein_as>630),]   

range(data$Protein_dupstart)
range(data$Protein_dupend)
data$Protein_dupstart<-as.numeric(data$Protein_dupstart)
data$Protein_dupend<-as.numeric(data$Protein_dupend)
data<-data[intersect(which(data$Protein_dupstart>0),which(data$Protein_dupend>0)),]

dat<-data
p<- order(dat[,"Protein_dupstart"])
dat_Sort<-dat[p,] 

dat_Sort$y<-c(1:dim(dat_Sort)[1])
dat_Sort$ITD<-paste("ITD",dat_Sort$y,sep = "")
dat_Sort$Ex_Protein<-as.character(dat_Sort$Ex_Protein)
s<-dat_Sort[,24:29]
s<-na.omit(s)  

dat_Sort<-dat_Sort[s$y,]
dat_Sort<-dat_Sort[-(nrow(dat_Sort)),]
dat_Sort<-dat_Sort[-(nrow(dat_Sort)),]

pdf("Fig1F_ITD_sequence_plot.pdf",width = 15,height = 20)
plot(x=1,y=1,xlim=c(570,630),ylim=c(1,dim(dat_Sort)[1]+6),type = "n",xlab = "",ylab = "",axes=F)

for (i in 1:nrow(dat_Sort)) {
  a=dat_Sort$Ex_Protein[i]
  b=dat_Sort$Protein_dupstart[i]
  c=dat_Sort$Protein_dupend[i]
  z=dat_Sort$y[i]
  lines(x=c(b:c),y=rep(z,c-b+1),type = "l",col="#00BFC4",lwd=4)
  print(i)
  lines(x=c((b-nchar(a)):b),y=rep(z,nchar(a)+1),type = "l",col="#F8766D",lwd=4)
  
}

{
  abline(v=572,col="lightgray",lwd=2)
  abline(v=578,col="lightgray",lwd=2)
  abline(v=592,col="lightgray",lwd=2)
  abline(v=603,col="lightgray",lwd=2)
  abline(v=609,col="lightgray",lwd=2)
  abline(v=615,col="lightgray",lwd=2)
  abline(v=623,col="lightgray",lwd=2)
  abline(v=630,col="lightgray",lwd=2)
  abline(h=dim(dat_Sort)[1]+6,col="gray",lwd=3)
  polygon(c(572,572,578,578),c(dim(dat_Sort)[1]+10,dim(dat_Sort)[1]+6.5,dim(dat_Sort)[1]+6.5,dim(dat_Sort)[1]+10),col='#F8766D',border='white',lwd=2)
  polygon(c(578,578,592,592),c(dim(dat_Sort)[1]+10,dim(dat_Sort)[1]+6.5,dim(dat_Sort)[1]+6.5,dim(dat_Sort)[1]+10),col='#F8766D',border='white',lwd=2)
  polygon(c(592,592,603,603),c(dim(dat_Sort)[1]+10,dim(dat_Sort)[1]+6.5,dim(dat_Sort)[1]+6.5,dim(dat_Sort)[1]+10),col='#F8766D',border='white',lwd=2)
  polygon(c(603,603,609,609),c(dim(dat_Sort)[1]+10,dim(dat_Sort)[1]+6.5,dim(dat_Sort)[1]+6.5,dim(dat_Sort)[1]+10),col='#F8766D',border='white',lwd=2)
  polygon(c(609,609,615,615),c(dim(dat_Sort)[1]+10,dim(dat_Sort)[1]+6.5,dim(dat_Sort)[1]+6.5,dim(dat_Sort)[1]+10),col='#F8766D',border='white',lwd=2)
  polygon(c(615,615,623,623),c(dim(dat_Sort)[1]+10,dim(dat_Sort)[1]+6.5,dim(dat_Sort)[1]+6.5,dim(dat_Sort)[1]+10),col='#F8766D',border='white',lwd=2)
  polygon(c(623,623,630,630),c(dim(dat_Sort)[1]+10,dim(dat_Sort)[1]+6.5,dim(dat_Sort)[1]+6.5,dim(dat_Sort)[1]+10),col='#F8766D',border='white',lwd=2)
  
  text(c(572,578,592,603,609,615,623,630),
       c(rep(dim(dat_Sort)[1]+5,8)),labels=c(572,578,592,603,609,615,623,630),col='black',cex=1)
  
  text(c(575,585,597.5,606,612,619,626.5),
       c(rep(dim(dat_Sort)[1]+8,7)),labels=c("JM-B","JM-S","JM-Z","HR","β1","NBL","β2"),col='black',cex=1)
  text(c(570),
       c(dim(dat_Sort)[1]+8),labels=c("FLT3"),col='black',cex=1)
  
  text(c(572:630),
       c(rep(dim(dat_Sort)[1]+6,630-572+1)),labels=c("Y","E","S","Q","L","Q","M","V","Q","V","T",
                                                     "G","S","S","D","N","E","Y","F","Y","V","D","F","R","E","Y","E","Y","D","L","K",
                                                     "W","E","F","P","R","E","N","L","E","F","G","K","V","L","G","S","G","A","F","G",
                                                     "K","V","M","N","A","T","A","Y"),col='black',font =2,cex=0.8)
}
dev.off()
