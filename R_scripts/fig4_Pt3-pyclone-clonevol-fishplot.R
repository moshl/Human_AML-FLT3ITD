
library(data.table)
library(supraHex)
##define funtion for reshape the loci
Loci_tsv_To_clonevol <- function(dt){
  dc <- dcast(dt, formula = mutation_id + cluster_id ~ sample_id, value.var = "variant_allele_frequency")
  dc<-as.data.frame(dc) # it's important to translation
  #rownames(dc)<-dc[,1]
  #dt_out <- dc[,-1]
  #return(as.matrix(dt_out))
  return(dc)
}


###for citup input
setwd("D:\\Project\\project_yangjy_NR\\FLT3-Project\\results\\target/Ampliseq/Results/")

#nr5
dt <- fread("clonal/1S1Pyclone/excnvNR5/tables/loci.tsv")
clusters <- Loci_tsv_To_clonevol(dt)
head(clusters)

T1<-read.table("clonal/1S1Pyclone/excnvNR5/nr5_T1_excnv.pycloneinput.txt",sep="\t",header = T)
T2<-read.table("clonal/1S1Pyclone/excnvNR5/nr5_T2_excnv.pycloneinput.txt",sep="\t",header = T)
T3<-read.table("clonal/1S1Pyclone/excnvNR5/nr5_T3_excnv.pycloneinput.txt",sep="\t",header = T)

q<-match(clusters$mutation_id,T1$mutation_id)

T1<-T1[q,]
T2<-T2[q,]
T3<-T3[q,]

#prepare inputfle for clonevol
data<-cbind(clusters[,2:dim(clusters)[2]],T1[,c(2:3)],T2[,c(2:3)],T3[,c(2:3)])
rownames(data)<-clusters$mutation_id
colnames(data)<-c("cluster","D0.vaf","R1.vaf","R2.vaf",
                  "D0.ref.count","D0.var.count","R1.ref.count","R1.var.count",
                  "R2.ref.count","R2.var.count")
head(data)
data[is.na(data)]<- -1
dat<-na.omit(data)
dat<-as.data.frame(dat)
dat<-dat[dat$cluster>=0,]

##first-final
data<-cbind(clusters[,c(2,3,5)],T1[,c(2:3)],T3[,c(2:3)])
rownames(data)<-clusters$mutation_id
colnames(data)<-c("cluster","D0.vaf","R2.vaf",
                  "D0.ref.count","D0.var.count",
                  "R2.ref.count","R2.var.count")
head(data)
data[is.na(data)]<- -1
dat<-na.omit(data)
dat<-as.data.frame(dat)
dat<-dat[dat$cluster>=0,]
#dat$cluster<-dat$cluster+1
##########
vaf.col.names <- grep(".vaf", colnames(dat), value=TRUE)
head(dat)
head(dat[which(dat$cluster==1),])

library(clonevol)
x <- infer.clonal.models(variants=dat,
                         cluster.col.name="cluster",
                         vaf.col.names=vaf.col.names,
                         model="monoclonal",
                         subclonal.test="bootstrap",
                         subclonal.test.model="non-parametric",
                         num.boots=1000,
                         founding.cluster=1,
                         cluster.center="mean",
                         #ignore.clusters=c(5,6),
                         #clone.colors=clone.colors,
                         min.cluster.vaf=0.00001,
                         # min probability that CCF(clone) is non-negative
                         sum.p = 0.9,
                         # alpha level in confidence interval estimate for CCF(clone)
                         alpha = 0.01)


plot.clonal.models(x$models,
                   matched=x$matched,
                   variants=dat,
                   clone.shape="bell",
                   box.plot=TRUE,
                   out.format="pdf",
                   overwrite.output=TRUE,
                   scale.monoclonal.cell.frac=TRUE,
                   cell.frac.ci=TRUE,
                   tree.node.shape="circle",
                   tree.node.size=40,
                   tree.node.text.size=0.65,
                   width=11, height=5,
                   out.dir="clonal/NR5_clonevol")



################### fishplot
library(fishplot)
##NR5
timepoints=c(0,44,233,294) ##44ÊÇ¼ÓµÄÐéÄâµã
parents = c(0,1,1,3,4)
frac.table = matrix(
  c(89.01,54.23,33.93,0.00,0.00,
    6.00,0.00,5.00,2.50,1.00,
    62.57,0.00,54.04,26.00,9.10,
    99.00,0.00,90.89,81.23,40.84),
  ncol=length(timepoints))

fish = createFishObject(frac.table,parents,timepoints=timepoints)
fish = layoutClones(fish) #FD8EE5
#fish = setCol(fish,col= c('#d0ced0','gray','#ffa5be','#ff9095','#F8766D',
#                          '#FFD94B','orange','#2CD0AB')) #»Ò£¬»Ò2£¬·Ûºì£¬Éî·Û,ºì£¬ÂÌ£¬×ÏÂÞÀ¼£¬×Ï£¬»ÆÉ«, '#8d4891', '#f8e356'

fish = setCol(fish,col= c('#999999','#3333FF','#FF3333','#FF6666','#FF9999')) ###ffa5be #9933FF #FF6666 #9933FF

pdf("clonal/2S2Clonevol/pt3.fishplot.pdf",width = 8,height = 4)
fishPlot(fish, shape="spline", vlab=c("Day1","Day44","Day233","Day294"), vlines=timepoints, 
         col.vline ='white' ,col.border= "darkgray",border = 0.1,
         title="Patient3", cex.title=1.5,cex.vlab=1,ramp.angle=1, pad.left=0.3)

dev.off()










