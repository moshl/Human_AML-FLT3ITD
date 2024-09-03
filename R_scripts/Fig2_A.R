################################################################################
### All genes VAF distribution ###
################################################################################
# ================ #
# Load packages ####
# ================ #
# Package names
packages <- c("ggbreak","ggplot2", "dplyr", "cometExactTest", "stringr", "maditr", "reshape2", "data.table", "epitools", "corrplot", "plyr", "muhaz", "reshape", "survival", "survivalAnalysis", "survMisc", "survminer", "ggsci", "vegan", "ggrepel", "ggforce", "rstatix", "effsize", "psych", "maxstat", "RCurl", "ggpubr", "UpSetR", "cowplot", "readxl", "scales", "rlist", "ggplotify", "ggridges", "gridExtra", "BradleyTerry2")
# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
# Packages loading
invisible(lapply(packages, library, character.only = TRUE))
rm(list=ls())

###############################  load data  #######################################
setwd("D:/ITD/2A")
sample_list <- readxl::read_excel("Sample list.xlsx")   
s_ITD <- sample_list$ID[sample_list$ITD_type=="Single"]
m_ITD <- sample_list$ID[sample_list$ITD_type=="Multiple"]
all_sample <- c(m_ITD,s_ITD)

sum_mat2 <- read.csv("2A_final_data_sum.csv")   
sum_mat2$Gene[sum_mat2$ITD_type %in% "Multiple" & sum_mat2$Gene %in% "FLT3_ITD"] <- "FLT3_ITD(m)"   
sum_mat2$Gene[sum_mat2$ITD_type %in% "Single" & sum_mat2$Gene %in% "FLT3_ITD"] <- "FLT3_ITD(s)"    

i=1   
for(i in 1:length(m_ITD)) {
  temp_filter <- which(sum_mat2$ID == m_ITD[i] & sum_mat2$Gene == "FLT3_ITD(m)")
  temp_major <- temp_filter[which.max(sum_mat2$VAF[temp_filter])]
  temp_minor <- temp_filter[-which.max(sum_mat2$VAF[temp_filter])]
  sum_mat2$Gene[temp_major] <- "Major_ITD(m)"
  sum_mat2$Gene[temp_minor] <- "Minor_ITD(m)"
}


DNA_Methylation <- c("DNMT3A","IDH1","IDH2","TET2")
Chromatin_Cohesin <- c("ASXL1","STAG2","EZH2","BCOR")
Signaling <- c("PTPN11","KRAS","NRAS","FLT3_TKD","CBL","NOTCH1","MPL","KIT","JAK2","JAK3","CSF3R")
Splicing <- c("SRSF2","U2AF1","ZRSR2")
Tumor_Suppressors <- c("TP53","PHF6","WT1","CDKN2A")
Transcription_Factor <- c("CEBPA","GATA2","RUNX1","ETV6","SETBP1") 
NPM1 <- c("NPM1")
FLT3_ITD <- c("Major_ITD(m)","Minor_ITD(m)","FLT3_ITD(s)")
Control <- c("Control")
Type_list <- c("Control","FLT3_ITD","NPM1",
               "DNA_Methylation","Transcription_Factor",
               "Tumor_Suppressors","Signaling",
               "Chromatin_Cohesin","Splicing")
values = c("DNA_Methylation" = "#54beaa",
           "Chromatin_Cohesin" = "#eca680",
           "Signaling" = "#7ac7e2",
           "FLT3_ITD"="#f8766d",
           "Splicing" = "burlywood4",
           "Transcription_Factor" = "#b0b0ff",
           "NPM1" = "cyan3",
           "Tumor_Suppressors" = "#f7df87")

Gene_list <- unlist(sapply(Type_list,function(x) {
  rep(x,length(get(x)))
}))
names(Gene_list) <- unlist(mget(Type_list))
Color_list <- values[Gene_list]
names(Color_list) <- names(Gene_list)   
Color_list <- Color_list[-1]   

###############################  select genes  #######################################
sum_mat <- sum_mat2[!is.na(sum_mat2$VAF),]   
TOP_gene_list <- names(sort(table(sum_mat$Gene),decreasing = T))     
TOP_gene_list <- TOP_gene_list[TOP_gene_list %in% names(Gene_list)]  
sum_mat <- sum_mat[sum_mat$Gene %in% TOP_gene_list,]
sum_mat <- sum_mat[sum_mat$Gene %in% names(table(sum_mat$Gene)[table(sum_mat$Gene) >= 9]),]  

temp <- unique(sum_mat$Gene)
temp 
#[1] "TET2"         "NOTCH1"       "IDH2"         "IDH1"         "DNMT3A"       "CEBPA"        "WT1"         
#[8] "NRAS"         "RUNX1"        "FLT3_TKD"     "NPM1"         "FLT3_ITD(s)"  "Minor_ITD(m)" "Major_ITD(m)"


###############################  the order in any two genes#######################################
genes=as.data.frame(t(combn(as.vector(temp),2)))  
cutoff <- 5

results_list <- list()
n=1  
i=1
for(i in 1:nrow(genes)){
  # select mutations of interest
  gene_1 <- as.character(genes$V1[i])
  gene_2 <- as.character(genes$V2[i])
   
  sub_gene_1 <- subset(sum_mat, sum_mat$Gene == gene_1)  
  sub_gene_2 <- subset(sum_mat, sum_mat$Gene == gene_2)  
  
  # select patients with both mutations   
  gene_1_and_2 <- inner_join(sub_gene_1, sub_gene_2, by = "ID")  
  
  n_cases <- as.numeric(nrow(gene_1_and_2))
  if(n_cases==0) {
    temp_dat <- data.frame(matrix(NA, nrow = 1, ncol = 4))
    colnames(temp_dat) = c("Mutation1.Gene", "Mutation2.Gene", "Gene1.wins", "Gene2.wins")
    temp_dat[n,1] <- (gene_1) 
    temp_dat[n,2] <- (gene_2)
    temp_dat[n,3] <- 0
    temp_dat[n,4] <- 0
    results_list[[n]] <- temp_dat
    n=n+1     
    next
  }  
  
  # add point color columns for visualizing clonal/subclonal trends
  gene_1_and_2$vaf_ratio <- NA
  gene_1_and_2$vaf_ratio <- as.numeric((gene_1_and_2$VAF.x - gene_1_and_2$VAF.y))
  gene_1_and_2$vaf_ratio <- as.numeric(gene_1_and_2$vaf_ratio)
  gene_1_and_2 <- gene_1_and_2[complete.cases(gene_1_and_2$vaf_ratio), ]   # 筛选出包含缺失值的所有数据行，TRUE包含，FALSE不包含
  
  # define order
  if(as.numeric(nrow(gene_1_and_2))==0) {
    temp_dat <- data.frame(matrix(NA, nrow = 1, ncol = 4))
    colnames(temp_dat) = c("Mutation1.Gene", "Mutation2.Gene", "Gene1.wins", "Gene2.wins")
    temp_dat[n,1] <- (gene_1) 
    temp_dat[n,2] <- (gene_2)
    temp_dat[n,3] <- 0
    temp_dat[n,4] <- 0
    results_list[[n]] <- temp_dat
    n=n+1     
    next
  }  
  
  gene_1_and_2$Clonality <- 0
  
  if(nrow(gene_1_and_2) > 0){  
    for(k in 1:nrow(gene_1_and_2)){
      if(gene_1_and_2$vaf_ratio[k] <= cutoff & gene_1_and_2$vaf_ratio[k] >= -cutoff){
        gene_1_and_2$Clonality[k] <- 0
      }
      if(gene_1_and_2$vaf_ratio[k] > cutoff){
        gene_1_and_2$Clonality[k] <- 1
      }
      if(gene_1_and_2$vaf_ratio[k] < -cutoff){
        gene_1_and_2$Clonality[k] <- 2
      }
    }
  }
  gene_1_and_2$Clonality <- as.numeric(gene_1_and_2$Clonality)
  
  # fraction of cases where gene x occurs before gene y
  n_1_before_2 <- as.numeric(length(which(gene_1_and_2$Clonality == 1)))
  n_2_before_1 <- as.numeric(length(which(gene_1_and_2$Clonality == 2)))
  
  temp_dat <- data.frame(matrix(NA, nrow = 1, ncol = 4))
  colnames(temp_dat) = c("Mutation1.Gene", "Mutation2.Gene", "Gene1.wins", "Gene2.wins")
  
  temp_dat[n,1] <- (gene_1) 
  temp_dat[n,2] <- (gene_2)
  temp_dat[n,3] <- (n_1_before_2) 
  temp_dat[n,4] <- (n_2_before_1) 
  
  results_list[[n]] <- temp_dat
  n=n+1     
  
}

temp_final = na.omit(as.data.frame(do.call(rbind, results_list)))


###############################  execute BT model  #######################################
temp <- TOP_gene_list
temp <- unique(c("NPM1",temp))   
temp_final$Mutation1.Gene = factor(temp_final$Mutation1.Gene,levels=temp)  
temp_final$Mutation2.Gene = factor(temp_final$Mutation2.Gene,levels=temp) 

# now that we have the wins for each mutation counted, apply the BT model
BT_results = BTm(cbind(Gene1.wins, Gene2.wins), Mutation1.Gene, Mutation2.Gene, data = temp_final)
BT_MetaAML_mutation_ordering=as.data.frame(summary(BT_results)$coefficients)
BT_MetaAML_mutation_ordering$Gene = rownames(BT_MetaAML_mutation_ordering)
library(stringr)
BT_MetaAML_mutation_ordering$Gene = substring(BT_MetaAML_mutation_ordering$Gene, 3)  
# manually add ASXL1 because it gets removed in the BTm for some reason
BT_MetaAML_mutation_ordering = add_row(BT_MetaAML_mutation_ordering)
# add functional category to the mutations for visualization purposes
BT_MetaAML_mutation_ordering$mutation_category <- NA
# temp_BT <- BT_MetaAML_mutation_ordering[!BT_MetaAML_mutation_ordering$`Std. Error` > 10,]
temp_BT <- BT_MetaAML_mutation_ordering[-nrow(BT_MetaAML_mutation_ordering),]
temp_BT <- temp_BT[order(temp_BT$Estimate),]

temp_BT[1,2] <- 0.0
BTabilities(BT_results)     

### lable mutation_category
gene_order <- unlist(substring(rownames(temp_BT),3))
sum_mat$mutation_category <- NA
sum_mat$mutation_category <- sum_mat$Gene # Gene_list[sum_mat$Gene]
sum_mat$mutation_category <- factor(sum_mat$mutation_category, levels = gene_order)
# unique(sum_mat$ID[sum_mat$Gene=="NPM1"])
sum_mat <- sum_mat[!is.na(sum_mat$mutation_category),]    


###############################  plot  #######################################
temp_BT$ypos <- temp_BT$Estimate - temp_BT$`Std. Error` - 0.4
temp_BT$ypos <- temp_BT$Estimate
temp_BT$num <- table(sum_mat2$Gene)[gene_order]   
temp_BT$text <- paste("n = ",temp_BT$num,sep="")
temp_BT$Category <- Gene_list[temp_BT$Gene]
temp_BT$ypos[1] <- temp_BT$ypos[1] + 0.5
values <- values[unique(rev(temp_BT$Category))]

library(ggplot2)
library(dplyr)
library(forcats)

gene_new_sort <- temp_BT[sort(temp_BT$Estimate,decreasing = F,index.return = TRUE)$ix,"Gene"] 
gene_order_2 <- gene_new_sort

s_b <- temp_BT %>%
  mutate(Gene = fct_relevel(Gene, gene_order_2)) %>%
  ggplot(aes(x=Gene, y=Estimate,
             ymin=(Estimate-`Std. Error`),ymax=(Estimate+`Std. Error`))) +
  geom_pointrange(size = 0.75, aes(x = Gene, y = Estimate, color = Category, 
                                   ymin = (Estimate-`Std. Error`), ymax = (Estimate+`Std. Error`))) +        
  scale_y_break(c(-10, -24), scales=5) +
  scale_color_manual(name = "Mutation Category",
                     values = values,
                     breaks = names(values)) +
  scale_y_continuous(position = "right") +
  ylab("Point Estimate + 95% CI")+
  theme_cowplot() +
  theme(
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.line.y = element_blank()) +
  theme(legend.position = c(0.45, 0.8))   + 
  coord_flip() + 
  scale_y_reverse() +
  geom_text(aes(y=ypos,label=paste("n = ",num,sep=""),color="azrue4"),
            vjust=-1,size=4)

s_b


s_c <- sum_mat %>% 
  mutate(mutation_category = fct_relevel(mutation_category, gene_order_2)) %>%
  ggplot(aes(y = mutation_category, x = VAF, fill = mutation_category, height = ..density..)) +
  geom_vline(xintercept = 50, linewidth = 0.3, linetype = 2) +
  stat_density_ridges(
    alpha = 1,
    panel_scaling=T,
    scale = 1.1,
    quantile_lines = TRUE, 
    quantiles = 2) +
  ylab(label = NULL) +
  xlab(label = "VAF") +
  xlim(0,100) +
  # ylim(0,1) +
  theme_cowplot() +
  theme(legend.title = element_text()) +
  scale_fill_manual(name = "",values = Color_list) +
  theme(legend.position = "none")

s_c

# merge s_b and s_c
s_order <- ggarrange(s_c,s_b,
                     ncol = 2, nrow = 1, widths = c(1,1.2))

s_order

ggsave(paste("Fig2A_BT model.pdf",sep=""),s_order,width=9,height=7.5)  
