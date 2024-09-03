library("ggplot2")
library("ggpubr")
library("patchwork")

# Extended_1D: diff of VAF(NGS) in sITD and mITD
setwd("D:/ITD/Extended_1D")
data <- read.csv("ITD_VAF_NGS_138_85.csv")

compare_means(VAF ~ ITD_type, data=data, method = "wilcox.test", paired = FALSE)
# A tibble: 1 × 8
  .y.   group1 group2          p    p.adj p.format p.signif method  
  <chr> <chr>  <chr>       <dbl>    <dbl> <chr>    <chr>    <chr>   
1 VAF   Single Multiple 3.29e-13 3.30e-13 3.3e-13  ****     Wilcoxon

pdf("Extended_1D_ITD_VAF_NGS_density.pdf",width = 5,height = 4)
p <- ggplot(data, aes(VAF)) +
  geom_histogram(aes(x=VAF, y=after_stat(density), fill = ITD_type), bins=30, position ='dodge', alpha=0.5) +
  geom_density(aes(fill = ITD_type, color = ITD_type), alpha = 0.3) +
  labs(x = "VAF") +
  theme_classic()
print(p)
dev.off()


#Extended_1E: diff of AA_Start, AA_End, Length, Domain in sITD and mITD
setwd("D:/ITD/Extended_1E")
info<-read.csv("ITD_input.csv", encoding = "UTF-8")
data <- info

compare_means(c(Protein_dupstart, Protein_dupend, Length) ~ ITD_type, data=data, method = "wilcox.test", paired = FALSE)
.y.              group1 group2       p p.adj p.format p.signif method  
<fct>            <chr>  <chr>    <dbl> <dbl> <chr>    <chr>    <chr>   
1 Protein_dupstart Single Multiple 0.527  0.53 0.53     ns       Wilcoxon
2 Protein_dupend   Single Multiple 0.743  0.74 0.74     ns       Wilcoxon
3 Length           Single Multiple 0.388  0.39 0.39     ns       Wilcoxon

p1 <- ggplot(data, aes(Protein_dupstart)) +
  geom_histogram(aes(x=Protein_dupstart, y=..density.., fill = ITD_type), bins=30, position ='dodge', alpha=0.5) +
  geom_density(aes(fill = ITD_type, color = ITD_type), alpha = 0.3) +
  labs(x = "AA_Start") +
  theme_classic() 

p2 <- ggplot(data, aes(Protein_dupend)) +
  geom_histogram(aes(x=Protein_dupend, y=..density.., fill = ITD_type), bins=30, position ='dodge', alpha=0.5) +
  geom_density(aes(fill = ITD_type, color = ITD_type), alpha = 0.3) +
  labs(x = "AA_End") +
  theme_classic() 

p3 <- ggplot(data, aes(Length)) +
  geom_histogram(aes(x=Length, y=..density.., fill = ITD_type), bins=30, position ='dodge', alpha=0.5) +
  geom_density(aes(fill = ITD_type, color = ITD_type), alpha = 0.3) +
  labs(x = "Length(bp)") +
  theme_classic() 

domain <- table(data$ITD_type, data$Insertion_domain)
domain <- as.data.frame(domain)
colnames(domain) <- c("Group","Domain","count")
domain <- domain[-c(1,2),]

#F8766D,#00BFC4
Custom.color <- c("#ECA680",  "#54BEAA", "#E3716E", "#F7DF87", "#BDB5E1", "#EFC0D2")
p4 <- domain %>%
  mutate(Group = factor(Group, levels=c("Single", "Multiple")))%>%  
  ggplot(aes(x = Group, y = count, fill=Domain)) +
  geom_bar(position="fill",stat="identity", width=0.5)+ 
  scale_fill_manual(values = Custom.color)+   
  scale_color_manual(values = Custom.color)+  
  ylab("Frequence") +   
  theme_classic()

p <- (p1 | p2)/(p3 | p4)   

pdf("Extended_1E_ITD_diversity.pdf",width = 10,height = 8)
print(p)
dev.off()
