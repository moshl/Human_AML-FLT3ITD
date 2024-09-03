library("ggplot2")
library("ggpubr")
library("patchwork")


# 1I: diff of VAF(PCR) in sITD and mITD
setwd("D:/ITD/1I")
data <- read.csv("ITD_VAF_PCR_130_81.csv")

compare_means(VAF ~ ITD_type, data=data, method = "wilcox.test", paired = FALSE)
# A tibble: 1 × 8
  .y.   group1 group2          p    p.adj p.format p.signif method  
  <chr> <chr>  <chr>       <dbl>    <dbl> <chr>    <chr>    <chr>   
1 VAF   Single Multiple 6.50e-16 6.5e-16 6.5e-16  ****     Wilcoxon

pdf("1I_ITD_VAF_PCR_density.pdf",width = 5,height = 4)
p <- ggplot(data, aes(VAF)) +
  geom_histogram(aes(x=VAF, y=after_stat(density), fill = ITD_type), bins=30, position ='dodge', alpha=0.5) +
  geom_density(aes(fill = ITD_type, color = ITD_type), alpha = 0.3) +
  labs(x = "VAF") +
  theme_classic()
print(p)
dev.off()

#the VAF of top peak in sITD
sITD_VAF <- data[which(data$ITD_type == "Single"), "VAF"]
y_peak <- which.max(density(sITD_VAF)$y) 
x_peak <- density(sITD_VAF)$x[y_peak]
> x_peak
[1] 39.92038

#the VAF of top peak in mITD
mITD_VAF <- data[which(data$ITD_type == "Multiple"), "VAF"]
plot(density(mITD_VAF))
density_mITD <- cbind(density(mITD_VAF)$x,density(mITD_VAF)$y)
write.csv(density_mITD, file="density_mITD_VAF_PCR.csv")
