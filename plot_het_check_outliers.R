args <- commandArgs(TRUE)
het_file_stem <- args[1]

library(data.table)
library(ggplot2)

data <- fread(paste0(het_file_stem, ".het"))

#define outliers
data$groups <- "<3SD"
data$groups[data$`F`<mean(data$`F`)-6*sd(data$`F`) | data$`F`>mean(data$`F`)+6*sd(data$`F`)] <- "6 SD"
data$groups[data$`F`<mean(data$`F`)-5*sd(data$`F`) | data$`F`>mean(data$`F`)+6*sd(data$`F`)] <- "5 SD"
data$groups[data$`F`<mean(data$`F`)-4*sd(data$`F`) | data$`F`>mean(data$`F`)+4*sd(data$`F`)] <- "4 SD"
data$groups[data$`F`<mean(data$`F`)-3*sd(data$`F`) | data$`F`>mean(data$`F`)+3*sd(data$`F`)] <- "3 SD"
levels(data$groups) <- c("<3 SD", "3 SD", "4 SD", "5 SD", "6 SD")

#plot heterozygosity
png(paste0(het_file_stem, ".png"), type = "cairo")
ggplot(data = data, aes(x = `O(HOM)`, y = `E(HOM)`, color = groups)) + 
  geom_point() + scale_color_manual(values = c("grey80", "grey70", "grey50", "grey40", "black")) +
  theme_bw() + labs(x = "Observed Homozygosity", y = "Expected Homozygosity", color = NULL)
dev.off()

#write out ids of 6 sd outliers if present
if(sum(data$groups == "6 SD")>0){
  
  outliers <- data[data$groups == "6 SD",c("FID", "IID")]
  print(paste(nrow(outliers), "heterozygosity outliers to remove. See plots for more details."))
  write.table(outliers, paste0(het_file_stem, "_outliers.txt"), row.names = F, col.names = F, sep = " ", quote = F)
  
} else {
  print("No heterozygosity outliers.")
}