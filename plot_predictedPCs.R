args <- commandArgs(TRUE)
inferanc_file <- args[1] #the predicted PCs and inferred ancestry estimates

library(ggplot2)
library(data.table)

file_stem <- gsub(".out", "", inferanc_file)

#read in PCs
data <- fread(inferanc_file, data.table = F)
names(data) <- c("id","Pop_label", "#SNPs", "pred_PC1", "pred_PC2", "pred_PC3", "pred_PC4", "%SAS", "%EAS", "%AMR", "%AFR", "%EUR")
data$FID <- sapply(strsplit(data$id, ":"), "[", 1)
data$IID <- sapply(strsplit(data$id, ":"), "[", 2)
data[c("id", "Pop_label", "#SNPs")] <- NULL

# #check for the existence of a hash file and if so read it in to convert the ids
# geno_file_path <- gsub('(.*)/\\w+', '\\1', file_stem)
# if(length(list.files(path = geno_file_path, pattern = "_dummy_famids.txt"))>0){
#   id_hash <- list.files(path = geno_file_path, pattern = "_dummy_famids.txt")
#   id_hash <- fread(paste0(geno_file_path, "/", id_hash), data.table = F)
#   names(id_hash) <- c("FID", "IID", "FID_hash", "IID_hash")
#   names(data)[names(data)== "FID"] <- "FID_hash"
#   names(data)[names(data)== "IID"] <- "IID_hash"
#   data <- merge(data, id_hash, by = c("FID_hash", "IID_hash"), all.x = T)
#   
#   #write out file to update ids
#   update_ids_filename <- paste0(file_stem, "_update_ids.txt")
#   update_ids_file <- id_hash[,c("FID_hash", "IID_hash", "FID", "IID")]
#   write.table(update_ids_file, update_ids_filename, row.names = F, col.names = F, sep = " ", quote = F)
#   print(paste("Hash list of FID/IIDs found! Please be sure to update your fam file with the correct IDs using", update_ids_filename, "before subsetting for race or moving to imputation!"))
# }

# assign ancestries based on inferences
# European: >= 80% EUR
# AA: >= 80% combination of AFR & EUR and < 3& Asian and AMR
# Latino/Hispanic: >= 85% EUR & AMR, <10% AFR and <3% Asian
data$inferred_ancestry <- "Unknown"
data$inferred_ancestry[data$`%EUR`>0.7] <- "EUR"
data$inferred_ancestry[data$`%EUR`<0.7 & data$`%EUR` + data$`%AFR`>0.8] <- "AA"
data$inferred_ancestry[data$`%EUR`>0.8 & data$`%AMR`>0.8 & data$`%AFR`<0.10] <- "Latino/Hispanic"

# get which clusters are present in this set
ancestries <- unique(data$inferred_ancestry)

#create plots
pdf(paste0(file_stem, ".pdf"))

# Plot in everyone
a <- ggplot(data = data, aes(x=pred_PC1, y=pred_PC2, color=inferred_ancestry)) + geom_point() 
print(a)
b <- ggplot(data = data, aes(x=pred_PC2, y=pred_PC3, color=inferred_ancestry)) + geom_point() 
print(b)
c <- ggplot(data = data, aes(x=pred_PC3, y=pred_PC4, color=inferred_ancestry)) + geom_point() 
print(c)

for(i in 1:length(ancestries)){
  
  # subset to the current group
  tmp_data <- data[data$inferred_ancestry==ancestries[i],]

  # generate plots, still coloring on ancestral group so that it's clear which is being presented
  a <- ggplot(data = tmp_data, aes(x=pred_PC1, y=pred_PC2, color=inferred_ancestry)) + geom_point()
  print(a)
  b <- ggplot(data = tmp_data, aes(x=pred_PC2, y=pred_PC3, color=inferred_ancestry)) + geom_point() 
  print(b)
  c <- ggplot(data = tmp_data, aes(x=pred_PC3, y=pred_PC4, color=inferred_ancestry)) + geom_point() 
  print(c)
  
}

dev.off()
