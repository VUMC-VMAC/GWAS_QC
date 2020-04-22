args <- commandArgs(TRUE)
pcs_file <- args[1]
race_file <- args[2]

library(ggplot2)
library(data.table)

pc_file_stem <- gsub(".pca.evec", "", pcs_file)

#read in PCs
#pcs <- read.table(pcs_file, skip = 1, header = F, stringsAsFactors = F)
pcs <- fread(pcs_file, data.table = F)
names(pcs) <- c("id", paste0("PC", 1:10), "pheno")
pcs$FID <- sapply(strsplit(pcs$id, ":"), "[", 1)
pcs$IID <- sapply(strsplit(pcs$id, ":"), "[", 2)
pcs[c("id", "pheno")] <- NULL

#check for the existence of a hash file and if so read it in to convert the ids
geno_file_path <- gsub('(.*)/\\w+', '\\1', pc_file_stem)
if(length(list.files(path = geno_file_path, pattern = "_dummy_famids.txt"))>0){
  id_hash <- list.files(path = geno_file_path, pattern = "_dummy_famids.txt")
#  id_hash <- read.table(paste0(geno_file_path, "/", id_hash), header = F, stringsAsFactors = F)
  id_hash <- fread(paste0(geno_file_path, "/", id_hash), data.table = F)
  names(id_hash) <- c("FID", "IID", "FID_hash", "IID_hash")
  names(pcs)[names(pcs)== "FID"] <- "FID_hash"
  names(pcs)[names(pcs)== "IID"] <- "IID_hash"
  pcs <- merge(pcs, id_hash, by = c("FID_hash", "IID_hash"))

  #write out file to update ids
  update_ids_filename <- paste0(pc_file_stem, "_update_ids.txt")
  update_ids_file <- id_hash[,c("FID_hash", "IID_hash", "FID", "IID")]
  write.table(update_ids_file, update_ids_filename, row.names = F, col.names = F, sep = " ", quote = F)
  print(paste("Hash list of FID/IIDs found! Please be sure to update your fam file with the correct IDs using", update_ids_filename, "before subsetting for race or moving to imputation!"))
}

#read in race
data <- fread(race_file, data.table = F, header = F)
names(data) <- c("FID", "IID", "race", "sex")
#merge
data <- merge(pcs, data, by = c("FID", "IID"))

pdf(paste0(pc_file_stem, ".pdf"))

#color just to see if they cluster by race
ggplot(data = data, aes(x=PC1, y=PC2, color=race)) + geom_point() +
  geom_vline(xintercept = c((mean(data$PC1)-5*sd(data$PC1)), (mean(data$PC1)+5*sd(data$PC1)))) +
  geom_hline(yintercept = c((mean(data$PC2)-5*sd(data$PC2)), (mean(data$PC2)+5*sd(data$PC2))))

ggplot(data = data, aes(x=PC2, y=PC3, color=race)) + geom_point() +
  geom_vline(xintercept = c((mean(data$PC2)-5*sd(data$PC2)), (mean(data$PC2)+5*sd(data$PC2)))) +
  geom_hline(yintercept = c((mean(data$PC3)-5*sd(data$PC3)), (mean(data$PC3)+5*sd(data$PC3))))

ggplot(data = data, aes(x=PC3, y=PC4, color=race)) + geom_point() +
  geom_vline(xintercept = c((mean(data$PC3)-5*sd(data$PC3)), (mean(data$PC3)+5*sd(data$PC3)))) +
  geom_hline(yintercept = c((mean(data$PC4)-5*sd(data$PC4)), (mean(data$PC4)+5*sd(data$PC4))))


#restrict to non-hispanic whites (or unknown)
ggplot(data = data[data$race=="White",], aes(x=PC1, y=PC2, color = race)) + geom_point()  +
  geom_vline(xintercept = c((mean(data$PC1)-5*sd(data$PC1)), (mean(data$PC1)+5*sd(data$PC1)))) +
  geom_hline(yintercept = c((mean(data$PC2)-5*sd(data$PC2)), (mean(data$PC2)+5*sd(data$PC2))))

ggplot(data = data[data$race=="White",], aes(x=PC2, y=PC3, color = race)) + geom_point()  +
  geom_vline(xintercept = c((mean(data$PC2)-5*sd(data$PC2)), (mean(data$PC2)+5*sd(data$PC2)))) +
  geom_hline(yintercept = c((mean(data$PC3)-5*sd(data$PC3)), (mean(data$PC3)+5*sd(data$PC3))))

ggplot(data = data[data$race=="White",], aes(x=PC3, y=PC4, color = race)) + geom_point()  +
  geom_vline(xintercept = c((mean(data$PC3)-5*sd(data$PC3)), (mean(data$PC3)+5*sd(data$PC3)))) +
  geom_hline(yintercept = c((mean(data$PC4)-5*sd(data$PC4)), (mean(data$PC4)+5*sd(data$PC4))))

dev.off()


#get NHW and non-outliers
ids_to_keep <- data[data$race == "White" &
                      data$PC1 > mean(data$PC1)-5*sd(data$PC1) & data$PC1 < mean(data$PC1)+5*sd(data$PC1) &
                      data$PC2 > mean(data$PC2)-5*sd(data$PC2) & data$PC2 < mean(data$PC2)+5*sd(data$PC2) &
                      data$PC3 > mean(data$PC3)-5*sd(data$PC3) & data$PC3 < mean(data$PC3)+5*sd(data$PC3),
                    c("FID", "IID")]
write.table(ids_to_keep, paste0(pc_file_stem, "_withoutoutliers.txt"), col.names = F, row.names = F, quote = F)
