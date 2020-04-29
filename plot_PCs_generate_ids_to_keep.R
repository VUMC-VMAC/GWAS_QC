args <- commandArgs(TRUE)
pcs_file <- args[1]
race_file <- args[2]
race_1000G_file <- args[3] #defines race file for 1000G; set to NA if 1000G were not included in these PCs
write_excl_file <- args[4] #defines whether or not to write out a file for default exclusion decisions


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
data$set <- "current"

if(race_1000G_file != "NA"){
	#read in race for 1000G
	data_1000G <- fread(race_1000G_file, header = F)
	names(data_1000G) <- c("FID", "IID", "race")
	data_1000G$set <- "1000G"
	data <- rbind(data[,c("FID", "IID", "race", "set")], data_1000G[,c("FID", "IID", "race", "set")])
}

#merge race and PCs
data <- merge(pcs, data, by = c("FID", "IID"))

#get NHW subset
data_nhw <- data[data$race %in% c("EUR", "White"),]

#print out how many rows
print(paste(nrow(data), "people with PCs, ",  nrow(data_nhw), "NHW"))

pdf(paste0(pc_file_stem, ".pdf"))

#color just to see if they cluster by 1000G race
ggplot(data = data, aes(x=PC1, y=PC2, color=race)) + geom_point()  +
  geom_vline(xintercept = c((mean(data_nhw$PC1)-5*sd(data_nhw$PC1)), (mean(data_nhw$PC1)+5*sd(data_nhw$PC1)))) +
  geom_hline(yintercept = c((mean(data_nhw$PC2)-5*sd(data_nhw$PC2)), (mean(data_nhw$PC2)+5*sd(data_nhw$PC2))))

ggplot(data = data, aes(x=PC2, y=PC3, color=race)) + geom_point()  +
  geom_vline(xintercept = c((mean(data_nhw$PC2)-5*sd(data_nhw$PC2)), (mean(data_nhw$PC2)+5*sd(data_nhw$PC2)))) +
  geom_hline(yintercept = c((mean(data_nhw$PC3)-5*sd(data_nhw$PC3)), (mean(data_nhw$PC3)+5*sd(data_nhw$PC3))))

ggplot(data = data, aes(x=PC3, y=PC4, color=race)) + geom_point()  +
  geom_vline(xintercept = c((mean(data_nhw$PC3)-5*sd(data_nhw$PC3)), (mean(data_nhw$PC3)+5*sd(data_nhw$PC3)))) +
  geom_hline(yintercept = c((mean(data_nhw$PC4)-5*sd(data_nhw$PC4)), (mean(data_nhw$PC4)+5*sd(data_nhw$PC4))))


#subset to just the current dataset
ggplot(data = data[data$set=="current",], aes(x=PC1, y=PC2, color = race)) + geom_point()  +
  geom_vline(xintercept = c((mean(data_nhw$PC1)-5*sd(data_nhw$PC1)), (mean(data_nhw$PC1)+5*sd(data_nhw$PC1)))) +
  geom_hline(yintercept = c((mean(data_nhw$PC2)-5*sd(data_nhw$PC2)), (mean(data_nhw$PC2)+5*sd(data_nhw$PC2))))

ggplot(data = data[data$set=="current",], aes(x=PC2, y=PC3, color = race)) + geom_point()  +
  geom_vline(xintercept = c((mean(data_nhw$PC2)-5*sd(data_nhw$PC2)), (mean(data_nhw$PC2)+5*sd(data_nhw$PC2)))) +
  geom_hline(yintercept = c((mean(data_nhw$PC3)-5*sd(data_nhw$PC3)), (mean(data_nhw$PC3)+5*sd(data_nhw$PC3))))

ggplot(data = data[data$set=="current",], aes(x=PC3, y=PC4, color = race)) + geom_point()  +
  geom_vline(xintercept = c((mean(data_nhw$PC3)-5*sd(data_nhw$PC3)), (mean(data_nhw$PC3)+5*sd(data_nhw$PC3)))) +
  geom_hline(yintercept = c((mean(data_nhw$PC4)-5*sd(data_nhw$PC4)), (mean(data_nhw$PC4)+5*sd(data_nhw$PC4))))


#restrict to non-hispanic whites 
ggplot(data = data[data$set=="current" & data$race == "White",], aes(x=PC1, y=PC2, color = race)) + geom_point()  +
  geom_vline(xintercept = c((mean(data_nhw$PC1)-5*sd(data_nhw$PC1)), (mean(data_nhw$PC1)+5*sd(data_nhw$PC1)))) +
  geom_hline(yintercept = c((mean(data_nhw$PC2)-5*sd(data_nhw$PC2)), (mean(data_nhw$PC2)+5*sd(data_nhw$PC2))))

ggplot(data = data[data$set=="current" & data$race == "White",], aes(x=PC2, y=PC3, color = race)) + geom_point()  +
  geom_vline(xintercept = c((mean(data_nhw$PC2)-5*sd(data_nhw$PC2)), (mean(data_nhw$PC2)+5*sd(data_nhw$PC2)))) +
  geom_hline(yintercept = c((mean(data_nhw$PC3)-5*sd(data_nhw$PC3)), (mean(data_nhw$PC3)+5*sd(data_nhw$PC3))))

ggplot(data = data[data$set=="current" & data$race == "White",], aes(x=PC3, y=PC4, color = race)) + geom_point()  +
  geom_vline(xintercept = c((mean(data_nhw$PC3)-5*sd(data_nhw$PC3)), (mean(data_nhw$PC3)+5*sd(data_nhw$PC3)))) +
  geom_hline(yintercept = c((mean(data_nhw$PC4)-5*sd(data_nhw$PC4)), (mean(data_nhw$PC4)+5*sd(data_nhw$PC4))))

dev.off()

if(write_excl_file == "yes"){
   #get NHW and non-outliers (according to the current sample mean)
   ids_to_keep <- data_nhw[data_nhw$set=="current" &
                      data_nhw$PC1 > mean(data_nhw$PC1)-5*sd(data_nhw$PC1) & data_nhw$PC1 < mean(data_nhw$PC1)+5*sd(data_nhw$PC1) &
                      data_nhw$PC2 > mean(data_nhw$PC2)-5*sd(data_nhw$PC2) & data_nhw$PC2 < mean(data_nhw$PC2)+5*sd(data_nhw$PC2) &
                      data_nhw$PC3 > mean(data_nhw$PC3)-5*sd(data_nhw$PC3) & data_nhw$PC3 < mean(data_nhw$PC3)+5*sd(data_nhw$PC3),
                    c("FID", "IID")]
   write.table(ids_to_keep, paste0(pc_file_stem, "_nooutliers.txt"), col.names = F, row.names = F, quote = F)
}
