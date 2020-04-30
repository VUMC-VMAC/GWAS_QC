args <- commandArgs(TRUE)
pcs_file <- args[1]
race_file <- args[2] #defines race file for main dataset; set to NA if plots should not be colored on race
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

if(race_file != "NA"){

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
	
		#merge race and PCs
		data <- merge(pcs, data, by = c("FID", "IID"))

		#get NHW subset
		data_nhw <- data[data$race %in% c("EUR", "White"),]
		
		#set outlier thresholds based on NHW in 1000G and current dataset
		PC1_thresh <- c((mean(data_nhw$PC1)-5*sd(data_nhw$PC1)), (mean(data_nhw$PC1)+5*sd(data_nhw$PC1)))
		PC2_thresh <- c((mean(data_nhw$PC2)-5*sd(data_nhw$PC2)), (mean(data_nhw$PC2)+5*sd(data_nhw$PC2)))
		PC3_thresh <- c((mean(data_nhw$PC3)-5*sd(data_nhw$PC3)), (mean(data_nhw$PC3)+5*sd(data_nhw$PC3)))
		PC4_thresh <- c((mean(data_nhw$PC4)-5*sd(data_nhw$PC4)), (mean(data_nhw$PC4)+5*sd(data_nhw$PC4)))

		#remove NHW dataframe since it's not needed anymore
		rm(data_nhw)
	} else {
		#merge race and PCs
		data <- merge(pcs, data, by = c("FID", "IID"))

		#set outlier thresholds based on whole sample
		PC1_thresh <- c((mean(data$PC1)-5*sd(data$PC1)), (mean(data$PC1)+5*sd(data$PC1)))
		PC2_thresh <- c((mean(data$PC2)-5*sd(data$PC2)), (mean(data$PC2)+5*sd(data$PC2)))
		PC3_thresh <- c((mean(data$PC3)-5*sd(data$PC3)), (mean(data$PC3)+5*sd(data$PC3)))
		PC4_thresh <- c((mean(data$PC4)-5*sd(data$PC4)), (mean(data$PC4)+5*sd(data$PC4)))
	}
} else {
	#if there's no race, just use PCs
	data <- pcs
	#set outlier thresholds based on whole sample
	PC1_thresh <- c((mean(data$PC1)-5*sd(data$PC1)), (mean(data$PC1)+5*sd(data$PC1)))
	PC2_thresh <- c((mean(data$PC2)-5*sd(data$PC2)), (mean(data$PC2)+5*sd(data$PC2)))
	PC3_thresh <- c((mean(data$PC3)-5*sd(data$PC3)), (mean(data$PC3)+5*sd(data$PC3)))
	PC4_thresh <- c((mean(data$PC4)-5*sd(data$PC4)), (mean(data$PC4)+5*sd(data$PC4)))
}	

#create plots
pdf(paste0(pc_file_stem, ".pdf"))

if(race_1000G_file != NA & race_file != NA){

	#color just to see if they cluster by 1000G race
	ggplot(data = data, aes(x=PC1, y=PC2, color=race)) + geom_point()  +
	  geom_vline(xintercept = PC1_thresh) + geom_hline(yintercept = PC2_thresh)
	ggplot(data = data, aes(x=PC2, y=PC3, color=race)) + geom_point()  +
	  geom_vline(xintercept = PC2_thresh) + geom_hline(yintercept = PC3_thresh)
	ggplot(data = data, aes(x=PC3, y=PC4, color=race)) + geom_point()  +
	  geom_vline(xintercept = PC3_thresh) + geom_hline(yintercept = PC4_thresh)

	#remove 1000G samples for the rest of the plots
	data <- data[data$set == "current",]
}
if(race_file != "NA"){
	#regardless of whether 1000G data is here or not, if race for this dataset is present, plot it with everyone
	ggplot(data = data, aes(x=PC1, y=PC2, color=race)) + geom_point()  +
	  geom_vline(xintercept = PC1_thresh) + geom_hline(yintercept = PC2_thresh)
	ggplot(data = data, aes(x=PC2, y=PC3, color=race)) + geom_point()  +
	  geom_vline(xintercept = PC2_thresh) + geom_hline(yintercept = PC3_thresh)
	ggplot(data = data, aes(x=PC3, y=PC4, color=race)) + geom_point()  +
	  geom_vline(xintercept = PC3_thresh) + geom_hline(yintercept = PC4_thresh)

	#plot in non-hispanic whites 
	data <- data[data$race == "White",]
	ggplot(data = data, aes(x=PC1, y=PC2, color=race)) + geom_point()  +
	  geom_vline(xintercept = PC1_thresh) + geom_hline(yintercept = PC2_thresh)
	ggplot(data = data, aes(x=PC2, y=PC3, color=race)) + geom_point()  +
	  geom_vline(xintercept = PC2_thresh) + geom_hline(yintercept = PC3_thresh)
	ggplot(data = data, aes(x=PC3, y=PC4, color=race)) + geom_point()  +
	  geom_vline(xintercept = PC3_thresh) + geom_hline(yintercept = PC4_thresh)
}
if(race_file == "NA"){
	#if no race information is available, plot everyone with lines based on whole sample
	ggplot(data = data, aes(x=PC1, y=PC2)) + geom_point()  +
	  geom_vline(xintercept = PC1_thresh) + geom_hline(yintercept = PC2_thresh)
	ggplot(data = data, aes(x=PC2, y=PC3)) + geom_point()  +
	  geom_vline(xintercept = PC2_thresh) + geom_hline(yintercept = PC3_thresh)
	ggplot(data = data, aes(x=PC3, y=PC4)) + geom_point()  +
	  geom_vline(xintercept = PC3_thresh) + geom_hline(yintercept = PC4_thresh)
}

dev.off()

if(write_excl_file == "yes"){
   #get NHW and non-outliers (according to the current sample mean)
   ids_to_keep <- data[data$PC1 > PC1_thresh[1] & data$PC1 < PC1_thresh[2] &
                      data$PC2 > PC2_thresh[1] & data$PC2 < PC2_thresh[2] &
                      data$PC3 > PC3_thresh[1] & data$PC3 < PC3_thresh[2],
                    c("FID", "IID")]
   write.table(ids_to_keep, paste0(pc_file_stem, "_nooutliers.txt"), col.names = F, row.names = F, quote = F)
}
