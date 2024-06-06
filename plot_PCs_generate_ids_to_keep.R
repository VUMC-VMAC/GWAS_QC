args <- commandArgs(TRUE)
pcs_file <- args[1]
race_file <- args[2] #defines race file for main dataset; set to "none" if plots should not be colored on self-report race
race_1000G_file <- args[3] #defines race file for 1000G; set to "none" if 1000G were not included in these PCs
write_excl_file <- args[4] #defines whether or not to write out a file for default exclusion decisions
dataset_label <- args[5] #defines label for the race categories in the current set; set to "none" if unspecified

library(ggplot2)
library(data.table)

###################### Read in PCs ########################

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
  pcs <- merge(pcs, id_hash, by = c("FID_hash", "IID_hash"), all.x = T)
  
  if(race_1000G_file != "none"){
    #replace the 1000G ids if that file is present
    pcs$FID[is.na(pcs$FID)] <- pcs$FID_hash[is.na(pcs$FID)]
    pcs$IID[is.na(pcs$IID)] <- pcs$IID_hash[is.na(pcs$IID)]
  }
  
  #write out file to update ids
  update_ids_filename <- paste0(pc_file_stem, "_update_ids.txt")
  update_ids_file <- id_hash[,c("FID_hash", "IID_hash", "FID", "IID")]
  write.table(update_ids_file, update_ids_filename, row.names = F, col.names = F, sep = " ", quote = F)
  print(paste("Hash list of FID/IIDs found! Please be sure to update your fam file with the correct IDs using", update_ids_filename, "before subsetting for race or moving to imputation!"))
}

###################### Set up plotting df ########################

# set up df for plotting by reading in self-report/external race categories and generating thresholds

## if cohort race file is present, read in. 
## Otherwise, generate df with IDs and blank columns
if(race_file != "none"){
  
  #read in race
  data <- fread(race_file, data.table = F, header = F)
  names(data) <- c("FID", "IID", "race")
  data$set <- "current"
  
} 

## if 1000G data was used in generating the PCs, read in race file and combine with the cohort df
if(race_1000G_file != "none"){
  
  #read in race for 1000G
  data_1000G <- fread(race_1000G_file, header = F)
  names(data_1000G) <- c("FID", "IID", "race")
  data_1000G$set <- "1000G"
  
  #add 1000G to race values
  data_1000G$race <- paste("1000G", data_1000G$race, sep = " ")
  
}

## if both are present, combine, but if just 1000G is present, rename
if(race_file != "none" && race_1000G_file != "none"){
  data <- rbind(data[,c("FID", "IID", "race", "set")], data_1000G[,c("FID", "IID", "race", "set")])
} else if(race_1000G_file != "none"){
  data <- data_1000G
  rm(data_1000G)
}

## if either is present, merge with PCs
if(race_file != "none" || race_1000G_file != "none"){
  
  # merge with PCs, keeping all PC values
  data <- merge(pcs, data, by = c("FID", "IID"), all.x = T)
  
  # set the NA values if there are any (which would be the case if there were no self-report race file)
  data$set[is.na(data$set)] <- "current"
  data$race[is.na(data$race)] <- "not reported"
  
  #get NHW subset 
  # Note that this will be based on 1000G only if there are no NHW in the set or if the self-report is not present
  data_nhw <- data[data$race %in% c("1000G EUR", "White", "EUR"),]
  
  #set outlier thresholds based on NHW in 1000G and current dataset
  PC1_thresh <- c((mean(data_nhw$PC1)-5*sd(data_nhw$PC1)), (mean(data_nhw$PC1)+5*sd(data_nhw$PC1)))
  PC2_thresh <- c((mean(data_nhw$PC2)-5*sd(data_nhw$PC2)), (mean(data_nhw$PC2)+5*sd(data_nhw$PC2)))
  PC3_thresh <- c((mean(data_nhw$PC3)-5*sd(data_nhw$PC3)), (mean(data_nhw$PC3)+5*sd(data_nhw$PC3)))
  PC4_thresh <- c((mean(data_nhw$PC4)-5*sd(data_nhw$PC4)), (mean(data_nhw$PC4)+5*sd(data_nhw$PC4)))
  
  #remove NHW dataframe since it's not needed anymore
  rm(data_nhw)
} else {
  
  ## if there's no race, just use PCs
  data <- pcs
  #set outlier thresholds based on whole sample
  PC1_thresh <- c((mean(data$PC1)-5*sd(data$PC1)), (mean(data$PC1)+5*sd(data$PC1)))
  PC2_thresh <- c((mean(data$PC2)-5*sd(data$PC2)), (mean(data$PC2)+5*sd(data$PC2)))
  PC3_thresh <- c((mean(data$PC3)-5*sd(data$PC3)), (mean(data$PC3)+5*sd(data$PC3)))
  PC4_thresh <- c((mean(data$PC4)-5*sd(data$PC4)), (mean(data$PC4)+5*sd(data$PC4)))
}

#regardless of whether 1000G data is present or not,
#add dataset label to this dataset's race categories
#if a label was supplied
if(dataset_label != "none"){
  data$race[data$set == "current"] <- paste(dataset_label, data$race[data$set == "current"], sep = " ")
}

###################### Generate plots ########################

#create plots
pdf(paste0(pc_file_stem, ".pdf"))

if(race_1000G_file != "none"){
  
  #color just to see if they cluster by 1000G race
  a <- ggplot(data = data) +
    geom_point(data = data[data$set == "1000G",], aes(x=PC1, y=PC2, color=race))  +
    geom_point(data = data[data$set == "current",], aes(x=PC1, y=PC2, color=race))  +
    geom_vline(xintercept = PC1_thresh) + geom_hline(yintercept = PC2_thresh)
  print(a)
  b <- ggplot(data = data) +
    geom_point(data = data[data$set == "1000G",], aes(x=PC2, y=PC3, color=race))  +
    geom_point(data = data[data$set == "current",], aes(x=PC2, y=PC3, color=race))  +
    geom_vline(xintercept = PC2_thresh) + geom_hline(yintercept = PC3_thresh)
  print(b)
  c <- ggplot(data = data) +
    geom_point(data = data[data$set == "1000G",], aes(x=PC3, y=PC4, color=race))  +
    geom_point(data = data[data$set == "current",], aes(x=PC3, y=PC4, color=race))  +
    geom_vline(xintercept = PC3_thresh) + geom_hline(yintercept = PC4_thresh)
  print(c)
  
  #remove 1000G samples for the rest of the plots
  data <- data[data$set == "current",]
}
if(race_file != "none"){
  #regardless of whether 1000G data is here or not, if race for this dataset is present, plot it with everyone
  a <- ggplot(data = data, aes(x=PC1, y=PC2, color=race)) + geom_point()  +
    geom_vline(xintercept = PC1_thresh) + geom_hline(yintercept = PC2_thresh)
  print(a)
  b <- ggplot(data = data, aes(x=PC2, y=PC3, color=race)) + geom_point()  +
    geom_vline(xintercept = PC2_thresh) + geom_hline(yintercept = PC3_thresh)
  print(b)
  c <- ggplot(data = data, aes(x=PC3, y=PC4, color=race)) + geom_point()  +
    geom_vline(xintercept = PC3_thresh) + geom_hline(yintercept = PC4_thresh)
  print(c) 
   
  if(write_excl_file == "yes"){
    #get NHW and non-outliers (according to the current sample mean)
    ids_to_keep <- data[data$PC1 > PC1_thresh[1] & data$PC1 < PC1_thresh[2] &
                            data$PC2 > PC2_thresh[1] & data$PC2 < PC2_thresh[2] &
                            data$PC3 > PC3_thresh[1] & data$PC3 < PC3_thresh[2],
                          c("FID", "IID")]
    write.table(ids_to_keep, paste0(pc_file_stem, "_nooutliers.txt"), col.names = F, row.names = F, quote = F)
    
  }
} else {
  #if no race information is available, plot everyone with lines based on whole sample
  a <- ggplot(data = data, aes(x=PC1, y=PC2)) + geom_point()  +
    geom_vline(xintercept = PC1_thresh) + geom_hline(yintercept = PC2_thresh)
  print(a)
  b <- ggplot(data = data, aes(x=PC2, y=PC3)) + geom_point()  +
    geom_vline(xintercept = PC2_thresh) + geom_hline(yintercept = PC3_thresh)
  print(b)
  c <- ggplot(data = data, aes(x=PC3, y=PC4)) + geom_point()  +
    geom_vline(xintercept = PC3_thresh) + geom_hline(yintercept = PC4_thresh)
  print(c)
  
  if(write_excl_file == "yes"){
    #get non-outliers (according to the whole sample mean)
    ids_to_keep <- data[data$PC1 > PC1_thresh[1] & data$PC1 < PC1_thresh[2] &
                          data$PC2 > PC2_thresh[1] & data$PC2 < PC2_thresh[2] &
                          data$PC3 > PC3_thresh[1] & data$PC3 < PC3_thresh[2],
                        c("FID", "IID")]
    write.table(ids_to_keep, paste0(pc_file_stem, "_nooutliers.txt"), col.names = F, row.names = F, quote = F)
  }
  
}

dev.off()