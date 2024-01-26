# load required packages
suppressMessages(library(SNPRelate))
suppressMessages(library(ggplot2))
suppressMessages(library(optparse))
suppressMessages(library(data.table))

## set up options
option_list = list(
  make_option(c("-i", "--inputgenotypes"), type="character", default=NULL, 
              help="full path and stem of input genotypes in plink format", metavar="character"),
  make_option(c("-l", "--label"), type="character", default=NULL, 
              help="short label for the cohort or 'no'", metavar="character"),
  make_option(c("-r", "--racefile"), type="character", default=NULL, 
              help="full path and name of the race/sex file", metavar="character"), 
  make_option(c("-R", "--racefile1000G"), type="character", default=NULL, 
              help="full path and name of the 1000G race/sex file", metavar="character"), 
  make_option(c("-o", "--output"), type = "character", default="no",
              help="indication to not output an exclusion file based on PC outlier status (>5 sd from the NHW mean)")
); 
## parse them
opt = parse_args(OptionParser(option_list=option_list))
inputgeno <- opt$inputgenotypes
datasetlabel <- opt$label
racefile <- opt$racefile
racefile1000G <- opt$racefile1000G
exclusion <- opt$output

## output files will have the same stem as the input genotypes

########################## convert plink to GDS #############################

print(paste0("Converting plink files to gds format... ", Sys.time(), "\n"))

#create the .gds file from plink file 
snpgdsBED2GDS(paste0(inputgeno, ".bed"), paste0(inputgeno, ".fam"), paste0(inputgeno, ".bim"),
                paste0(inputgeno, ".gds"),
                cvt.chr = "int", family = TRUE)

# now open it
genofile <- snpgdsOpen(paste0(inputgeno, ".gds"))

## decided against doing the pruning with the GDS file since it is much less efficient than plink
# # generate a pruned list of variants
# snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2, slide.max.n = 200)

########################## calculate PCs #############################

print(paste0("GDS file generated. Now calculating PCs... ", Sys.time(), "\n"))

# calculate PCs
pccalc <- snpgdsPCA(genofile, eigen.cnt = 10)

# extract just PC values to dt
pcs <- as.data.table(cbind(pccalc$sample.id, pccalc$eigenvect))
names(pcs) <- c("ID", paste0("PC", 1:10))

#update col classes
pcs[, paste0("PC", 1:10)] <- pcs[, lapply(.SD, as.numeric), .SDcols=paste0("PC", 1:10)]

# split FID/IID
pcs$FID <- sapply(strsplit(pcs$ID, ":"), "[", 1)
pcs$IID <- sapply(strsplit(pcs$ID, ":"), "[", 2)
pcs[,ID:=NULL]

# remove as much as possible from env
# rm(genofile, pccalc)

########################## determine outliers #############################

print(paste0("PCs calculated. Now calculating outliers and plotting... ", Sys.time(), "\n"))

if(!is.null(racefile)){
  
  #read in race
  data <- fread(racefile, data.table = F, header = F, select = 1:3)
  names(data) <- c("FID", "IID", "race")
  data$set <- "current"
  
  if(!is.null(racefile1000G)){
    
    #read in race
    data_1000G <- fread(racefile1000G, data.table = F, header = F)
    names(data_1000G) <- c("FID", "IID", "race")
    data_1000G$set <- "1000G"
    
    #combine
    data <- rbindlist(list(data, data_1000G))
    rm(data_1000G)
    
  }
  
  #merge race and PCs
  data <- merge(pcs, data, by = c("FID", "IID"))
  
  #get NHW subset (including 1000G and the current dataset)
  data_nhw <- data[data$race %in% c("EUR", "White"),]
  
  #set outlier thresholds based on NHW in 1000G and current dataset
  PC1_thresh <- c((mean(data_nhw$PC1)-5*sd(data_nhw$PC1)), (mean(data_nhw$PC1)+5*sd(data_nhw$PC1)))
  PC2_thresh <- c((mean(data_nhw$PC2)-5*sd(data_nhw$PC2)), (mean(data_nhw$PC2)+5*sd(data_nhw$PC2)))
  PC3_thresh <- c((mean(data_nhw$PC3)-5*sd(data_nhw$PC3)), (mean(data_nhw$PC3)+5*sd(data_nhw$PC3)))
  PC4_thresh <- c((mean(data_nhw$PC4)-5*sd(data_nhw$PC4)), (mean(data_nhw$PC4)+5*sd(data_nhw$PC4)))
  
  #add dataset label to this dataset's race categories (if the label was supplied) for the sake of a clearer legend
  if(!is.null(datasetlabel)){
    data$race[data$set == "current"] <- paste(datasetlabel, data$race[data$set == "current"], sep = " ")
  }
  
  #get non-outliers from just nhw in the current dataset
  ids_to_keep <- data_nhw[data_nhw$PC1 > PC1_thresh[1] & data_nhw$PC1 < PC1_thresh[2] &
                            data_nhw$PC2 > PC2_thresh[1] & data_nhw$PC2 < PC2_thresh[2] &
                            data_nhw$PC3 > PC3_thresh[1] & data_nhw$PC3 < PC3_thresh[2] &
                            data_nhw$set == "current",
                          c("FID", "IID")]
  if(exclusion == "yes"){
    write.table(ids_to_keep, paste0(inputgeno, "_NHW_nooutliers.txt"), col.names = F, row.names = F, quote = F)
    print(paste0("List of IDs for samples who are not outliers written to ", inputgeno, "_NHW_nooutliers.txt.\n"))
  }
  
  #remove NHW dataframe since it's not needed anymore
  rm(data_nhw)
  
} else {
  
  #if there's no race, just use PCs to determin outlier thresholds
  data <- pcs
  
  #set outlier thresholds based on whole sample
  PC1_thresh <- c((mean(data$PC1)-5*sd(data$PC1)), (mean(data$PC1)+5*sd(data$PC1)))
  PC2_thresh <- c((mean(data$PC2)-5*sd(data$PC2)), (mean(data$PC2)+5*sd(data$PC2)))
  PC3_thresh <- c((mean(data$PC3)-5*sd(data$PC3)), (mean(data$PC3)+5*sd(data$PC3)))
  PC4_thresh <- c((mean(data$PC4)-5*sd(data$PC4)), (mean(data$PC4)+5*sd(data$PC4)))
  
  #get non-outliers
  ids_to_keep <- data[data$PC1 > PC1_thresh[1] & data$PC1 < PC1_thresh[2] &
                        data$PC2 > PC2_thresh[1] & data$PC2 < PC2_thresh[2] &
                        data$PC3 > PC3_thresh[1] & data$PC3 < PC3_thresh[2],
                      c("FID", "IID")]
  if(exclusion == "yes"){
    write.table(ids_to_keep, paste0(inputgeno, "_nooutliers.txt"), col.names = F, row.names = F, quote = F)
    print(paste0("List of IDs for samples who are not outliers written to ", inputgeno, "_nooutliers.txt.\n"))
  }
}

########################## plot #############################

#create plots
pdf(paste0(inputgeno, ".pdf"))

if(!is.null(racefile1000G) & !is.null(racefile)){
  
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

if(!is.null(racefile)){
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
  
  #plot in non-hispanic whites
  if(!is.null(datasetlabel)){
    data <- data[data$race == paste(datasetlabel, "White", sep = " "),]
  } else {
    data <- data[data$race == "White",]
  }
  if(nrow(data)>0){
    a <- ggplot(data = data, aes(x=PC1, y=PC2, color=race)) + geom_point()  +
      geom_vline(xintercept = PC1_thresh) + geom_hline(yintercept = PC2_thresh)
    print(a)
    b <- ggplot(data = data, aes(x=PC2, y=PC3, color=race)) + geom_point()  +
      geom_vline(xintercept = PC2_thresh) + geom_hline(yintercept = PC3_thresh)
    print(b)
    c <- ggplot(data = data, aes(x=PC3, y=PC4, color=race)) + geom_point()  +
      geom_vline(xintercept = PC3_thresh) + geom_hline(yintercept = PC4_thresh)
    print(c)
    
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
  
}

dev.off()

print(paste0("PC calculation and plotting complete at ", Sys.time(), ". See ", inputgeno, ".pdf for PC plots to make decisions about outliers.\n\n"))