###Script to find A1 matches and MAF differences < 0.1 for any number of datasets

#arguments
args <- commandArgs(TRUE)
dataset <- args[1] #name of dataset
bim_files <- args[2:length(args)] #path with file name for bim files (including .bim)

#packages
library(data.table)

#limit the threads it can use
setDTthreads(4)

#read in bim files
bim <- lapply(bim_files, function(x) fread(x, select = c(2)))

#merge into one df, only keeping overlapping variants
bim <- Reduce(merge, bim)
print(paste0("The number of overlapping variants is: ", nrow(bim)))

#writing out the variant column the df
write.table(bim$V2, paste0(dataset,"_overlap.txt"), row.names=F,col.names=F,quote=F)
