###Script to find A1 matches and MAF differences < 0.1 for any number of datasets

#arguments
args <- commandArgs(TRUE)
dataset <- args[1] #name of dataset
bim_files <- args[2:length(args)] #path with file name for bim files (including .bim)

#packages
library(data.table)

#limit the threads it can use
setDTthreads(4)

#read in frq files
bim <- lapply(bim_files, function(x) fread(x, select = c(2,5)))

#merge into one df, only keeping overlapping variants
bim <- Reduce(function(...) merge(..., by = c("V2", "V5"), all = F), bim)
print(paste0("The number of variants with A1 matches is: ", nrow(bim)))

#writing out the variant column the df
write.table(bim$V2, paste0(dataset,"_A1_matches.txt"), row.names=F,col.names=F,quote=F)
