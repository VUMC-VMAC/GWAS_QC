###Script to find A1 matches and MAF differences < 0.1 for any number of datasets

#arguments
args <- commandArgs(TRUE)
dataset <- args[1] #name of dataset
frq_files <- args[2:length(args)] #path with file name for frq files (including .frq)

#packages
library(data.table)

#read in frq files
frq <- lapply(frq_files, fread)

#subset to just SNP and MAF cols
frq <- lapply(frq, function(x) x[,c("SNP","A1","MAF")])
frq_maf_names <- paste0("MAF", 1:length(frq))
frq <- lapply(seq(frq_maf_names), function(j) {
  names(frq[[j]])[3] <- frq_maf_names[j]
  return(frq[j])
})

#merge into one df, only keeping overlapping variants
frq <- Reduce(function(...) merge(..., by = c("SNP", "A1"), all = F), frq)
print(paste0("The number of variants with A1 matches is: ", nrow(frq)))

#check that mafs do not differ by more than 0.1 for these variants
frq_combos <- combn(frq_maf_names, 2)
for(i in 1:ncol(frq_combos)){
  frq1 <- frq_combos[1,i]
  frq2 <- frq_combos[2,i]
  MAF_diff <- frq[,frq1] - frq[,frq2]
  frq <- frq[abs(MAF_diff)<0.1,]
}
print(paste0("The number of variants with A1 matches and MAF differences less than 0.1 is: ", nrow(frq)))

#writing out the variant column the df
write.table(frq$SNP, paste0(dataset,"_A1_MAF_matches.txt"), row.names=F,col.names=F,quote=F)
