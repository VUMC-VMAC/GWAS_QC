args <- commandArgs(TRUE)
file_stem <- args[1]

library(data.table)

#check ids in bim
bim <- fread(paste0(file_stem, ".bim"), data.table = F)

print("Checking variant and person id lengths for input into smartpca...")

#if there are too long ids, then try to create sensibly named shorter ids
if(sum(sapply(bim$V2, nchar)>30)>1){
  print("Too long variant ids present in the bim file. Attempting to update.")
  to_update_bim <- bim[sapply(bim$V2, nchar)>30,]
  to_update_bim$new_id <- paste(to_update_bim$V1, to_update_bim$V4, sep = "_")
  
  #kill script if this doesn't work
  if(sum(sapply(to_update_bim$new_id, nchar)>30)>0){
    stop("IDs in bim file are too long and cannot create simple shorter ids! Please fix!")
  }
  
  #if it works write it out and tell the audience
  to_update_bim <- to_update_bim[,c("V2", "new_id")]
  write.table(to_update_bim, paste0(file_stem, "_shorter_bimids.txt"), row.names = F, col.names =  F, sep = " ", quote = F)
} else {
  print("Variant ids ok.")
}

#check ids in fam
fam <- fread(paste0(file_stem, ".fam"), data.table = F)

#if there are too long ids, then create a hash file with dummy ids
if(sum(sapply(fam$V1, nchar)>30 | sapply(fam$V2, nchar)>30)>1){
  
  print("Too long sample ids present in the fam file. Updating all with dummy ids for PC calculation...")
  fam$FID_dummy <- 1:nrow(fam)
  fam$IID_dummy <- 1:nrow(fam)
  
  fam <- fam[,c("V1", "V2", "FID_dummy", "IID_dummy")]
  write.table(fam, paste0(file_stem, "_dummy_famids.txt"), row.names = F, col.names =  F, sep = " ", quote = F)
} else {
  print("Sample ids ok.")
}


