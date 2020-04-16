args <- commandArgs(TRUE)
relatedness_file_stem <- args[1]

library(data.table)

#read
genome <- fread(paste0(relatedness_file_stem, ".genome"), select = c("PI_HAT", "FID1", "IID1", "FID2", "IID2"))

#list pairs of individuals with PI_HAT > 0.90 because they are probably sample mix-ups
dups <- genome[genome$PI_HAT >= 0.90,] 
dups1 <- dups[,c("FID1","IID1")]
names(dups1) <- c("FID","IID")
dups2 <- dups[,c("FID2","IID2")]
names(dups2) <- c("FID","IID")
dups <- rbind(dups1,dups2)
dups <- unique(dups) 

#pull out one individual from each pair with PI_HAT < 0.90 & >0.25
top <- genome[genome$PI_HAT > 0.25 & genome$PI_HAT < 0.90,]
top1 <- top[,c("FID1","IID1")]
top2 <- top[,c("FID2","IID2")]

#Combine the individuals with PI_HAT>0.90 with
#the smallest number of unique individuals from each pair with PI_HAT < 0.90 & >0.25
if ( length(unique(top1$FID)) < length(unique(top2$FID)) ) {
  names(top1) <- c("FID","IID")
  related <- rbind(dups,top1)
} else {
  names(top2) <- c("FID","IID")
  related <- rbind(dups,top2)
}

#Get unique IDs
related <- unique(related)

#save
write.table(related,paste0(relatedness_file_stem, "_related_ids.txt"), col.names = FALSE,row.names=FALSE,quote = FALSE,sep="\t")