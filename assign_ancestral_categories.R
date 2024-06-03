args <- commandArgs(TRUE)
inferanc <- args[1] #name of the inferred ancestry file with full path

library(data.table)

# get the stem without the extension
filestem <- gsub(".out", "", inferanc)

# read in inferred ancestry
anc <- fread(inferanc)

# label column labels
names(anc) <- c("ID", "Pop_label", "N_SNPs", "pred_PC1", "pred_PC2", "pred_PC3", "pred_PC4", "SAS", "EAS", "AMR", "AFR", "EUR")

# generate additive columns
anc$AFR_EUR <- anc$AFR + anc$EUR
anc$AMR_EUR <- anc$AMR + anc$EUR
anc$AFR_EUR_AMR <- anc$AFR + anc$EUR + anc$AMR

# categorize based on thresholds
anc$European <- ifelse(anc$EUR >= 0.75, 1, 0)
anc$AA <- ifelse(anc$AFR_EUR >= 0.80 & anc$European == 0, 1, 0)
anc$Lat_Hisp <- ifelse(anc$AMR_EUR >= 0.80 & anc$AFR <= 0.10 & anc$European == 0, 1, 0)
anc$CaribHisp <- ifelse(anc$AFR_EUR_AMR >= 0.85 & anc$AA <= 0.3 & anc$Lat_Hisp == 0 & anc$AA == 0 & anc$European == 0, 1, 0)

# make one ancestry column
anc$Ancestry <- "Admixed"
anc$Ancestry[anc$European=="1"] <- "European"
anc$Ancestry[anc$AA=="1"] <- "AA"
anc$Ancestry[anc$Lat_Hisp=="1"] <- "Latino_Hispanic"
anc$Ancestry[anc$CaribHisp=="1"] <- "Caribbean_Hispanic"

# print out how many in each group
print(paste("There are", sum(anc$Ancestry == "Admixed"), "admixed individuals,", sum(anc$Ancestry == "European"), "European,", sum(anc$Ancestry == "AA"), "African American,", sum(anc$Ancestry == "Latino_Hispanic"), "Latino/Hispanic, and", sum(anc$Ancestry == "Caribbean_Hispanic"), "Caribbean Hispanic ancestry."))

# split ID back into FID and IID
anc$FID <- sapply(strsplit(anc$ID, ":"), "[", 1)
anc$IID <- sapply(strsplit(anc$ID, ":"), "[", 2)

# split into separate df
aa <- anc[anc$Ancestry == "AA",c("FID", "IID")]
eur <- anc[anc$Ancestry == "European",c("FID", "IID")]
lat <- anc[anc$Ancestry == "Latino_Hispanic",c("FID", "IID")]
carib <- anc[anc$Ancestry == "Caribbean_Hispanic",c("FID", "IID")]

# write out lists of ids to keep
write.table(aa, paste0(filestem, "_AA_keep.txt"), quote=F, row.names=F, col.names=F, sep="\t")
write.table(eur, paste0(filestem, "_EUR_keep.txt"),  quote=F, row.names=F, col.names=F, sep="\t")
write.table(lat, paste0(filestem, "_LatHisp_keep.txt"),  quote=F, row.names=F, col.names=F, sep="\t")
write.table(carib, paste0(filestem, "_CaribHisp_keep.txt"),  quote=F, row.names=F, col.names=F, sep="\t")
print(paste0("Lists of IDs in AA, EUR, Latino/Hispanic, and Caribbean Hispanic subsets have been written to ", filestem, "_*_keep.txt. Use these to subset the cohort before finishing post-imputation QC."))
