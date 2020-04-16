args <- commandArgs(TRUE)
geno_file_stem <- args[1]
ref_ambig_file <- args[2]

library(data.table)

# read in the freq file
frq <- fread(paste0(geno_file_stem, ".frq"), header=TRUE, data.table=FALSE)
bim <- fread(paste0(geno_file_stem, ".bim"), header=F, data.table=FALSE)
bim <- bim[,c("V2", "V4")]
names(bim) <- c("SNP", "BP")
frq <- merge(frq, bim, by = "SNP")

# read in reference file
ref_annot <- fread(ref_ambig_file)
names(ref_annot) <- c("CHR", "BP", "ref", "alt", "maf")

# get all ambiguous SNPs
ambig_snps <- frq[(frq$A1=="A" & frq$A2=="T") | (frq$A1=="T" & frq$A2=="A") | (frq$A1=="G" & frq$A2=="C") | (frq$A1=="C" & frq$A2=="G") | (frq$A1=="a" & frq$A2=="t") | (frq$A1=="t" & frq$A2=="a") | (frq$A1=="g" & frq$A2=="c") | (frq$A1=="c" & frq$A2=="g"),]
ambig_snps$id <- paste(ambig_snps$CHR, ambig_snps$BP, sep = "_")
print(paste(nrow(ambig_snps), "ambiguous SNPs present in current dataset"))

#subset to those present in the reference panel
#not_present <- ambig_snps[!paste(ambig_snps$CHR, ambig_snps$BP, sep = "_") %in% paste(ref_annot$CHR, ref_annot$BP, sep = "_"),]
ambig_snps <- merge(ambig_snps, ref_annot, by = c("CHR", "BP"))
print(paste(nrow(ambig_snps), "ambiguous SNPs present in the reference panel"))

## get ambiguous SNPs which will need to be updated or dropped
#list and remove those with a MAF < 0.4
high_maf_ambig <- ambig_snps[ambig_snps$MAF>0.4,]
ambig_snps <- ambig_snps[ambig_snps$MAF<0.4,]
# need ref allele switch
SNPs_recode_ambig <- ambig_snps[ambig_snps$A1==ambig_snps$alt,c("SNP", "ref")]
# matched SNPS, non-ambiguous
SNPs_match_ambig <- as.character(ambig_snps$SNP[ambig_snps$A1==ambig_snps$ref_vcf])
#no match
SNPs_no_match <- ambig_snps[!ambig_snps$SNP %in% SNPs_match_ambig & !ambig_snps$SNP %in% SNPs_recode_ambig$SNP,]

#print numbers for ambiguous SNPs to log
print(paste(length(SNPs_match_ambig), "non-problematic ambiguous SNPs"))
print(paste(nrow(SNPs_recode_ambig), "ambiguous SNPs to force reference allele"))
print(paste0("Removing ",nrow(high_maf_ambig), " ambiguous SNPs with a minor allele frequency greater than 40%"))
print(paste(nrow(SNPs_no_match), "SNPs which do not match the reference panel."))

#write out
to_drop <- rbind(high_maf_ambig, SNPs_no_match)
write.table(to_drop[,c("SNP"), drop = F], paste0(geno_file_stem, "_ambig_to_drop.txt"), row.names = F, col.names = F, sep = " ", quote = F)

write.table(SNPs_recode_ambig, paste0(geno_file_stem, "_ambig_to_force_ref.txt"), row.names = F, col.names = F, sep = " ", quote = F)


