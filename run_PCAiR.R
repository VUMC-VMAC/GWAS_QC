# PC AiR
args <- commandArgs(TRUE)
input_geno <- args[1]

library(GWASTools)
library(SNPRelate)
library(GENESIS)
library(ggplot2)
library(GGally)

# outline:
## Kinship calculation
## PC and kinship calculation
## Second round of PC and kinship calculation
## Remove outliers and write out the final PC and kinship calculations

data.table::setDTthreads(4)

# required inputs
output_stem <- paste0(input_geno, "_PCAiR")

#Create the .gds file
snpgdsBED2GDS(paste0(input_geno, ".bed"),
              paste0(input_geno, ".fam"), 
              paste0(input_geno, ".bim"),
              paste0(input_geno, ".gds"),
              cvt.chr = "int")


# now open it
genofile <- snpgdsOpen(paste0(input_geno, ".gds"))

############ First, run IBD calculation ##############
#Run the ibd calculation
ibd.robust <- snpgdsIBDKING(genofile)

#get the kinship matrix
KINGmat <- ibd.robust$kinship

#set the col/row names to be the sample names
row.names(KINGmat) <- ibd.robust$sample.id
colnames(KINGmat) <- ibd.robust$sample.id

pdf(paste0(output_stem, "_rel0.pdf"))
#Make plot (with kinship lines)
kinship <- snpgdsIBDSelection(ibd.robust)
ggplot(kinship, aes(IBS0, kinship)) +
  geom_hline(yintercept=2^(-seq(3,9,2)/2), linetype="dashed", color="grey") +
  geom_point(alpha=0.5) +
  ylab("kinship estimate") +
  theme_bw()
dev.off()

#Close the snpgdsOpen connection since PC-air wants it another way
snpgdsClose(genofile)

############ Next, calculate PCs using those relatedness estimates and recalculate relatedness ##############

############### PCs
#now use the other function to open the connection
genofile <- GdsGenotypeReader(paste0(input_geno, ".gds"))
genodata <- GenotypeData(genofile)

#Run PC-AiR, setting the kinship threshold to 
PC_calc1 <- pcair(gdsobj =genodata, kinobj = KINGmat, divobj = KINGmat, kin.thresh = 2^(-7/2), div.thresh = -2^(-7/2), num.cores = 2)

#Get pcs in df
pcs <- as.data.frame(PC_calc1$vectors)
names(pcs) <- paste0("PC", 1:32)
pcs$FID <- row.names(pcs)


#Get ancestry informative PCs
anc_inf_pcs <- pcs[,c("FID", "PC1", "PC2", "PC3")]
row.names(anc_inf_pcs) <- anc_inf_pcs$FID
anc_inf_pcs$FID <- NULL
anc_inf_pcs <- as.matrix(anc_inf_pcs)

############### recalculate relatedness

#Create the genotype iterator object, can specify the number of variants per block here with the snpBlock argument
iterator <- GenotypeBlockIterator(genodata)

# recalculate kinship
new_kinship <- pcrelate(iterator, pcs = anc_inf_pcs)

#Make kinship matrix using built-in function
## Default to be aware of scaleKin = 2
relatedness_matrix <- pcrelateToMatrix(new_kinship)

############### Plot those PCs and relatedness estimates


pdf(paste0(output_stem, "_pcs1.pdf"))

# make plot of the first 12
ggparcoord(pcs, columns=1:12, scale="uniminmax") +
  xlab("PC") + ylab("")

# make plots of the first 4 without thresholds and without color since this is a first-pass
a <- ggplot(data = pcs, aes(x=PC1, y=PC2)) + geom_point()
print(a)
b <- ggplot(data = pcs, aes(x=PC2, y=PC3)) + geom_point() 
print(b)
c <- ggplot(data = pcs, aes(x=PC3, y=PC4)) + geom_point() 
print(c)

dev.off()


pdf(paste0(output_stem, "_rel1.pdf"))
#Make plot (with kinship lines)
ggplot(new_kinship$kinBtwn, aes(k0, kin)) +
  geom_hline(yintercept=2^(-seq(3,9,2)/2), linetype="dashed", color="grey") +
  geom_point(alpha=0.5) +
  ylab("kinship estimate") +
  theme_bw()
dev.off()



############ Then, recalculate PCs and relatedness ############


############### recalculate PCs for the last time

#Now recalculate the PCs using these kinship calculations
new_pcair <- pcair(gdsobj =genodata, kinobj =relatedness_matrix, divobj =  KINGmat, kin.thresh = 2^(-7/2), div.thresh = -2^(-7/2), num.cores = 2)

#get pcs in df
pcs <- as.data.frame(new_pcair$vectors)
names(pcs) <- paste0("PC", 1:32)
pcs$FID <- row.names(pcs)

#get ancestry informative PCs
anc_inf_pcs <- pcs[,c("FID", "PC1", "PC2", "PC3")]
row.names(anc_inf_pcs) <- anc_inf_pcs$FID
anc_inf_pcs$FID <- NULL
anc_inf_pcs <- as.matrix(anc_inf_pcs)

############### recalculate relatedness for the last time

#recalculate relatedness again, using these pcs
new_kinship2 <- pcrelate(iterator, pcs = anc_inf_pcs, training.set = new_pcair$unrels, sample.include = row.names(anc_inf_pcs), num.cores = 2)

############ Finally, make final plots and write out corrected PCs and relatedness ############

# will color the PC plots on kinship groups to see whether the PCs are picking up on more distant ancestry rather than kinship groups
not_related <- c(new_kinship$kinBtwn$ID1[new_kinship$kinBtwn$kin<0.04419], new_kinship$kinBtwn$ID2[new_kinship$kinBtwn$kin <0.04419])
fourth_degree <- c(new_kinship$kinBtwn$ID1[new_kinship$kinBtwn$kin>0.04419], new_kinship$kinBtwn$ID2[new_kinship$kinBtwn$kin>0.04419])
third_degree <- c(new_kinship$kinBtwn$ID1[new_kinship$kinBtwn$kin>0.088388], new_kinship$kinBtwn$ID2[new_kinship$kinBtwn$kin>0.088388])
second_degree <- c(new_kinship$kinBtwn$ID1[new_kinship$kinBtwn$kin>0.17678], new_kinship$kinBtwn$ID2[new_kinship$kinBtwn$kin>0.17678])
first_degree <- c(new_kinship$kinBtwn$ID1[new_kinship$kinBtwn$kin>0.35355], new_kinship$kinBtwn$ID2[new_kinship$kinBtwn$kin>0.35355])
pcs$related_groups[pcs$FID %in% not_related] <- "not_related"
pcs$related_groups[pcs$FID %in% fourth_degree] <- "fourth_degree"
pcs$related_groups[pcs$FID %in% third_degree] <- "third_degree"
pcs$related_groups[pcs$FID %in% second_degree] <- "second_degree"
pcs$related_groups[pcs$FID %in% first_degree] <- "first_degree"

# set outlier thresholds 
PC1_thresh <- c((mean(pcs$PC1)-5*sd(pcs$PC1)), (mean(pcs$PC1)+5*sd(pcs$PC1)))
PC2_thresh <- c((mean(pcs$PC2)-5*sd(pcs$PC2)), (mean(pcs$PC2)+5*sd(pcs$PC2)))
PC3_thresh <- c((mean(pcs$PC3)-5*sd(pcs$PC3)), (mean(pcs$PC3)+5*sd(pcs$PC3)))
PC4_thresh <- c((mean(pcs$PC4)-5*sd(pcs$PC4)), (mean(pcs$PC4)+5*sd(pcs$PC4)))


# make PC plots
pdf(paste0(output_stem, "_pcs2.pdf"))

## plot of the first 12 all together
ggparcoord(pcs, columns=1:12, scale="uniminmax") +
  xlab("PC") + ylab("")

## plots of the first 4 with 
ggplot(pcs, aes(PC1, PC2, color=related_groups)) + geom_point() +
  geom_vline(xintercept = PC1_thresh) + geom_hline(yintercept = PC2_thresh)

ggplot(pcs, aes(PC2, PC3, color=related_groups)) + geom_point() +
  geom_vline(xintercept = PC2_thresh) + geom_hline(yintercept = PC3_thresh) 

ggplot(pcs, aes(PC3, PC4, color=related_groups)) + geom_point() +
  geom_vline(xintercept = PC3_thresh) + geom_hline(yintercept = PC4_thresh)

dev.off()

#plot kinship estimates
pdf(paste0(output_stem, "_rel2.pdf"))
ggplot(new_kinship2$kinBtwn, aes(k0, kin)) +
  geom_hline(yintercept=2^(-seq(3,9,2)/2), linetype="dashed", color="grey") +
  geom_point(alpha=0.5) +
  ylab("kinship estimate") +
  theme_bw()
dev.off()

############### write out those PCs and relatedness calculations since they should be corrected enough

# save PC output
write.table(pcs, paste0(output_stem, "_pcs2.txt"), sep = " ", quote = F, row.names = F, col.names = T)

# save these relatedness estimates
save(new_kinship2, file = paste0(output_stem, "_rel2.RData"))

