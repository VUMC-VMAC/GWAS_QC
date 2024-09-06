args <- commandArgs(TRUE)
inferanc <- args[1] #name of the inferred ancestry file with full path

suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(forcats))

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
anc$Ancestry <- "Other"
anc$Ancestry[anc$EUR >= 0.75] <- "EUR"
anc$Ancestry[anc$AFR_EUR >= 0.80 & !(anc$Ancestry %in% c("EUR"))] <- "AFR"
anc$Ancestry[anc$AMR >= 0.60 & anc$AFR <= 0.05 & !(anc$Ancestry %in% c("AFR","EUR"))] <- "AMR"
anc$Ancestry[anc$AMR_EUR >= 0.80 & anc$AFR <= 0.10 & !(anc$Ancestry %in% c("AMR","EUR"))] <- "AMR2admx"
anc$Ancestry[anc$AFR_EUR_AMR >= 0.85 & anc$AFR <= 0.45 & !(anc$Ancestry %in% c("AMR2admx", "AMR", "AFR","EUR"))] <- "AMR3admx"

# print out how many in each group
print(paste("There are", sum(anc$Ancestry == "Other"), "other admixed,", sum(anc$Ancestry == "EUR"), "EUR,", sum(anc$Ancestry == "AFR"), "AFR,", sum(anc$Ancestry == "AMR"), "AMR,", sum(anc$Ancestry == "AMR2admx"), "2-way admixed AMR, and", sum(anc$Ancestry == "AMR3admx"), "3-way admixed AMR individuals based on genetic similarity to 1000G reference superpopulations."))

# split ID back into FID and IID
anc$FID <- sapply(strsplit(anc$ID, ":"), "[", 1)
anc$IID <- sapply(strsplit(anc$ID, ":"), "[", 2)

# split into separate df
afr <- anc[anc$Ancestry == "AFR",c("FID", "IID")]
eur <- anc[anc$Ancestry == "EUR",c("FID", "IID")]
amr <- anc[anc$Ancestry == "AMR",c("FID", "IID")]
amr2admx <- anc[anc$Ancestry == "AMR2admx",c("FID", "IID")]
amr3admx <- anc[anc$Ancestry == "AMR3admx",c("FID", "IID")]
other <- anc[anc$Ancestry == "Other",c("FID", "IID")]

# write out lists of ids to keep
write.table(other, paste0(filestem, "_other_keep.txt"), quote=F, row.names=F, col.names=F, sep="\t")
write.table(afr, paste0(filestem, "_AFR_keep.txt"), quote=F, row.names=F, col.names=F, sep="\t")
write.table(eur, paste0(filestem, "_EUR_keep.txt"),  quote=F, row.names=F, col.names=F, sep="\t")
write.table(amr, paste0(filestem, "_AMR_keep.txt"), quote=F, row.names=F, col.names=F, sep="\t")
write.table(amr2admx, paste0(filestem, "_AMR2admx_keep.txt"),  quote=F, row.names=F, col.names=F, sep="\t")
write.table(amr3admx, paste0(filestem, "_AMR3admx_keep.txt"),  quote=F, row.names=F, col.names=F, sep="\t")
print(paste0("Lists of IDs in AFR, EUR, AMR, AMR2admx, and AMR3admx subsets have been written to ", filestem, "_*_keep.txt. Use these to subset the cohort before finishing post-imputation QC."))


###################################################
##########Create Ancestry Estimates Graph##########
###################################################
#Select variables
ancplot <- anc %>%
  dplyr::rename(percent_SAS = SAS, percent_EAS = EAS, percent_AMR = AMR, percent_AFR = AFR, percent_EUR = EUR) %>%
  dplyr::select("IID", "percent_SAS", "percent_EAS", "percent_AFR", "percent_AMR", "percent_EUR", "Ancestry")

#Create variables for graphs
ancplot <- ancplot %>% 
  pivot_longer(c("percent_EUR", "percent_AFR", "percent_AMR", "percent_SAS", "percent_EAS"), names_to = "population")

ancplot$population <- factor(ancplot$population, levels = c("percent_EUR", "percent_AFR", "percent_AMR", "percent_SAS", "percent_EAS"))

ancplot <- ancplot %>%
  mutate(
    Ancestry = 
      case_when(
        Ancestry == "EUR"  ~ 1, 
        Ancestry == "AFR" ~ 2,
	Ancestry == "AMR" ~ 3,
        Ancestry == "AMR2admx" ~ 4,
        Ancestry == "AMR3admx" ~ 5,
        Ancestry == "Other" ~ 6
      )
  )

#Create ancestry names for plot
ancestry_names <- c(
  '1'="EUR",
  '2'="AFR",
  '3'="AMR",
  '4'="AMR2admx",
  '5'="AMR3admx",
  '6'="Other"
)
#ancestry_labeller <- function(variable,value){
#  return(ancestry_names[value])
#}

#Create plot with ancestry estimates
options(bitmapType='cairo')
png(paste0(filestem, ".ancestry_estimates_graph.png"), width = 1200, height = 720)
ancplot %>%
  arrange(population, value) %>%
  mutate(IID = forcats::fct_inorder(IID)) %>%
  ggplot() +
  geom_col(aes(y = reorder(IID, Ancestry), x = value, fill = population, group = Ancestry), width = 1) +
  scale_x_continuous(labels = scales::percent_format(scale = 100), name = "Estimates (%)") +
  scale_y_discrete(name = "Individuals") +
  coord_flip() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.spacing.x=unit(0.025, "lines"),
        strip.text.x = element_text(angle = 270),
        panel.grid = element_blank()) +
  facet_grid(~Ancestry, scales="free_x", space="free_x", switch="x", drop = TRUE, labeller = labeller(Ancestry = ancestry_names))
dev.off()
