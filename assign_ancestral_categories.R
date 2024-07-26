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
anc$European <- ifelse(anc$EUR >= 0.75, 1, 0)
anc$AA <- ifelse(anc$AFR_EUR >= 0.80 & anc$European == 0, 1, 0)
anc$NatAmer <- ifelse(anc$AMR >= 0.60 & anc$AFR <= 0.05 & anc$AA == 0 & anc$European == 0, 1, 0)
anc$Lat_Hisp <- ifelse(anc$AMR_EUR >= 0.80 & anc$AFR <= 0.10 & anc$NatAmer == 0 & anc$European == 0, 1, 0)
anc$CaribHisp <- ifelse(anc$AFR_EUR_AMR >= 0.85 & anc$AFR <= 0.45 & anc$Lat_Hisp == 0 & anc$NatAmer == 0 & anc$AA == 0 & anc$European == 0, 1, 0)

# make one ancestry column
anc$Ancestry <- "Admixed"
anc$Ancestry[anc$European=="1"] <- "European"
anc$Ancestry[anc$AA=="1"] <- "AA"
anc$Ancestry[anc$NatAmer=="1"] <- "Native_American"
anc$Ancestry[anc$Lat_Hisp=="1"] <- "Latino_Hispanic"
anc$Ancestry[anc$CaribHisp=="1"] <- "Caribbean_Hispanic"

# print out how many in each group
print(paste("There are", sum(anc$Ancestry == "Admixed"), "admixed individuals,", sum(anc$Ancestry == "European"), "European,", sum(anc$Ancestry == "AA"), "African American,", sum(anc$Ancestry == "Native_American"), "Native American,", sum(anc$Ancestry == "Latino_Hispanic"), "Latino/Hispanic, and", sum(anc$Ancestry == "Caribbean_Hispanic"), "Caribbean Hispanic ancestry."))

# split ID back into FID and IID
anc$FID <- sapply(strsplit(anc$ID, ":"), "[", 1)
anc$IID <- sapply(strsplit(anc$ID, ":"), "[", 2)

# split into separate df
aa <- anc[anc$Ancestry == "AA",c("FID", "IID")]
eur <- anc[anc$Ancestry == "European",c("FID", "IID")]
nat <- anc[anc$Ancestry == "Native_American",c("FID", "IID")]
lat <- anc[anc$Ancestry == "Latino_Hispanic",c("FID", "IID")]
carib <- anc[anc$Ancestry == "Caribbean_Hispanic",c("FID", "IID")]

# write out lists of ids to keep
write.table(aa, paste0(filestem, "_AA_keep.txt"), quote=F, row.names=F, col.names=F, sep="\t")
write.table(eur, paste0(filestem, "_EUR_keep.txt"),  quote=F, row.names=F, col.names=F, sep="\t")
write.table(nat, paste0(filestem, "_NatAmer_keep.txt"), quote=F, row.names=F, col.names=F, sep="\t")
write.table(lat, paste0(filestem, "_LatHisp_keep.txt"),  quote=F, row.names=F, col.names=F, sep="\t")
write.table(carib, paste0(filestem, "_CaribHisp_keep.txt"),  quote=F, row.names=F, col.names=F, sep="\t")
print(paste0("Lists of IDs in AA, EUR, Native American, Latino/Hispanic, and Caribbean Hispanic subsets have been written to ", filestem, "_*_keep.txt. Use these to subset the cohort before finishing post-imputation QC."))


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
        Ancestry == "European"  ~ 1, 
        Ancestry == "AA" ~ 2,
	Ancestry == "Native_American" ~ 3,
        Ancestry == "Latino_Hispanic" ~ 4,
        Ancestry == "Caribbean_Hispanic" ~ 5,
        Ancestry == "Admixed" ~ 6
      )
  )

#Create ancestry names for plot
ancestry_names <- c(
  '1'="European",
  '2'="AA",
  '3'="NatAmer",
  '4'="LatHisp",
  '5'="CaribHisp",
  '6'="Admixed"
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
