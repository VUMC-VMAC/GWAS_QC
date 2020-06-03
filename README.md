# GWAS QC and TOPMed Imputation Pipeline Documentation

Here are the general steps in this pipeline:
1.	Pre-Imputation QC Part 1
2.	Pre-Imputation QC Part 2
3.	TOPMed Imputation
4.	Post-Imputation QC

Each step (except the TOPMed imputation) is run as a shell script in a singularity. The usage message for each of the scripts can be shown by running the script with the -h flag. 

Depending on the size of the input dataset, it may be convenient to run the QC scripts as slurm jobs. There is an example slurm script for running the pre-imputation part 2 QC in this repository: run_gwas_qc_part2.slurm

Software used for QC:
- plink, versions 1.9 and 2.0: https://www.cog-genomics.org/plink/
- smartpca from Eigensoft: https://www.hsph.harvard.edu/alkes-price/software/ 
- R 3.6.3: https://www.r-project.org/
- bcftools: http://samtools.github.io/bcftools/bcftools.html
- htslib: http://www.htslib.org/ 

To be able to upload files for imputation to the TOPMed reference panel, please create an account on the TOPMed imputation server: https://imputation.biodatacatalyst.nhlbi.nih.gov/#!pages/home 

Supplementary files not included in this github:
- 1000G genetic data: https://www.internationalgenome.org/data/
- TOPMed SNP information files
- Singularity to run the GWAS QC

## Pre-Imputation QC Part 1

Files needed for the Pre-Imputation QC Part 1 script:
-	Raw genotype files 
-	Singularity container (current version: gwas_qc_singularity_v2.1.simg)
-	Files with FID, IID, race, and sex for each set
  - The race column should have non-Hispanic whites coded as "White" and none of the values can have spaces. Sex should be coded as 
-	1000G files to use in calculating PCs (optional)


To run the part 1 script, run a singularity command similar to the one below, binding the paths to the necessary input files. This script will take between 40 minutes to an hour to run, depending on the size of your dataset. In the following examples, COHORT will be the stem for all the input and output files with rest of the stems 
```
singularity exec --containall --bind /path/to/genotypes/:/inputs/ \
	    --bind /path/to/1000G/data/:/inputs2/ \
	    /data/h_vmac/GWAS_QC/singularity/gwas_qc_singularity_v2.1.simg \
	    /bin/bash -c "cd $SCRIPTS_DIR ; sh gwas_qc_part1.sh \
	    	      -i /inputs/COHORT_raw_genotypes \
		      	 -o /inputs/COHORT_cleaned \
			       -r /inputs/COHORT_race_sex.txt \
			       -G /inputs2/1000G_filestem |& tee /inputs/COHORT_gwas_qc_part1.log"
```

Pre-Imputation QC Part 1 output files:
-	Main log file
-	Heterozygosity plot: file ending in *_pruned_hetcheck.png
-	PC plots: 
  - 2 pdfs (if you supplied 1000G data), one ending in *_1000G_merged_geno01_pruned.pdf and the other one ending in *_pruned.pdf (no 1000G)

Once the part 1 script completes, view the resulting PC plots and make decisions regarding PC outliers. There will be a file output by default that includes on non-Hispanic white ids of individuals who did not fall more than 5 standard deviations from the mean of 1000G European samples and whites in the input dataset. If you decide to use this default inclusion file, run a plink command to keep those samples (using the --keep plink flag).


## Pre-Imputation QC Part 2

We run the part 2 pre-imputation QC and every step following in samples of all races as well as in just non-Hispanic whites.

Files needed for the Pre-Imputation QC Part 2 script:
-	Singularity container (current version: gwas_qc_singularity_v2.1.simg)
-	Genotype files resulting from the part 1 script
- Reference panel files


To run the part 2 script, run a command similar to the one below:
```
singularity exec --containall \
  --bind /scratch/mahone1/GWAS_QC_data/BIOCARD_recap/:/inputs/ \
  --bind /data/h_vmac/GWAS_QC/topmed/:/ref/ \
  /data/h_vmac/GWAS_QC/singularity/gwas_qc_singularity_v2.1.simg \
  /bin/bash -c "cd $SCRIPTS_DIR ; sh gwas_qc_part2.sh \
  -i /inputs/COHORT_genotyped_geno05_maf01_mind01_norelated_nomismatchedsex_keep_autosomes_NHW\
  -o /inputs/COHORT_NHW_genotyped_cleaned \
  -R /ref/topmed_imputed_ref_snps -n"

```

Note the -n flag which specifies to not exclude variants based on the reference panel. Additionally, one could add the -b flag followed by the genome build of the input dataset (ie b36, b37, or b38); the default is b37.

Pre-Imputation QC Part 2 primary output files:
-	Main log file
-	.vcf.gz files for each chromosome


## TOPMed Imputation

Upload the .vcf.gz files from the part 2 script to the TOPMed Imputation Server (https://imputation.biodatacatalyst.nhlbi.nih.gov/#!). Once imputation is complete, download the results and save the decryption password (received in an email from the imputation server) to a file called pass.txt in the folder with the imputation results. 
Alternatively, the imputation results can be unzipped before running the post-imputation QC script as long as the -z flag is excluded from the command.

## Post-Imputation QC

Files needed for the Pre-Imputation QC Part 2 script:
-	Singularity container (current version: gwas_qc_singularity_v2.1.simg)
-	Imputation results and password to unzip them
- File with race and sex for the dataset (used in the pre-imputation part 1 step)
- Files to convert variant ids back to rs number

To run the post-imputation QC, run a command similar to the one below.
```
singularity exec --contain --bind /path/to/genotype/data/:/inputs/ \
  --bind /path/to/rsnumber/conversion/files/:/ref/ \
  /data/h_vmac/GWAS_QC/singularity/gwas_qc_singularity_v2.1.simg \
  /bin/bash -c  "cd $SCRIPTS_DIR ; sh gwas_qc_postimputation.sh -z \
  -i /inputs/Imputed/Raw/ \
  -o /inputs/Imputed/Cleaned/COHORT_imputed_NHW \
  -r /inputs/COHORT_race_sex.txt \
  -s /ref/topmed_snp_names |& tee /inputs/Imputed/Cleaned/COHORT_NHW_qc_postimputation.log"
```

Post-imputation QC main output files:
- Main log file (COHORT_NHW_qc_postimputation.log)
- PCs and PC plots

After the post-imputation QC runs, check the principal components for outliers, removing any if necessary and recalculating PCs. 

To recalculate PCs, run a command similar to the one below. Adding the -n flag will suppress the creation of a file for default inclusions (ie only whites who fall within 5 SD from the mean).
```
singularity exec --contain --bind /path/to/genotype/data/:/inputs/ \
  /data/h_vmac/GWAS_QC/singularity/gwas_qc_singularity_v1.3.simg \
  /bin/bash -c  "cd $SCRIPTS_DIR ; sh calc_plot_PCs.sh \
  -i /inputs/Imputed/Raw/COHORT_imputed_NHW_ids_sex_maf01_hwe6_ids"
```
