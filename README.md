# GWAS QC and TOPMed Imputation Pipeline Documentation

Here are the general steps in this pipeline:
1. Pre-Imputation QC Part 1
   - Variant filtering excluding for missingness (>5%) and minor allele frequency (<0.01)
   - Sample missingness filtering (5%)
   - Relatedness check (removes 1 of each pair of second degree relatives and both of each pair with a pi-hat >0.9)
   - Sex check (removes individuals with discrepancies)
   - Heterozygosity check
   - PC calculation and plotting
2. Pre-Imputation QC Part 2
   - Hardy-Weinberg filter (p<1e-6)
   - Removing palindromic, multi-allelic, and duplicate variants
   - Lifts to build 38 if necessary
   - Prepares files for uploading to the imputation server
3. TOPMed Imputation
4. Post-Imputation QC
   - Variant filtering excluding for.
      - imputation quality (<0.8)
      - multi-allelic variants
   - Hardy-Weinberg equilibrium (p<1e-6) and minor allele frequency (<0.01)
   - Merging in genotyped variants with the imputed data
   - PC calculation


Each step (except the TOPMed imputation) is run as a shell script in a singularity. The usage message for each of the scripts can be shown by running the script with the -h flag.

Depending on the size of the input dataset, it may be convenient to run the QC scripts as slurm jobs. There is an example slurm script for running the pre-imputation part 2 QC in this repository: run_gwas_qc_part2.slurm

Software used for QC:
- plink, versions 1.9 and 2.0: https://www.cog-genomics.org/plink/
- smartpca from Eigensoft: https://www.hsph.harvard.edu/alkes-price/software/
- R 3.6.3: https://www.r-project.org/
- bcftools: http://samtools.github.io/bcftools/bcftools.html
- htslib: http://www.htslib.org/
- Imputation checking script from the McCarthy group: https://www.well.ox.ac.uk/~wrayner/tools/

To be able to upload files for imputation to the TOPMed reference panel, please create an account on the TOPMed imputation server: https://imputation.biodatacatalyst.nhlbi.nih.gov/#!pages/home

Supplementary files not included in this github:
- 1000G genetic data: https://www.internationalgenome.org/data/
- TOPMed SNP information files
- Singularity to run the GWAS QC

## Pre-Imputation QC Part 1

Files needed for the Pre-Imputation QC Part 1 script:
- Raw genotype files
- Singularity container (current version: CNT_genomic_processing_v2.simg)
- Files with FID, IID, race, and sex for each set
  - The race column should have non-Hispanic whites coded as "White" and none of the values can have spaces. Sex should be coded as 1=male, 2=female, 0=unknown.
- 1000G files to use in calculating PCs (optional)


To run the part 1 script, run a singularity command similar to the one below, binding the paths to the necessary input files. This script will take between 40 minutes to an hour to run, depending on the size of your dataset. As in the example below, -i specifies the input genotypes stem, -o specifies the output stem, -r specifies the race/sex file, and -G specifies the 1000G data.
```
singularity exec --containall --bind /nfs:/nfs --bind  /data/h_vmac/GWAS_QC/:/ref/ --bind  /scratch/:/scratch/   \ 
	    /data/h_vmac/GWAS_QC/singularity/CNT_genomic_processing_v2.simg /bin/bash -c "cd /scripts/GWAS_QC/ ; \
	    sh gwas_qc_part1.sh -i   /scratch/mahone1/GWAS_QC_data/A4/Genotyped/Raw/A4_ids_sex \
	    -o   /scratch/mahone1/GWAS_QC_data/A4/Genotyped/Cleaned/part1/A4_genotyped \
	    -r   /scratch/mahone1/GWAS_QC_data/A4/A4_race_sex.txt \
	    -G   /ref/1000G_data/1000G_cleaned |& tee /scratch/mahone1/GWAS_QC_data/A4/Genotyped/Cleaned/part1/A4_gwas_qc_part1.log"
```

Pre-Imputation QC Part 1 output files:
- Main log file
- Heterozygosity plot: file ending in *_pruned_hetcheck.png
- PC plots:
  - 2 pdfs (if you supplied 1000G data), one ending in *_1000G_merged_geno01_pruned.pdf and the other one ending in *_pruned.pdf (no 1000G)
- Cleaned genotype files, ending in *_keep_autosomes (unless there were heterozygosity outliers)

Once the part 1 script completes, view the resulting PC plots and make decisions regarding PC outliers. There will be a file output by default that includes on non-Hispanic white ids of individuals who did not fall more than 5 standard deviations from the mean of 1000G European samples and whites in the input dataset. If you decide to use this default inclusion file, run a plink command to keep those samples (using the --keep plink flag).

## Pre-Imputation QC Part 2

We run the part 2 pre-imputation QC and every step following in samples of all races as well as in just non-Hispanic whites.

Files needed for the Pre-Imputation QC Part 2 script:
- Singularity container (current version: CNT_genomic_processing_v2.simg)
- Genotype files resulting from the part 1 script
- Reference panel files


To run the part 2 script, run a command similar to the one below:
```
singularity exec --containall \
	    --bind /nfs:/nfs --bind  /data/h_vmac/GWAS_QC/:/ref/ --bind  /scratch/:/scratch/  \
	    /data/h_vmac/GWAS_QC/singularity/CNT_genomic_processing_v2.simg /bin/bash -c   "cd /scripts/GWAS_QC/ ;\
	    sh gwas_qc_part2.sh -i /scratch/mahone1/GWAS_QC_data/A4/Genotyped/Cleaned/part1/A4_genotyped_geno05_maf01_mind01_norelated_sex_nomismatchedsex_keep_autosomes \
	    -o /scratch/mahone1/GWAS_QC_data/A4/Genotyped/Cleaned/part2_allraces/A4_AllRaces_genotyped_cleaned \
	    -R /ref/topmed/topmed_imputed_ref_snps |& tee  /scratch/mahone1/GWAS_QC_data/A4/Genotyped/Cleaned/part2_allraces/A4_gwas_qc_part2_allraces.log"
```
The -n flag would specify to not exclude variants based on the reference panel. The -b flag can be used to specify genome build (b36, b37, b38); the default is b37.

Pre-Imputation QC Part 2 primary output files:
- Main log file
- .vcf.gz files for each chromosome


## TOPMed Imputation

Upload the .vcf.gz files from the part 2 script to the TOPMed Imputation Server (https://imputation.biodatacatalyst.nhlbi.nih.gov/#!). Once imputation is complete, download the results and save the decryption password (received in an email from the imputation server) to a file called pass.txt in the folder with the imputation results.
Alternatively, the imputation results can be unzipped before running the post-imputation QC script as long as the -z flag is excluded from the command.

## Post-Imputation QC

Files needed for the Pre-Imputation QC Part 2 script:
- Singularity container (current version: CNT_genomic_processing_v2.simg)
- Imputation results and password to unzip them
- File with race and sex for the dataset (used in the pre-imputation part 1 step)
  - Note that this *must* include every sample that is in the imputation files to be able to update FID/IID and sex correctly as well as to merge in the genotype data correctly.
- Files to convert variant ids back to rs number
- Cleaned genotype files to merge back in
  - This will be the last file from pre-imputation QC, before prep for the imputation server (ending in *_chr_bplifted_nosamepos or *_chr_bplifted, depending on whether there are multi-allelic or duplicated variants to remove)

To run the post-imputation QC, run a command similar to the one below.
```
singularity exec --contain --bind /scratch:/scratch --bind /data/h_vmac/GWAS_QC/topmed/:/ref/ \
	    /data/h_vmac/GWAS_QC/singularity/CNT_genomic_processing_v2.simg /bin/bash -c  "cd /scratch/mahone1/GWAS_QC/ ; \
	    sh gwas_qc_postimputation.sh -i /scratch/mahone1/GWAS_QC_data/A4/Imputed/Raw/ \
	    -o /scratch/mahone1/GWAS_QC_data/A4/Imputed/Cleaned/AllRaces/A4_AllRaces_imputed  \
	    -r /scratch/mahone1/GWAS_QC_data/A4/A4_race_sex.txt \
	    -g /scratch/mahone1/GWAS_QC_data/A4/Genotyped/Cleaned/part2_allraces/A4_genotyped_geno05_maf01_mind01_norelated_sex_nomismatchedsex_keep_autosomes_hwe6_nopal_noprobSNPs_chr_bplifted_nosamepos \
	    -s /ref/topmed_snp_names -x"
```

Post-imputation QC main output files:
- Main log file
- PCs and PC plots
- Cleaned, imputed+genotyped plink files

After the post-imputation QC runs, check the principal components for outliers, removing any if necessary and recalculating PCs.

To recalculate PCs, run a command similar to the one below. Adding the -n flag will suppress the creation of a file for default inclusions (ie only whites who fall within 5 SD from the mean).
```
singularity exec --contain --bind /path/to/genotype/data/:/inputs/ \
  /data/h_vmac/GWAS_QC/singularity/CNT_genomic_processing_v2.simg \
  /bin/bash -c  "cd /scripts/GWAS_QC/ ; sh calc_plot_PCs.sh \
  -i /inputs/Imputed/Raw/COHORT_imputed_NHW_ids_sex_maf01_hwe6_ids"
```

If you get an error because there is incomplete overlap between the genotyped and imputed files, it is likely there was some issue with converting the imputed ids back to FID and IID. To resolve this issue, ensure that the supplied race/sex file has every sample in the imputed files.
