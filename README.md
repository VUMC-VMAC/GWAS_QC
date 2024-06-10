# GWAS QC and TOPMed Imputation Pipeline Documentation

Here are the general steps in this pipeline:
1. Pre-Imputation QC
	- Variant filtering excluding for missingness (>5%) and minor allele frequency (<0.01)
	- Sample missingness filtering (5%)
	- Relatedness check (removes 1 of each pair of second degree relatives, 0.9> pi-hat >0.25, and both of each pair with a pi-hat >=0.9)
	- Sex check (removes individuals with discrepancies)
	- Heterozygosity check
	- Hardy-Weinberg filter (p<1e-6)
	- Removing palindromic, multi-allelic, and duplicate variants
	- Lifts to build 38 if necessary
	- Prepares files for uploading to the imputation server
2. TOPMed Imputation
3. Post-Imputation QC
	- Part 1
		- Variant filtering excluding for:
			- imputation quality (<0.8)
			- multi-allelic variants
		- Merging in genotyped variants with the imputed data
		- Infer ancestry using pre-calculated weights
	- Part 2
		- Hardy-Weinberg equilibrium (p<1e-6) and minor allele frequency (<0.01)
		- PC calculation

Each step (except the TOPMed imputation) is run as a shell script in a singularity. The usage message for each of the scripts can be shown by running the script with the -h flag.

In general, it is best to run the QC scripts as slurm jobs. There is an example slurm script for running pre-imputation QC in this repository: run_gwas_qc_preimputation.slurm

Software used for QC:
- plink, versions 1.9 and 2.0: https://www.cog-genomics.org/plink/
- R 3.6.3: https://www.r-project.org/
- bcftools: http://samtools.github.io/bcftools/bcftools.html
- htslib: http://www.htslib.org/
- Imputation checking script from the McCarthy group: https://www.well.ox.ac.uk/~wrayner/tools/
- Eigensoft: https://github.com/DReichLab/EIG
- SNPWeights: https://www.hsph.harvard.edu/alkes-price/software/

To be able to upload files for imputation to the TOPMed reference panel, please create an account on the TOPMed imputation server: https://imputation.biodatacatalyst.nhlbi.nih.gov/#!pages/home

Supplementary files not included in this github:
- 1000G genetic data: https://www.internationalgenome.org/data/
- TOPMed SNP information files
- Singularity to run the GWAS QC
- Weights to infer ancestry

## Pre-Imputation QC

Files needed for the Pre-Imputation QC script:
- Raw genotype files
- Singularity container (current version: CNT_genomic_processing_v3.6.simg)
- File with FID, IID, and sex for each set
  - Sex should be coded as 1=male, 2=female, 0=unknown.
- File with FID, IID, and race for each set (optional)
  - Non-Hispanic whites coded as "White" or "EUR" and none of the values can have spaces.
- 1000G files to use in calculating PCs (optional)
- Reference panel files

To run the part 1 script, run a singularity command similar to the one below, binding the paths to the necessary input files. This script will take between 40 minutes to an hour to run, depending on the size of your dataset. As in the example below, -i specifies the input genotypes stem, -o specifies the output stem, -r specifies the race file, -f specifies the sex file, -G specifies the 1000G data, the -R specifies the reference file stem, and -b specifies the input build.
```
singularity exec --containall --bind /nobackup/h_vmac/mahone1/GWAS_QC/:/scratch \
	    --bind  /data/h_vmac/:/data/ /data/h_vmac/GWAS_QC/singularity/CNT_genomic_processing_v3.5.simg /bin/bash -c "cd /data/mahone1/GWAS_QC/ ; \
	    sh gwas_qc_preimputation.sh -i /scratch/VMAP/Genotyped/Raw/VMAP2_mapids \
	       				-o /scratch/VMAP/Genotyped/Cleaned/VMAP2_mapids \
					-f /scratch/VMAP/Genotyped/Raw/VMAP2_sex.txt \
					-G /data/GWAS_QC/1000G_data/1000G_final_b37 \
					-R /data/GWAS_QC/topmed/topmed_imputed_ref_snps -b b37 \
					|& tee /scratch/VMAP/Genotyped/Cleaned/VMAP2_genotyping_QC_preimputation.log"
```
The -n flag would specify to not exclude variants based on the reference panel. The -c flag would specify to skip the clean-up step at the end of the script which will remove bed and vcf files prior to the final ones. 

Pre-Imputation QC output files:
- Main log file
- Heterozygosity plot: file ending in *_pruned_hetcheck.png
- .vcf.gz files for each chromosome

After this script completes, before uploading to the imputation server, check the number of samples and variants filtered at each step as well as the heterozygosity plot to ensure quality of the data. If the numbers of samples/variants removed at each step look acceptable, then upload the resulting vcf files to the imputation server. 

## TOPMed Imputation

Upload the .vcf.gz files from the preimputation script to the TOPMed Imputation Server (https://imputation.biodatacatalyst.nhlbi.nih.gov/#!). On the imputation upload, select build 38 (since the results will have been lifted to b38 by the preimputation script), phasing using Eagle, and the TOPMed panel as the population (which will check for significant MAF differences). Optionally, you can select an Rsq filter which will automatically drop variants below the selected Rsq threshold, making the imputation results files smaller. 

Once imputation is complete, download the results and save the decryption password (received in an email from the imputation server) to a file called pass.txt in the folder with the imputation results. Alternatively, the imputation results can be unzipped before running the post-imputation QC script as long as the -z flag is excluded from the command.

## Post-Imputation QC

Post-imputation QC begins with merging imputed and genotyped variants and then estimating genetic ancestry using SNPWeights. SNPWeights allows for inferring ancestry and subsetting people into groups with greater genetic similarity rather than simply self-report. Running post-imputation with SNPWeights will take place in two steps. The first step is to merge together the imputed variants and the genotyped variants and calculate ancestral cateogries using SNPWeights. Before the next step, samples will be split into ancestral categories. The second is to perform final filtering (MAF, HWE, heterozygosity, and final PCs) on each ancestral category.

Files needed for the Post-imputation QC:
- Singularity container (current version: CNT_genomic_processing_v3.6.simg)
- Imputation results and password to unzip them
- File with sex for the dataset (used in the pre-imputation part 1 step)
  - Note that this *must* include every sample that is in the imputation files to be able to update FID/IID and sex correctly as well as to merge in the genotype data correctly.
- Files to convert variant ids back to rs number
- Cleaned genotype files to merge back in
  - This will be the last file from pre-imputation QC, before prep for the imputation server (ending in *_chr_bplifted_nosamepos or *_chr_bplifted, depending on whether there are multi-allelic or duplicated variants to remove)

If you get an error because there is incomplete overlap between the genotyped and imputed files, it is likely there was some issue with converting the imputed ids back to FID and IID. To resolve this issue, ensure that the supplied sex file has every sample in the imputed files.

### Post-Imputation QC Part 1

Additional files required for the Post-imputation SNPWeights Part 1 script:
- File with pre-calculated SNP weights from which to infer ancestry (current version: /data/h_vmac/GWAS_QC/1000G_data/1000G_snpwt; note that there must be an additional file with the same stem ending in *_snp.txt that lists the variants in the weights file). 

The first step will run all the initial post-imputation steps up to inferring ancestry including:
- Unzipping results from the imputation server (if -z is supplied)
- Filter for R2>0.8
- Removing multi-allelic variants
- Merging back in genotyped variants
- Infer ancestry

To run the first step in the SNPWeights post-imputation step, run something like:
```
singularity exec --contain --bind /scratch:/scratch --bind /data/h_vmac/GWAS_QC/topmed/:/ref/ \
	    --bind /data/h_vmac/GWAS_QC/1000G_data/:/data/ \
	    /data/h_vmac/GWAS_QC/singularity/CNT_genomic_processing_v3.6.simg /bin/bash -c  "cd /scripts/GWAS_QC/ ; \
	    sh gwas_qc_postimputation_part1.sh -i /scratch/mahone1/GWAS_QC_data/A4/Imputed/Raw/ \
	    -o /scratch/mahone1/GWAS_QC_data/A4/Imputed/Cleaned/AllRaces/A4_AllRaces_imputed  \
	    -r /scratch/mahone1/GWAS_QC_data/A4/A4_race_sex.txt \
	    -g /scratch/mahone1/GWAS_QC_data/A4/Genotyped/Cleaned/part2_allraces/A4_genotyped_geno05_maf01_mind01_norelated_sex_nomismatchedsex_keep_autosomes_hwe6_nopal_noprobSNPs_chr_bplifted_nosamepos \
	    -w /data/1000G_snpwt
	    -s /ref/topmed_snp_names"
```

The outputs from the first step will be a file (ending in *.out) which contains predicted ancestry and a pdf with plots of the predicted PCs. Before moving to the next step, subset the data into whatever ancestral categories are present and run the second step on each set. 

### Post-Imputation QC Part 2

The second script will run the post-imputation steps and should be run in the dataset after separating the data into the genetically similar groups. Specifically, the following steps will be done:
- Filter out variants with HWE p<0.000001 and MAF<0.01
- Check heterozygosity and remove outliers if present
- Calculate PCs

To run the second step in the SNPWeights post-imputation step, run something like:
```
singularity exec --contain --bind /scratch:/scratch --bind /data/h_vmac/GWAS_QC/topmed/:/ref/ \
	    /data/h_vmac/GWAS_QC/singularity/CNT_genomic_processing_v3.6.simg /bin/bash -c  "cd /scripts/GWAS_QC/ ; \
	    sh gwas_qc_postimputation_part2.sh \
	    -o /scratch/mahone1/GWAS_QC_data/A4/Imputed/Cleaned/AllRaces/A4_AllRaces_imputed  \
	    -r /scratch/mahone1/GWAS_QC_data/A4/A4_race.txt"
```


### Appendix: Non-SNPWeights Post-Imputation QC

To run the post-imputation QC without SNPWeights, run a command similar to the one below.
```
singularity exec --contain --bind /scratch:/scratch --bind /data/h_vmac/GWAS_QC/topmed/:/ref/ \
	    /data/h_vmac/GWAS_QC/singularity/CNT_genomic_processing_v3.6.simg /bin/bash -c  "cd /scripts/GWAS_QC/ ; \
	    sh gwas_qc_postimputation_old.sh -i /scratch/mahone1/GWAS_QC_data/A4/Imputed/Raw/ \
	    -o /scratch/mahone1/GWAS_QC_data/A4/Imputed/Cleaned/AllRaces/A4_AllRaces_imputed  \
	    -r /scratch/mahone1/GWAS_QC_data/A4/A4_race.txt \
	    -f /scratch/mahone1/GWAS_QC_data/A4/A4_sex.txt \
	    -g /scratch/mahone1/GWAS_QC_data/A4/Genotyped/Cleaned/part2_allraces/A4_genotyped_geno05_maf01_mind01_norelated_sex_nomismatchedsex_keep_autosomes_hwe6_nopal_noprobSNPs_chr_bplifted_nosamepos \
	    -s /ref/topmed_snp_names"
```

Post-imputation QC main output files:
- Main log file
- PCs and PC plots
- Cleaned, imputed+genotyped plink files

After the post-imputation QC runs, check the principal components for outliers, removing any if necessary and recalculating PCs.


## Additional notes on PC calculation

To recalculate PCs using smartpca from eigensoft, run a command similar to the one below. Adding the -n flag will suppress the creation of a file for default inclusions (ie only whites who fall within 5 SD from the mean).
```
singularity exec --contain --bind /path/to/genotype/data/:/inputs/ \
  /data/h_vmac/GWAS_QC/singularity/CNT_genomic_processing_v3.6.simg \
  /bin/bash -c  "cd /scripts/GWAS_QC/ ; sh calc_plot_PCs.sh \
  -i /inputs/Imputed/Raw/COHORT_imputed_NHW_ids_sex_maf01_hwe6_ids"
```

An additional method to calculate PCs is also available using SNPRelate in R. This method is significantly faster than the smartpca method. To calculate PCs using this method, run a command similar to this:
```
singularity exec --contain --bind /path/to/genotype/data/:/inputs/ \
  /data/h_vmac/GWAS_QC/singularity/CNT_genomic_processing_v3.6.simg \
  /bin/bash -c  "cd /scripts/GWAS_QC/ ; sh calc_plot_PCs_snprelate.sh \
  -i /inputs/Imputed/Raw/COHORT_imputed_NHW_ids_sex_maf01_hwe6_ids"
```
