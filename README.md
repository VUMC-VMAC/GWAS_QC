# GWAS QC and TOPMed Imputation Pipeline Documentation

This README describes the CNT methods for QC and imputation of genotype data. Specifically, it describes how to run the main autosomal pipeline, how to QC the X chromosome, how to perform multi-set merges, as well as how to run PC calculations on individual sets. 

## Autosomal Pipeline

Here are the general steps in this pipeline:
1. Pre-Imputation QC
	- Variant filtering excluding for missingness (>5%) and minor allele frequency (<0.01)
	- Sample missingness filtering (1%)
	- Relatedness check (removes 1 of each pair of second degree relatives, 0.9> pi-hat >0.25, and both of each pair with a pi-hat >=0.9)
	- Sex check (removes individuals with discrepancies)
	- Heterozygosity check
	- Hardy-Weinberg filter (p<1e-6)
	- Removing palindromic variants
	- Lifts to build 38 if necessary
	- Removing one of each duplicate variant and multi-allelic variants
	- Compare the files with the reference panel
	- Prepares files for uploading to the imputation server
2. TOPMed Imputation
3. Post-Imputation QC
	- Part 1
		- Unzip if needed
		- Variant filtering excluding for:
			- imputation quality (<0.8)
			- multi-allelic variants
		- Update sample IDs and sex
		- Infer ancestry using pre-calculated weights
	- Part 2
		- Subset to the supplied subset and filter variants for missingness (5%) in genotype and imputed variants
		- Merging in genotyped variants with the imputed data
		- Hardy-Weinberg equilibrium (p<1e-6) and minor allele frequency (1%)
		- Heterozygosity check
		- PC calculation

Each step (except the TOPMed imputation) is run as a bash script in a singularity. The usage message for each of the scripts can be shown by running the script with the -h flag.

In general, it is best to run the QC scripts as slurm jobs. There is an example slurm script for running pre-imputation QC in this repository: run_gwas_qc_preimputation.slurm However, if it is desired to run the scripts locally, most scripts have the -m flag which can be used to limit the memory allotments of plink commands making it possible to run locally. 

Current version of the singularity: CNT_genomic_processing_v3.10.simg

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

### Pre-Imputation QC

Pre-imputation QC can be run in a single script which runs all the filters and finishes by formatting the genotype data for upload to the TOPMed imputation server. 

Files needed for the Pre-Imputation QC script:
- Raw genotype files
- Singularity container 
- File with FID, IID, and sex for each set
  - Sex should be coded as 1=male, 2=female, 0=unknown.
- Reference panel files

To run the preimputation script, run a singularity command similar to the one below, binding the paths to the necessary input files. 
```
singularity exec --contain \
	--bind /nobackup/h_vmac/mahone1/GWAS_QC/:/input \
	--bind /data/h_vmac/:/ref/ /data/h_vmac/GWAS_QC/singularity/ \
	CNT_genomic_processing_v3.10.simg /bin/bash -c "cd /scripts/GWAS_QC/ ; \
		./gwas_qc_preimputation.sh \
			-i /input/VMAP/Genotyped/Raw/VMAP2_mapids \ 
			-o /input/VMAP/Genotyped/Cleaned/VMAP2_mapids \
			-f /input/VMAP/Genotyped/Raw/VMAP2_sex.txt \
			-R /ref/GWAS_QC/topmed/topmed_imputed_ref_snps \
			-b b37 \
			-m 15000 |& tee /input/VMAP/Genotyped/Cleaned/VMAP2_genotyping_QC.log"
```
Explanation of flags:
- -i : input genotypes stem 
- -o : output stem
- -f : file with sex information
- -R : reference file stem
- Optional flags:
	- -b : input build (default = b37)
	- -m : value for plink --memory flag to limit resource usage
	- -t : include related but not identical samples (see "GWAS QC including related individuals" section for more details)
	- -n : skip exclusion of variants not matching the reference panel
	- -c : skip the clean-up step at the end 

Pre-Imputation QC output files:
- Main log file
- Heterozygosity plot: file ending in *_pruned_hetcheck.png
- .vcf.gz files for each chromosome

After this script completes, before uploading to the imputation server, check the number of samples and variants filtered at each step as well as the heterozygosity plot to ensure quality of the data. If the numbers of samples/variants removed at each step look acceptable, then upload the resulting vcf files to the imputation server. 

### TOPMed Imputation

Upload the .vcf.gz files from the preimputation script to the TOPMed Imputation Server (https://imputation.biodatacatalyst.nhlbi.nih.gov/#!). On the imputation upload, select build 38 (since the results will have been lifted to b38 by the preimputation script), phasing using Eagle, and the TOPMed panel as the population (which will check for significant MAF differences). Optionally, you can select an Rsq filter which will automatically drop variants below the selected Rsq threshold, making the imputation results files smaller. 

Once imputation is complete, download the results and save the decryption password (received in an email from the imputation server) to a file called pass.txt in the folder with the imputation results. Alternatively, the imputation results can be unzipped before running the post-imputation QC script as long as the -z flag is excluded from the command.

### Post-Imputation QC

Running post-imputation QC will take place in two parts. The first part involves doing initial filtering of the imputed variants and then estimating genetic ancestry using SNPWeights. SNPWeights allows for inferring ancestry and subsetting people into groups with greater genetic similarity rather than simply subsetting by self-report race. The second step involves subsetting samples into the ancestry groups defined by SNPWeights, merging imputed and genotyped variants together, and then performing final variant (MAF and HWE) and sample (heterozygosity and PCs) filtering. 

Files needed for the Post-imputation QC:
- Singularity container 
- Imputation results and password to unzip them
- File with sex for the dataset (used in the pre-imputation script)
  - Note that this *must* include every sample that is in the imputation files to be able to update FID/IID and sex correctly as well as to merge in the genotype data correctly.
- Files to convert variant ids back to rs number
- Cleaned genotype files to merge back in
  - This will be the last file from pre-imputation QC, before prep for the imputation server (ending in *_chr_bplifted_nosamepos or *_chr_bplifted, depending on whether there are multi-allelic or duplicated variants to remove)

If you get an error because there is incomplete overlap between the genotyped and imputed files, it is likely there was some issue with converting the imputed ids back to FID and IID. To resolve this issue, ensure that the supplied sex file has every sample in the imputed files.

#### Post-Imputation QC Part 1

The first part of post-imputation QC performs initial quality filtering and assigns ancestry categories based on SNPWeights. 

Files needed for part 1 of post-imputation QC:
- Singularity container 
- Imputation results and password to unzip them
- File with sex for the dataset (used in the pre-imputation part 1 step)
  - Note that this *must* include every sample that is in the imputation files to be able to update FID/IID and sex correctly as well as to merge in the genotype data correctly.
- Files to convert variant ids back to rs number
- File with pre-calculated SNP weights from which to infer ancestry (current version: /data/h_vmac/GWAS_QC/1000G_data/1000G_snpwt; note that there must be an additional file with the same stem ending in *_snp.txt that lists the variants in the weights file). 

To run the first step in the SNPWeights post-imputation step, run something like:
```
singularity exec --contain \
	--bind /nobackup/h_vmac/mahone1/GWAS_QC/VMAP/:/input/ \
	--bind /data/h_vmac/GWAS_QC/:/ref/ \
	/data/h_vmac/GWAS_QC/singularity/CNT_genomic_processing_v3.10.simg /bin/bash -c "cd /scripts/GWAS_QC/ ; \
		./gwas_qc_postimputation_part1.sh \
			-i /input/Imputed/Raw/VMAP2/ \
			-o /input/Imputed/Cleaned/VMAP2/VMAP2_imputed \
			-f /input/Genotyped/Raw/VMAP2_sex.txt \
			-w /ref/1000G_data/1000G_snpwt \
			-s /ref/topmed/topmed_snp_names \
			-m 15000 \
			-c |& tee /input/Imputed/Cleaned/VMAP2/VMAP2_postimputation_QC_part1.log"
```
Explanation of flags:
- -i : input folder with raw imputation results
- -o : output stem 
- -f : file with sex information
- -w : reference file with pre-calculated SNP weights
- -s : reference file with rsIDs
- Optional flags:
	- -z : unzip the files from the imputation server
	- -x : skip the first post-imputation filters (helpful if these have already been done)
	- -m : value for plink --memory flag to limit resource usage
	- -c : skip cleanup at the end


The outputs from the first step will be a plink file set with imputed variants for all samples (*_IDs_sex.{bed,bim,fam}), a set of files for each ancestry subset in which there are >10 samples (*keep.txt), and a *.png plot of the ancestry weights. Before moving to the next step, look at the plot and scan the numbers from this step to ensure nothing went awry. If everything looks fine, proceed in running part 2 post-imputation on each of the subsets with sufficient samples (>10).

#### Post-Imputation QC Part 2

The second postimputation script must be run in each ancestry subset separately. It will perform the rest of the post-imputation filtering steps, other than PC outlier removal which must be done manually. 

Files needed for part 2 of post-imputation QC:
- Singularity container 
- Plink fileset from part 1 (*_IDs_sex.{bed,bim,fam})
- List of samples to keep for this subset from part 1 (*keep.txt)
- Final pre-imputation genotyped plink fileset
- File with sex for the dataset (used in pre-imputation)
  - Note that this *must* include every sample that is in the imputation files to be able to update FID/IID and sex correctly as well as to merge in the genotype data correctly.
- Files to convert variant ids back to rs number (for genotype files)

To run the second post-imputation step, run something like the following command:
```
singularity exec --contain \
	--bind /nobackup/h_vmac/mahone1/GWAS_QC/VMAP/:/input/ \
	--bind /data/h_vmac/GWAS_QC/:/ref/ \
	/data/h_vmac/GWAS_QC/singularity/CNT_genomic_processing_v3.10.simg /bin/bash -c "cd /scripts/GWAS_QC/ ; \
		./gwas_qc_postimputation_part2.sh \
			-i /input/Imputed/Cleaned/VMAP2/VMAP2_imputed_IDs_sex \
			-g /input/Genotyped/Cleaned/VMAP2_mapids_forimputation-updated \
			-r /input/Imputed/Cleaned/VMAP2/VMAP2_imputed_IDs_sex_overlap_acgt_pruned_AFR_keep.txt \
			-l AFR -m 15000 -c -p \
			-s /ref/topmed/topmed_snp_names |& tee /input/Imputed/Cleaned/VMAP2/VMAP2_postimputation_QC_AFR_part2.log"
```
Explanation of flags:
- -i : file stem of the plink set output from part 1 post-imputation
- -g : file stem for preimputation genotype files to merge back in
- -r : name of the file to subset to this ancestry group (the *keep.txt file from part 1 post-imputation)
- -l : label for this ancestry group (EUR, AFR, AMR2admx, AMR3admx, AMR)
- -s : reference file with rsIDs
- Optional flags:
	- -p : skip PC calculation (only for sets that will be merged or that have related individuals and should be processed with PC-AiR)
	- -m : value for plink --memory flag to limit resource usage
	- -c : skip cleanup at the end

After this script is complete, check the PC plots from the final step to check for outliers and make decisions on samples which need to be removed. If there are samples to be removed, drop those samples from the genotype files and recalculate PCs. (See the section "Additional notes on PC calculation" below for instruction on doing so in the singularity.) If these resulting PCs look acceptable, then the QC of this set is complete. If additional outliers must be removed, do so and recalculate. 

If the set you are QC'ing is one of several for a given cohort, PC calculation can be skipped at the end of this step (with the -p) flag. In that case, PCs should be calculated in the whole merged set. (See "Multi-set Merges" below for instruction on how to merge multiple genotype sets together.)

#### Appendix: Non-SNPWeights Post-Imputation QC

Note: This section is retained for backward compatibility, but should not be run. 

To run the post-imputation QC without SNPWeights, run a command similar to the one below.
```
singularity exec --contain \
	    --bind /nobackup/h_vmac/mahone1/:/input \
	    --bind /data/h_vmac/GWAS_QC/topmed/:/ref/ \
	    /data/h_vmac/GWAS_QC/singularity/CNT_genomic_processing_v3.10.simg /bin/bash -c  "cd /scripts/GWAS_QC/ ; \
	    ./gwas_qc_postimputation_old.sh -i /input/GWAS_QC_data/A4/Imputed/Raw/ \
			-o /input/GWAS_QC_data/A4/Imputed/Cleaned/AllRaces/A4_AllRaces_imputed  \
			-r /input/GWAS_QC_data/A4/A4_race.txt \
			-f /input/GWAS_QC_data/A4/A4_sex.txt \
			-g /input/GWAS_QC_data/A4/Genotyped/Cleaned/part2_allraces/A4_genotyped_geno05_maf01_mind01_norelated_sex_nomismatchedsex_keep_autosomes_hwe6_nopal_noprobSNPs_chr_bplifted_nosamepos \
			-s /ref/topmed_snp_names"
```

Post-imputation QC main output files:
- Main log file
- PCs and PC plots
- Cleaned, imputed+genotyped plink files

After the post-imputation QC runs, check the principal components for outliers, removing any if necessary and recalculating PCs.

## X Chromosome GWAS QC

The X chromosome pipeline is similar to the autosomal pipeline with a few key differences. A few of the variant filters occur in a sex-specific manner. Additionally, some of the sample filters (like heterozygosity, ancestry subsetting, and PC outlier removal) occur based on autosomal thresholds. Thus, X chromosome QC must be run separately from autosomal QC. Below are instructions for running the X chromosome QC pipeline. 

Here are the general steps in this pipeline:
1. Pre-Imputation QC
	- Variant filtering excluding for missingness (>5%) and minor allele frequency (<0.01)
	- Sample missingness filtering (1%)
	- Relatedness check (removes 1 of each pair of second degree relatives, 0.9> pi-hat >0.25, and both of each pair with a pi-hat >=0.9)
	- Sex check (removes individuals with discrepancies)
	- Subset to X chr variants
	- Remove variants with high differential missingness between males and females (p<1e-7)
	- Hardy-Weinberg filter (p<1e-6 based on females)
	- Removing palindromic variants
	- Lifts to build 38 if necessary
	- Removing one of each duplicate variant and multi-allelic variants
	- Compare the files with the reference panel
	- Prepares files for uploading to the imputation server
2. TOPMed Imputation
3. Post-Imputation QC
	- Unzip if needed
	- Variant filtering excluding for:
		- imputation quality (<0.8)
		- multi-allelic variants
	- Update sample IDs and sex
	- Filters by ancestry group
		- Subset to the supplied subset and filter variants for missingness (5%) in genotype and imputed variants
		- Merging in genotyped variants with the imputed data
		- Hardy-Weinberg equilibrium (p<1e-6 based on females) and minor allele frequency (1%)
		- Set heterozygous haploid variants missing
		- Remove PC outliers based on the final autosomal plink file

### Pre-imputation QC (X Chromosome)

Files needed for the Pre-Imputation QC script:
- Raw genotype files
- Singularity container 
- File with FID, IID, and sex for each set
  - Sex should be coded as 1=male, 2=female, 0=unknown.
- Reference panel files

To run the pre-imputation script, run a singularity command similar to the one below, binding the paths to the necessary input files. 
```
singularity exec --contain \
	--bind /nobackup/h_vmac/mahone1/GWAS_QC/VMAP/:/input/ \
	--bind /data/h_vmac/GWAS_QC/:/data/ \
	/data/h_vmac/GWAS_QC/singularity/CNT_genomic_processing_v3.10.simg \
	/bin/bash -c "cd /scripts/GWAS_QC/ ; \
		./gwas_qc_preimputation_chrX.sh \
		-o /input/Genotyped/Cleaned/Xchr/VMAP2_X \
		-i /input/Genotyped/Raw/VMAP2_mapids \
		-f /input/Genotyped/Raw/VMAP2_sex.txt \
		-m 10000 \
		-R /data/topmed/X.PASS.Variants.TOPMed_freeze5_hg38_dbSNP.tab.gz \
		-b b37 |& tee /input/Genotyped/Cleaned/Xchr/VMAP2_Xchr_preimputation_qc.log"
```
Explanation of flags:
- -i : input genotypes stem 
- -o : output stem
- -f : file with sex information
- -R : reference file name
- Optional flags:
	- -b : input build (default = b37)
	- -m : value for plink --memory flag to limit resource usage
	- -t : include related but not identical samples (see "GWAS QC including related individuals" section for more details)
	- -n : skip exclusion of variants not matching the reference panel
	- -c : skip the clean-up step at the end 

Pre-Imputation QC output files:
- Main log file
- .vcf.gz files for chromosome X

After this script completes, before uploading to the imputation server, check the number of samples and variants filtered at each step. If the numbers of samples/variants removed at each step look acceptable, then upload the resulting vcf files to the imputation server. 

Note that if there were samples removed for heterozygosity in the autosomal pipeline they must be removed before uploading to the imputation server. However, this is not likely to occur. 

### TOPMed Imputation (X Chromosome)

Upload in the same manner as autosomal variant files. 

### Post-imptuation QC (X Chromosome)

Unlike autosomal QC, post-imputation QC for the X chromosome is run in one script, using outputs from the autosomal QC. Thus, it must be run after the autosomal QC is complete. Also unlike autosomal post-imputation QC, this script can process all the ancestry subgroups in one run. 

Files needed for post-imputation QC:
- Singularity container 
- Imputation results and password to unzip them
- File with sex for the dataset (used in the pre-imputation part 1 step)
  - Note that this *must* include every sample that is in the imputation files to be able to update FID/IID and sex correctly as well as to merge in the genotype data correctly.
- Files to convert variant ids back to rs number
- Subset-specific files (comma-separated, in the same order, with no spaces):
	- Lists of samples for the initial sample subsetting
	- Labels for each subset
	- File stem to the final autosomal plink file set (with PC outliers removed)

In order to run X chromosome post-imputation QC, run a command similar to the following:
```
singularity exec --contain \
	--bind /data/h_vmac/GWAS_QC/topmed/:/ref/ \
	--bind /nobackup/h_vmac/mahone1/GWAS_QC/VMAP/:/input/ \
	/data/h_vmac/GWAS_QC/singularity/CNT_genomic_processing_v3.10.simg \
	/bin/bash -c "cd /scripts/GWAS_QC/ ; \
		./gwas_qc_postimputation_chrX.sh \
		-o /input/Imputed/Cleaned/VMAP2/VMAP2_Xchr \
		-i /input/Imputed/Raw/VMAP2/Xchr/ \
		-f /input/Genotyped/Raw/VMAP2_sex.txt \
		-s /ref/topmed_snp_names \
		-g /input/Genotyped/Cleaned/Xchr/VMAP2_X-updated-chr23_TOPMED_varID \
		-c /input/Imputed/Cleaned/VMAP2/VMAP2_imputed_IDs_sex_overlap_acgt_pruned_test_EUR_keep.txt,/input/Imputed/Cleaned/VMAP2/VMAP2_imputed_IDs_sex_overlap_acgt_pruned_test_AFR_keep.txt,/input/Imputed/Cleaned/VMAP2/VMAP2_imputed_IDs_sex_overlap_acgt_pruned_test_AMR2admx_keep.txt \
		-p /input/Imputed/Cleaned/VMAP2/VMAP2_imputed_IDs_sex_nogeno_merged_EUR_hwe6_maf01_geno01,/input/Imputed/Cleaned/VMAP2/VMAP2_imputed_IDs_sex_nogeno_merged_AFR_hwe6_maf01_geno01,/input/Imputed/Cleaned/VMAP2/VMAP2_imputed_IDs_sex_nogeno_merged_AMR2admx_hwe6_maf01_geno01 \
		-l EUR,AFR,AMR2admx -d \
		-m 18000 |& tee /input/Imputed/Cleaned/VMAP2/Xchr/VMAP2_Xchr_postimputation_qc.log"
```
Explanation of flags:
- -i : input folder with raw imputation results
- -o : output stem 
- -f : file with sex information
- -g : file stem for preimputation genotype files to merge back in
- -s : reference file with rsIDs
- Subset specific arguments (comma-separated and in order if there are multiple)
	- -c : names of subset files for each ancestry group (the *keep.txt files from part 1 of autosomal post-imputation QC)
	- -p : names of the final autosomal files for each ancestry group
	- -l : labels for ancestry groups (EUR, AFR, AMR2admx, AMR3admx, AMR)
- Optional flags:
	- -z : unzip the files from the imputation server
	- -x : skip the first post-imputation filters (helpful if these have already been done)
	- -m : value for plink --memory flag to limit resource usage
	- -c : skip cleanup at the end

Once this runs, check the outputs of this script, confirming that the variants and samples removed at each step are reasonable. If so, X chromosome GWAS QC for this set is complete. 

## GWAS QC including related individuals

Some cohorts specifically recruited related individuals or have a high level of known relatedness among samples. For these, we will perform GWAS QC in a manner that allows inclusion of related individuals. While these individuals should not generally be included in most analyses together, it is helpful to have the genotype data for the highest number of samples possible so that if an analysis is being run on a phenotype that was collected in only some, the samples which have data are able to be retained and then filtered for relatedness. This will hypothetically allow for maximizing sample sizes in analyses with these cohorts. 

The process for QC'ing genotypes including related individuals is very similar to the standard pipeline with two significant modifications:
- Related individuals are not dropped from the dataset prior to imputation. Only identical samples are removed. 
- After subsetting to ancestry groups post-imputation, PC calculation must be performed using PC-AiR (PC Analysis in Related samples). 

To run pre-imputation QC including related samples for autosomes and X chromosome, run the same pre-imputation scripts as above simply adding the -t flag to exclude only samples which are identical but not those which are related. Run imputation and post-imputation QC as above, skipping PC calculation in post-imputation part 2 for autosomes (add the -p flag). Then, calculate PCs using PC-AiR. There is an example script in this repository (run_PCAiR.R), but note that it cannot be run in the singularity. Since this is an iterative process and requires evaluation of the plots, it may be easiest to run this process in RStudio, but it can be run on accre as well. 


## Multi-set Merges

Some cohorts have genotype data from multiple genotyping runs. This is often the case if the genotype data were collected at multiple sites (as in NACC) or at different times (as in VMAP). Typically, the genotype data must be processed by batch all the way through post-imputation QC and then merged at the end of the QC pipeline. The only step that can be skipped for these sets is PC outlier removal at the end of post-imputation QC. This should be performed in the merged genotype data for the whole cohort. However, when merging sets that have small sample sizes, there is a slightly different method in order to retain as many variants as possible. Below are the steps to perform a multi-set merge followed by the alternative method for sets with a small sample size. These remain the same no matter how many sets are being merged. These are also advisable if genotype data is being merged between cohorts for analysis. 

The basic process for merging multiple sets of GWAS data for a given cohort is the following:
1. Begin with the the sets after genotype data has been merged in and variant filters have been applied but before calculating PCs. Check for sample overlap between sets.
   	- If there is overlap, check for concordance of data between the sets for those overlapping samples to ensure that they are genetically identical.
   	- In general, keep samples in the set which has the most variants in the final preimputation file.
   	- There will likely be variants that cause the first attempt at a merge fail which will be saved to a *.missnp file by plink. Remove those variants in order to complete the merge. 
2. Check A1/MAF between sets and remove variants with differences > 10% (using the pre_merge_A1_MAF_check.R script).
3. Merge the sets.
4. Check relatedness in the merged set and remove individuals with relatedness >0.25.
5. Calculate PCs and remove outliers.
	- It is sometimes helpful to color PC plots on set to ensure that the sets are not clustering together. 

Here is the modified process for merging sets with small sample sizes:
1. Begin with the the sets after genotype data has been merged in but before variant filters have been applied or PCs have been calculated. Check for sample overlap between sets.
   	- If there is overlap, check for concordance of data between the sets for those overlapping samples to ensure that they are genetically identical.
   	- In general, keep samples in the set which has the most variants in the final preimputation file.
2. Subset to overlapping variants. 
3. Merge the sets.
4. Filter to overlapping variants in the merged set using a missingness filter.
5. Filter variants for HWE p<0.000001 and MAF<0.01. 
6. Check relatedness in the merged set and remove individuals with relatedness >0.25.
7. Calculate PCs and remove outliers.
	- It is sometimes helpful to color PC plots on set to ensure that the sets are not clustering together. 

## Additional notes on PC calculation

To recalculate PCs using smartpca from eigensoft, run a command similar to the one below. Adding the -n flag will suppress the creation of a file for default inclusions (ie only whites who fall within 5 SD from the mean).
```
singularity exec --contain \
	--bind /path/to/genotype/data/:/inputs/ \
  	/data/h_vmac/GWAS_QC/singularity/CNT_genomic_processing_v3.10.simg \
  	/bin/bash -c  "cd /scripts/GWAS_QC/ ; ./calc_plot_PCs.sh \
  		-i /inputs/Imputed/Raw/COHORT_imputed_NHW_ids_sex_maf01_hwe6_ids"
```

An additional method to calculate PCs is also available using SNPRelate in R. This method is significantly faster than the smartpca method. To calculate PCs using this method, run a command similar to this:
```
singularity exec --contain --bind /path/to/genotype/data/:/inputs/ \
  	/data/h_vmac/GWAS_QC/singularity/CNT_genomic_processing_v3.10.simg \
	/bin/bash -c  "cd /scripts/GWAS_QC/ ; ./calc_plot_PCs_snprelate.sh \
  		-i /inputs/Imputed/Raw/COHORT_imputed_NHW_ids_sex_maf01_hwe6_ids"
```
