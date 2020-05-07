# Example QC


To run the part 1 script: 
```
singularity exec --contain --bind /nfs/DATA/BIOCARD/GWAS/Genotyped/:/inputs/ \
	    --bind /data/h_vmac/GWAS_QC/1000G_data/:/inputs2/ \
	    --bind /data/h_vmac/GWAS_QC/:/scripts/ \
	    /data/h_vmac/GWAS_QC/singularity/gwas_qc_singularity_v1.3.simg \
	    /bin/bash -c "cd /scripts/ ; sh /scripts/gwas_qc_part1.sh \
	    	      -i /inputs/Cleaned/QC_topmed/BIOCARD_genotyped_ids \
		      	 -o /inputs/Cleaned/QC_topmed_recap/BIOCARD_genotyped \
			    -r /inputs/Cleaned/QC_topmed/BIOCARD_race_sex.txt \
			       -G /inputs2/all_1000G_maf001 |& tee /inputs/Cleaned/QC_topmed_recap/BIOCARD_qc_p1.log"
```

Then, make your decision regarding the PC plots. In this case, decided to go with the default decision.
```
plink --bfile  /nfs/DATA/BIOCARD/GWAS/Genotyped/Cleaned/QC_topmed_recap/BIOCARD_genotyped_geno05_maf01_mind01_norelated_nomismatchedsex_keep_autosomes --keep  /nfs/DATA/BIOCARD/GWAS/Genotyped/Cleaned/QC_topmed_recap/BIOCARD_genotyped_1000G_geno05_pruned_no1000G_nooutliers.txt --make-bed --out  /scratch/mahone1/GWAS_QC_data/BIOCARD_recap/BIOCARD_genotyped_geno05_maf01_mind01_norelated_nomismatchedsex_keep_autosomes_NHW
```

To run the part 2 script, create a slurm script like this:
```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=28G
#SBATCH --output=/scratch/mahone1/GWAS_QC/BIOCARD_recap/BIOCARD_genotyped_NHW_cleaned_gwas_qc_part2_slurm.log
#SBATCH --mail-user=emily.mahoney@vumc.org
#SBATCH --mail-type=ALL
#SBATCH --job-name=BIOCARD_genotyped_NHW_cleaned_gwas_qc_part2
#SBATCH --time=06:00:00
#SBATCH --account=h_vmac

singularity exec --contain --bind /scratch/mahone1/GWAS_QC_data/BIOCARD_recap/:/inputs/ \
    --bind /data/h_vmac/GWAS_QC/topmed/:/ref/ \
    --bind /data/h_vmac/GWAS_QC/:/scripts/ \
    /data/h_vmac/GWAS_QC/singularity/gwas_qc_singularity_v1.3.simg \
/bin/bash -c "cd /scripts/ ; sh gwas_qc_part2.sh \
-i /inputs/BIOCARD_genotyped_geno05_maf01_mind01_norelated_nomismatchedsex_keep_autosomes_NHW \
-o /inputs/BIOCARD_genotyped_NHW_cleaned \
-R /ref/topmed_imputed_ref_snps -n |& tee /inputs/VMAP_qc_p2.log"
```

Submit that using slurm:
```
sbatch /scratch/mahone1/GWAS_QC_data/BIOCARD_recap/BIOCARD_NHW_run_gwas_qc_part2.slurm
```

