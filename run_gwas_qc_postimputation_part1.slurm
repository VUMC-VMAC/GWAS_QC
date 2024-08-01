#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=28G
#SBATCH --output=INPUT_PATH/STEM_gwas_qc_part2_slurm.log
#SBATCH --mail-user=EMAIL
#SBATCH --mail-type=ALL
#SBATCH --job-name=STEM_gwas_qc_part2
#SBATCH --time=06:00:00
#SBATCH --account=h_vmac

singularity exec --contain \
    --bind /nobackup/h_vmac/mahone1/GWAS_QC/VMAP/:/input/ \
    --bind /data/h_vmac/mahone1/GWAS_QC/:/data/ \
    --bind /data/h_vmac/GWAS_QC/:/ref/ \
    /data/h_vmac/GWAS_QC/singularity/CNT_genomic_processing_v3.10.simg /bin/bash -c "cd /scripts/GWAS_QC/ ; \
sh gwas_qc_postimputation_part1.sh \
-i /input/Imputed/Raw/VMAP2/ \
-o /input/Imputed/Cleaned/VMAP2/VMAP2_imputed \
-f /input/Genotyped/Raw/VMAP2_sex.txt \
-w /ref/1000G_data/1000G_snpwt \
-s /ref/topmed/topmed_snp_names \
-c |& tee /input/Imputed/Cleaned/VMAP2/VMAP2_postimputation_QC_part1.log"
