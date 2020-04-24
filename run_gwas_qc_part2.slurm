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

module purge 
module load PLINK/1.9b_5.2 GCC/5.4.0-2.26 GSL/2.1 OpenBLAS/0.2.18-LAPACK-3.6.1 OpenMPI/1.10.3 R/3.3.3 tabix/0.2.6 VCFtools/0.1.14 Pango/.1.40.3

#BCFtools
PATH=$PATH:/home/clarklc/Programs/bcftools-1.3.1/bin
export PATH

cd /data/h_vmac/GWAS_QC
sh gwas_qc_part2.sh STEM INPUT_FILESET REF_FILE noexclude
