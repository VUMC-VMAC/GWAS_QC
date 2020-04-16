#!/bin/sh
filestem=$1 #give whole path and file stem, ie /nfs/DATA/A4/GWAS/Imputed/QC/A4_maf....

#This script will assume that the file path is also where the PCs should be calculated. Give the file stem and the file name. PC file names will be ${filename}_pccalc.

set -e

path=$(  echo ${filestem%/*} )
filename=$( echo ${filestem##*/} )

##module purge
module purge
module load PLINK/1.9b_5.2 GCC/5.4.0-2.26 GSL/2.1 OpenBLAS/0.2.18-LAPACK-3.6.1

printf "Calculating PCs for the ${filename} dataset in ${path}\n\n"

#move to right directory
cd $path

#prune and set phenotypes to 1
plink --bfile $filestem --indep-pairwise 200 100 0.2 --allow-no-sex --out  ${filename}_forprune
plink --bfile $filestem --output-missing-phenotype 1 --extract  ${filename}_forprune.prune.in --make-bed --out ${filename}_pruned

echo "Pruning complete. Running smartpca for raw ${cohort}..."

printf "genotypename: ${filename}_pruned.bed
snpname: ${filename}_pruned.bim
indivname: ${filename}_pruned.fam
evecoutname: ${filename}.pca.evec
evaloutname: ${filename}.eigenvalues
altnormstyle: NO
numoutevec: 10
numoutlieriter: 0
numoutlierevec: 10
outliersigmathresh: 6
qtmode: 0" > ${filename}_pccalc.par

smartpca -p ${filename}_pccalc.par > ${filename}_pccalc.log

echo "smartpca for the ${filename} dataset complete!"

rm *forprune*
