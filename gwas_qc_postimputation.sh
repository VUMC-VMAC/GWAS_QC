#/bin/sh
stem=$1 
imputation_results_folder=$2 
cleaned_imputed_folder=$3 
race_sex_file=$4 

#fail on error
set -e

#define usage
display_usage() {
    printf "GWAS QC post-imputation script
This script will unzip the imputation results (assuming password is saved in pass.txt in the same folder as the imputation results files and will perform standard post-imputation QC for our common variant pipeline. This includes filtering for R2, removing multi-allelic variants and filtering out variants for low MAF or HWE disequilibrium. Finally, PCs will be calculated on the final file-set.

Usage:
SCRIPTNAME.sh [stem] [imputation_results_folder] [cleaned_imputed_folder] [race_sex_file]

stem = the beginning part of all QC'ed files
imputation_results_folder = the folder to which the imputation results have been downloaded which also contains a file called pass.txt with the password to unzip the imputation results files
cleaned_imputed_folder = the folder to write the cleaned imputation files to
race_sex_file = a file with FID and IID (corresponding to the fam file), 1 column indicating both race and ethnicity for PC plots, and another indicating sex for the sex check (1 for males, 2 for females, 0 if unknown), with NO header. Non-hispanic whites need to be indicated with 'White.' No other values in the race column must be fixed. This will be used to update sex in the fam file. 
"
        }
#check if there are at least 3 arguments
if [  $# -lt 4 ]
        then
                display_usage
                exit 1
fi

#print out inputs
echo -e "GWAS QC Post-imputation Script\n"
echo "stem : "$stem
echo "imputation results folder : "$imputation_results_folder
echo "folder where cleaned, imputed files will be saved : "$cleaned_imputed_folder
echo "race/sex information file : "$race_sex_file

#check to make sure this is being run in the scripts folder (checking if necessary script is present)
if test ! -f get_related_ids.R ;
then
        printf "\nCurrently in $PWD, but must be in a folder with the necessary scripts to run the GWAS QC! Please move to that folder and run this script again.
Necessary scripts:
get_related_ids.R
check_id_length.R
plot_het_check_outliers.R
plot_PCs_generate_ids_to_keep.R
check_ambig_snps.R
HRC-1000G-check-bim-NoReadKey.pl\n"
        exit 1
fi

for i in $( ls ${imputation_results_folder}/*.zip );
do
    unzip -P $( cat ${imputation_results_folder}/pass.txt $i 
done


#load the most recent version of plink to have vcf conversion capabilities
module load PLINK/2.00-alpha1

for i in $(seq 1 22); do
    #restrict to variants with R2>=0.80
    plink2 --vcf ${imputation_results_folder}/chr${i}.dose.vcf.gz --exclude-if-info "R2<0.8" --make-bed --out ${cleaned_imputed_folder}/chr${i}_temp > /dev/null

    #get list of duplicates SNPs (they have more than 1 row in the bim file) and remove
    awk '{ print $2 }' ${cleaned_imputed_folder}/chr${i}_temp.bim | uniq -d > ${cleaned_imputed_folder}/chr${i}_temp.dups ;
    plink2 --bfile ${cleaned_imputed_folder}/chr${i}_temp --exclude ${cleaned_imputed_folder}/chr${i}_temp.dups --make-bed --out ${cleaned_imputed_folder}/chr${i}_temp_nodups > /dev/null
done

#print out numbers of variants
total_var=$(( $(grep "out of" ${cleaned_imputed_folder}/chr*_temp.log | awk 'BEGIN { ORS="+" } { print $4 }' | sed 's/\(.*\)+/\1 /' ) ))
afterR2=$( cat ${cleaned_imputed_folder}/chr*_temp_nodups.bim | wc -l )
nomulti=$( cat ${cleaned_imputed_folder}/chr*_temp.bim | wc -l )
echo -e "$total_var variants after imputation, $afterR2 variants with R2>0.8, and $nomulti biallelic variants\n"

#switch back to old plink
module unload PLINK/2.00-alpha1
module load PLINK/1.9b_5.2

#create merge file (emptying old version if present)
echo "" | tee ${cleaned_imputed_folder}/merge_list.txt
for i in $( seq 1 22 ); do echo "${cleaned_imputed_folder}/chr${i}_temp_nodups.bed  ${cleaned_imputed_folder}/chr${i}_temp_nodups.bim  ${cleaned_imputed_folder}/chr${i}_temp_nodups.fam" >> ${cleaned_imputed_folder}/merge_list.txt ; done

#merge individual chromosome files to be one large file
output=${cleaned_imputed_folder}/${cohort_stem}
plink --merge-list ${cleaned_imputed_folder}/merge_list.txt --make-bed --out $output

#update SNP names and add sex back into fam file
output_last=$output
output=${output}_names
plink --bfile ${output_last} --update-sex ${race_sex_file} 2 --update-name ${snp_names} --make-bed --out ${output}

# SNP filters
output_last=$output
output=${output}_maf01_hwe6
plink --bfile ${output_last} --maf 0.01 --hwe 0.000001 --make-bed --out ${output}

#run PC calculation

#prune
plink --bfile ${output} --indep-pairwise 200 100 0.2 --allow-no-sex --out ${output}_prune > /dev/null
plink --bfile ${output} --extract ${output}_prune.prune.in --make-bed --out ${output}_pruned > /dev/null

#Run smartpca
printf "genotypename: ${output}_pruned.bed
snpname: ${output}_pruned.bim
indivname: ${output}_pruned.fam
evecoutname: ${output}_pruned.pca.evec
evaloutname: ${output}_pruned.eigenvalues
altnormstyle: NO
numoutevec: 10
numoutlieriter: 0
numoutlierevec: 10
outliersigmathresh: 6
qtmode: 0" > ${output}_pruned.par
smartpca -p ${output}_pruned.par > ${output}_pruned_pccalc.log

#plot
Rscript plot_PCs_generate_ids_to_keep.R ${output}_pruned.pca.evec $race_sex_file

echo -e "PCs have been calculated and plots are saved here: ${output}_pruned.pdf and a list of NHW samples with outliers removed is here: ${output}_pruned_withoutoutliers.txt. Please check PC plots for outliers. If there are no outliers to remove, GWAS QC for $stem is complete!\n"

