#/bin/bash

#fail on error
set -e

#define usage
display_usage() {
    printf "GWAS QC post-imputation script

This script will unzip the imputation results (assuming password is saved in pass.txt in the same folder as the imputation results files and will perform standard post-imputation QC for our common variant pipeline. This includes filtering for R2, removing multi-allelic variants and filtering out variants for low MAF or HWE disequilibrium. Finally, PCs will be calculated on the final file-set.

Usage:
SCRIPTNAME.sh -o [output_stem] -i [imputation_results_folder] -r [race_sex_file] -s [snp_names_file]

stem = the beginning part of all QC'ed files including the full path to the folder in which they should be created
imputation_results_folder = the folder to which the imputation results have been downloaded which also contains a file called pass.txt with the password to unzip the imputation results files

race_sex_file = a file with FID and IID (corresponding to the fam file), 1 column indicating both race and ethnicity for PC plots, and another indicating sex for the sex check (1 for males, 2 for females, 0 if unknown), with NO header. Non-hispanic whites need to be indicated with 'White.' No other values in the race column must be fixed. This will be used to update sex in the fam file. 

snp_names_file = a file for converting the SNP names from imputation results to rs numbers. Must have 2 columns: imputation result SNP ids and rs numbers. Can have header but it will be ignored.

-h will display this message
"
        }

#parse options
while getopts 'o:i:r:s:h' flag; do
  case "${flag}" in
    o) output_stem="${OPTARG}" ;;
    i) input_fileset="${OPTARG}" ;;
    r) race_sex_file="${OPTARG}" ;;
    s) snp_names_file="${OPTARG}" ;;
    h) display_usage ; exit ;;
    \?|*) display_usage
       exit 1;;
  esac
done

#check to make sure necessary arguments are present                                                                                                    
if [ -z "$output_stem" ] || [ -z "$input_fileset" ] || [ -z "$race_sex_file" ] ;
then
    printf "Error: Necessary arguments not present!\n\n"
    display_usage
    exit 1
fi

#print out inputs
printf "GWAS QC Post-imputation Script

Output file path and stem for cleaned imputed files : $output_stem
Imputation results folder : $imputation_results_folder
Race/sex information file : $race_sex_file
File for SNP name conversion : $snp_names_file
"

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

#check to make sure password file is present
if [ ! -f "${imputation_results_folder}/pass.txt" ];
then
    printf "Imputation password file not found! Please save the password to unzip imputation results (found in the email notification from the imputation server of the job's completion) to ${imputation_results_folder}/pass.txt\n"
    exit 1
fi

#unzip imputation results using password
printf "Step 1 : Unzipping imputation results\n\n"
for i in $( ls ${imputation_results_folder}/*.zip );
do
    unzip -P $( cat ${imputation_results_folder}/pass.txt $i 
done

#filter imputation results for R2<0.8 and remove multi-allelic variants (multiple rows in vcf->bim)
printf "Step 2 : Filtering imputation results for R2<0.8 and multi-allelic variants\n\n"
for i in $(seq 1 22); do
    #restrict to variants with R2>=0.80
    plink2 --vcf ${imputation_results_folder}/chr${i}.dose.vcf.gz --exclude-if-info "R2<0.8" --make-bed --out ${output_stem}_chr${i}_temp > /dev/null

    #get list of duplicates SNPs (they have more than 1 row in the bim file) and remove
    awk '{ print $2 }' ${output_stem}_chr${i}_temp.bim | uniq -d > ${output_stem}_chr${i}_temp.dups ;
    plink2 --bfile ${output_stem}_chr${i}_temp --exclude ${output_stem}_chr${i}_temp.dups --make-bed --out ${output_stem}_chr${i}_temp_nodups > /dev/null
done

#print out numbers of variants
total_var=$(( $(grep "out of" ${output_stem}_chr*_temp.log | awk 'BEGIN { ORS="+" } { print $4 }' | sed 's/\(.*\)+/\1 /' ) ))
afterR2=$( cat ${output_stem}_chr*_temp_nodups.bim | wc -l )
nomulti=$( cat ${output_stem}_chr*_temp.bim | wc -l )
printf "$total_var variants after imputation, $afterR2 variants with R2>0.8, and $nomulti biallelic variants\n"

#Merge all chromosomes into one file
printf "Step 3 : Merging all chromosomes into one file\n\n"
#create merge file (emptying old version if present)
echo "" | tee ${output_stem}_merge_list.txt
for i in $( seq 1 22 );
do
    printf "${output_stem}_chr${i}_temp_nodups.bed  ${output_stem}_chr${i}_temp_nodups.bim  ${output_stem}_chr${i}_temp_nodups.fam\n" >> ${output_stem}_merge_list.txt ;
done

#merge individual chromosome files to be one large file
output=${output_stem}
plink --merge-list ${output_stem}_merge_list.txt --make-bed --out $output > /dev/null
grep 'pass filters and QC' ${output}.log

#update SNP names, sex, and perform standard SNP filtering
printf "Step 4 : Updating variant ids, sex in the fam file, and applying standard variant filters\n\n"
#update SNP names and add sex back into fam file
output_last=$output
output=${output}_names
plink --bfile ${output_last} --update-sex ${race_sex_file} 2 --update-name ${snp_names_file} --make-bed --out ${output} > /dev/null
grep -e "people updated" -e "values updated" ${output}.log

# SNP filters
output_last=$output
output=${output}_maf01_hwe6
plink --bfile ${output_last} --maf 0.01 --hwe 0.000001 --make-bed --out ${output} > /dev/null
grep -e "hwe: " -e "removed due to minor allele threshold" ${output}.log

# prune for heterozygosity check
printf "\n\nStep 5 : Pruning and running heterozygosity check\n"
plink --bfile ${output} --indep-pairwise 200 100 0.2 --allow-no-sex --out ${output}_prune > /dev/null
plink --bfile ${output} --output-missing-phenotype 1 --extract ${output}_prune.prune.in --make-bed --out ${output}_pruned > /dev/null
rm ${output}_prune.*
printf "$( wc -l < ${output}_pruned.bim ) variants out of $( wc -l < ${output}.bim ) left after pruning.\n"

##### heterozygosity check #####
plink --bfile ${output}_pruned --het --out ${output}_pruned_hetcheck > /dev/null
#plot and check for outliers
Rscript plot_het_check_outliers.R ${output}_pruned_hetcheck

#if there are outliers >6 sd from the F stat mean, remove them
if [ -f "${output}_pruned_hetcheck_outliers.txt" ];
then
    output_last=${output}
    output=${output}_nohetout
    plink --bfile $output_last --remove ${output}_pruned_hetcheck_outliers.txt --make-bed --out ${output} > /dev/null
    grep -e ' people pass filters and QC' ${output}.log
    printf "Output file: $output \n"
fi

#run PC calculation
printf "\n\nStep 5: Calculating post-imputation PCs\n\n"
sh calc_plot_PCs.sh -i $output -r $race_sex_file

printf "\nPlease check PC plots for outliers, remove if present, and recalculate PCs. If there are no outliers to remove, GWAS QC for this dataset is (probably) complete!
"

