#/bin/bash

#fail on error
set -e

#define usage
display_usage() {
    printf "GWAS QC post-imputation script

This script will unzip the imputation results (assuming password is saved in pass.txt in the same folder as the imputation results files and will perform standard post-imputation QC for our common variant pipeline. This includes filtering for R2, removing multi-allelic variants and filtering out variants for low MAF or HWE disequilibrium. Finally, PCs will be calculated on the final file-set.

Usage:
SCRIPTNAME.sh -o [output_stem] -i [imputation_results_folder] -r [race_sex_file] -s [snp_names_file] -z

output_stem = the beginning part of all QC'ed files including the full path to the folder in which they should be created

imputation_results_folder = the folder to which the imputation results have been downloaded which also contains a file called pass.txt with the password to unzip the imputation results files

race_sex_file = a file with FID and IID (corresponding to the fam file), 1 column indicating both race and ethnicity for PC plots, and another indicating sex for the sex check (1 for males, 2 for females, 0 if unknown), with NO header. Non-hispanic whites need to be indicated with 'White.' No other values in the race column must be fixed. This will be used to update sex in the fam file. 

snp_names_file = the file stem for converting the SNP names from imputation results to rs numbers. There should be one for each chromosome and each must have 2 columns: imputation result SNP ids and rs numbers. Can have header but it will be ignored.

-z indicates that the imputation results will need to be unzipped.

-h will display this message
"
        }

#parse options
do_unzip='false'
while getopts 'o:i:r:s:zh' flag; do
  case "${flag}" in
    o) output_stem="${OPTARG}" ;;
    i) imputation_results_folder="${OPTARG}" ;;
    r) race_sex_file="${OPTARG}" ;;
    s) snp_names_file="${OPTARG}" ;;
    z) do_unzip='true' ;;
    h) display_usage ; exit ;;
    \?|*) display_usage
       exit 1;;
  esac
done

#check to make sure necessary arguments are present                                                                                                    
if [ -z "$output_stem" ] || [ -z "$imputation_results_folder" ] || [ -z "$race_sex_file" ] || [ -z "$snp_names_file" ] ;
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
Stem for files for SNP name conversion : $snp_names_file
"
if [ "$do_unzip" = 'false' ];
then 
    printf "The -z was not specified so imputation results must already be unzipped.\n\n"
fi

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

################# Start the post-imputation QC ######################

if [ "$do_unzip" = 'true' ];
then
    #unzip imputation results using password
    printf "Step 1 : Unzipping imputation results\n\n"
    for i in $( ls ${imputation_results_folder}/*.zip );
    do
	unzip -P $( cat ${imputation_results_folder}/pass.txt ) $i -d $imputation_results_folder
    done
else
    #make sure the imputation results are actually unzipped
    if test ! -f ${imputation_results_folder}/chr1.dose.vcf.gz ;
    then 
	printf "The -z was not specified but cannot find the imputation results (ie ${imputation_results_folder}/chr1.dose.vcf.gz for each chromosome )! Please check whether imputation results have been unzipped. If they haven't, make sure the decryption password is saved in ${imputation_results_folder}/pass.txt and specify the -z flag.\n"
	exit 1
    else
	printf "Step 1 (unzipping the imputation results) is already complete. Proceeding the step 2...\n\n"
    fi
fi

#filter imputation results for R2<0.8 and remove multi-allelic variants (multiple rows in vcf->bim)
printf "Step 2 : Filtering imputation results for R2<0.8 and multi-allelic variants\n"
for i in $(seq 1 22); do
    #restrict to variants with R2>=0.80
    plink2 --vcf ${imputation_results_folder}/chr${i}.dose.vcf.gz --const-fid 0 --exclude-if-info "R2<0.8" --make-bed --out ${output_stem}_chr${i}_temp > /dev/null

    #get list of variants which have the same position (by base pair since this is just 1 chromosome) and remove
    dup_pos=$( awk '{ print $4 }' ${output_stem}_chr${i}_temp.bim | uniq -d )
    for j in $dup_pos ; do grep -w $j ${output_stem}_chr${i}_temp.bim | awk '{ print $2 }' >> ${output_stem}_chr${i}.dups ; done
    plink2 --bfile ${output_stem}_chr${i}_temp --exclude ${output_stem}_chr${i}.dups --make-bed --out ${output_stem}_chr${i}_temp_nodups > /dev/null

    #update variant names with rs numbers
    plink2 --bfile ${output_stem}_chr${i}_temp_nodups --update-name ${snp_names_file}_chr${i}.txt --make-bed --out ${output_stem}_chr${i}_temp_nodups_names > /dev/null
done

#print out numbers of variants
total_var=$(( $(grep "out of" ${output_stem}_chr*_temp.log | awk 'BEGIN { ORS="+" } { print $4 }' | sed 's/\(.*\)+/\1 /' ) ))
afterR2=$( cat ${output_stem}_chr*_temp.bim | wc -l )
nomulti=$( cat ${output_stem}_chr*_temp_nodups.bim | wc -l )
printf "$total_var variants after imputation, $afterR2 variants with R2>0.8, and $nomulti variants with unique positions\n\n"

#Merge all chromosomes into one file
printf "Step 3 : Merging all chromosomes into one file\n"
#create merge file (emptying old version if present)
echo "" | tee ${output_stem}_merge_list.txt
for i in $( seq 1 22 );
do
    printf "${output_stem}_chr${i}_temp_nodups_names\n" >> ${output_stem}_merge_list.txt ;
done

#merge individual chromosome files to be one large file
output=${output_stem}
plink --merge-list ${output_stem}_merge_list.txt --make-bed --out $output > /dev/null
grep 'pass filters and QC' ${output}.log

#update SNP names, sex, and perform standard SNP filtering
printf "Step 4 : Updating sex in the fam file and applying standard variant filters\n\n"
#update person ids using the race and sex file
output_last=$output
output=${output}_IDs
awk '{ print "0 "$1"_"$2" "$1" "$2 }' $race_sex_file > ${output_last}_update_ids.txt
plink --bfile $output_last --update-ids ${output_last}_update_ids.txt --make-bed --out $output > /dev/null

#update SNP names and add sex back into fam file
output_last=$output
output=${output}_sex
plink --bfile ${output_last} --update-sex ${race_sex_file} 2 --make-bed --out ${output} > /dev/null
grep -e "people updated" ${output}.log

# SNP filters
output_last=$output
output=${output}_maf01_hwe6
plink --bfile ${output_last} --maf 0.01 --hwe 0.000001 --make-bed --out ${output} > /dev/null
grep -e "hwe: " -e "removed due to minor allele threshold" -e 'pass filters and QC' ${output}.log

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

##### PC calculation ####
printf "\nStep 5: Calculating post-imputation PCs\n\n"

# Calculate PCs
sh calc_plot_PCs.sh -i $output -r $race_sex_file

printf "\nPlease check PC plots for outliers, remove if present, and recalculate PCs. If there are no outliers to remove, GWAS QC for this dataset is (probably) complete!\n"

