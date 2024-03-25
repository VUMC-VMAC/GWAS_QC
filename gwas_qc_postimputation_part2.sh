#/bin/bash

#fail on error
set -e

#define usage
display_usage() {
    printf "GWAS QC post-imputation script, part 2

This script will conclude the post-imputation process, performing the MAF and HWE filters, checking sample heterozygosity, and calculating PCs on the final, cleaned files. Please check the PC plots for outliers.

Usage:
SCRIPTNAME.sh -o [input_file_stem] -r [race_sex_file] -c -p

input_file_stem = the full path and name (without bed/bim/fam) of the plink file set which resulted from the first stage of post-imputation QC. The files generated in this script will be saved in the same folder. 

race_sex_file (optional) = a file with FID and IID (corresponding to the fam file), 1 column indicating both race and ethnicity for PC plots, and another indicating sex for the sex check (1 for males, 2 for females, 0 if unknown), with NO header. Non-hispanic whites need to be indicated with 'White.' No other values in the race column must be fixed. This will be used to update sex in the fam file. 

-c will skip the clean-up at the end of the script which removes intermediate *.bed files.
-p will skip the PC calculation at the end (only should be done if this set will be merged with other sets in the same cohort).
-h will display this message
"
        }

#parse options
do_unzip='false'
skip_first_filters='false'
skip_cleanup='false'
skip_pccalc='false'
while getopts 'o:r:cph' flag; do
  case "${flag}" in
    o) output_stem="${OPTARG}" ;;
    r) race_sex_file="${OPTARG}" ;;
    c) skip_cleanup='true' ;;
    p) skip_pccalc='true' ;;
    h) display_usage ; exit ;;
    \?|*) display_usage
       exit 1;;
  esac
done

#check to make sure necessary arguments are present                                                                                              
if [ -z "$output_stem" ];
then
    printf "Error: Please supply an input file stem!! Note that this should be the stem of the genotype files after performing SNPweights and applying any appropriate subsetting.\n\n"
    display_usage
    exit 1
fi

#print out inputs
printf "GWAS QC Post-imputation Script, part 2

Output file path and stem for cleaned imputed files : $output_stem
Race/sex information file : $race_sex_file
"

printf "\nStep 1: Applying Hardy-Weinberg equilibrium and MAF filters.\n"

# SNP filters
output=${output_stem}_hwe6_maf01
plink --bfile ${output_stem} --hwe 0.000001 --maf 0.01 --make-bed --out ${output} > /dev/null
grep -e "removed due to minor allele threshold" -e "hwe: " -e 'pass filters and QC' ${output}.log
printf "Output file: $output"

# prune for heterozygosity check
printf "\n\nStep 2: Pruning and running heterozygosity check\n"
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
if [ "$skip_pccalc" = 'true' ];
then 
    printf "\nStep 3: Calculating post-imputation PCs\n\n"

    # Calculate PCs, coloring based on ancestry categories if present
    if [ ! -z "$race_sex_file" ]; 
    then
	sh calc_plot_PCs.sh -i $output -r $race_sex_file
    else
	sh calc_plot_PCs.sh -i $output 
    fi

    printf "\nPlease check PC plots for outliers, remove if present, and recalculate PCs. If there are no outliers to remove, GWAS QC for this dataset is (probably) complete!\n"

else
    printf "\nSkipping calculation of post-imputation PCs because the -p flag was supplied.\n\n"
fi

if [ "$skip_cleanup" = 'true' ];
then 
    #do some clean-up
    printf "PC calculation complete. Doing some cleanup...\n"
    files_to_remove=$( find ${output_stem}*.bed | grep -v "${output}.bed" )
    rm $files_to_remove
else
    #tell the user to do the cleanup
    printf "PC calculation complete. The -c flag was specified to skip clean-up. Please remove any unnecessary *.bed files once complete to conserve space!\n"
fi
