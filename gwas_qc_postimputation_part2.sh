#/bin/bash

#fail on error
set -e

#define usage
display_usage() {
    printf "GWAS QC post-imputation script, part 2

This script will conclude the post-imputation process, performing the MAF and HWE filters, checking sample heterozygosity, and calculating PCs on the final, cleaned files. Please check the PC plots for outliers.

Usage:
SCRIPTNAME.sh -o [input_file_stem] -r [race_file] -c -p -m [plink_memory_flag]

input_file_stem = the full path and name (without bed/bim/fam) of the plink file set which resulted from the first stage of post-imputation QC. The files generated in this script will be saved in the same folder. 

race_file = a file with FID and IID (corresponding to the fam file) and 1 column indicating both race and ethnicity for PC plots with NO header. Values in this file will typically be the SNPWeights-derived ancestral groups. 

-c will skip the clean-up at the end of the script which removes intermediate *.bed files.

-p will skip the PC calculation at the end (only should be done if this set will be merged with other sets in the same cohort).

plink_memory_limit (optional) = argument indicating the memory limit for plink to use rather than the default of half the RAM. This is useful for running this step of QC locally.

-h will display this message
"
        }

#parse options
do_unzip='false'
skip_first_filters='false'
skip_cleanup='false'
skip_pccalc='false'
while getopts 'o:r:cpm:h' flag; do
  case "${flag}" in
    o) output_stem="${OPTARG}" ;;
    r) race_file="${OPTARG}" ;;
    c) skip_cleanup='true' ;;
    m) plink_memory_limit="${OPTARG}";;
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
Race information file : $race_file
"

if [ "$plink_memory_limit" ];
then 
    printf "Memory limit for plink calls: $memory_limit \n"
    plink_memory_limit=$( echo "--memory $plink_memory_limit" )
fi

printf "\nStep 1: Applying Hardy-Weinberg equilibrium and MAF filters.\n"

# SNP filters
output=${output_stem}_hwe6_maf01
plink --bfile ${output_stem} --hwe 0.000001 --maf 0.01 --make-bed --out ${output} $plink_memory_limit > /dev/null
grep -e "removed due to minor allele threshold" -e "hwe: " -e 'pass filters and QC' ${output}.log
printf "Output file: $output"

# prune for heterozygosity check
printf "\n\nStep 2: Pruning and running heterozygosity check\n"
plink --bfile ${output} --indep-pairwise 200 100 0.2 --allow-no-sex --out ${output}_prune $plink_memory_limit > /dev/null
plink --bfile ${output} --output-missing-phenotype 1 --extract ${output}_prune.prune.in --make-bed --out ${output}_pruned $plink_memory_limit > /dev/null
rm ${output}_prune.*
printf "$( wc -l < ${output}_pruned.bim ) variants out of $( wc -l < ${output}.bim ) left after pruning.\n"

##### heterozygosity check #####
plink --bfile ${output}_pruned --het --out ${output}_pruned_hetcheck $plink_memory_limit > /dev/null
#plot and check for outliers
Rscript plot_het_check_outliers.R ${output}_pruned_hetcheck

#if there are outliers >6 sd from the F stat mean, remove them
if [ -f "${output}_pruned_hetcheck_outliers.txt" ];
then
    output_last=${output}
    output=${output}_nohetout
    plink --bfile $output_last --remove ${output}_pruned_hetcheck_outliers.txt --make-bed --out ${output} $plink_memory_limit > /dev/null
    grep -e ' people pass filters and QC' ${output}.log
    printf "Output file: $output \n"
fi

##### PC calculation ####
if [ "$skip_pccalc" = 'false' ];
then 
    printf "\nStep 3: Calculating post-imputation PCs\n\n"

    # edit the memory limit variable for the purposes of supplying to the PC calculation script
    if [ "$plink_memory_limit" ]; then plink_memory_limit=$( echo $plink_memory_limit | sed 's/-memory/m/' ); fi 

    # Calculate PCs, coloring based on ancestry categories if present
    if [ ! -z "$race_file" ]; 
    then
	sh calc_plot_PCs.sh -i $output -r $race_file $plink_memory_limit
    else
	sh calc_plot_PCs.sh -i $output $plink_memory_limit
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
