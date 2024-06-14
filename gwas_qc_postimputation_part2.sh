#/bin/bash

#fail on error
set -e

#define usage
display_usage() {
    printf "GWAS QC post-imputation script, part 2

This script will conclude the post-imputation process, performing the MAF and HWE filters, checking sample heterozygosity, and calculating PCs on the final, cleaned files. Please check the PC plots for outliers.

Usage:
SCRIPTNAME.sh -o [input_file_stem] -g [preimputation_geno] -r [sample_ids] -l [subset_label] -c -p -m [plink_memory_flag]

input_file_stem = the full path and name (without bed/bim/fam) of the plink file set which resulted from the first stage of post-imputation QC. The files generated in this script will be saved in the same folder. 

preimputation_geno = the full path and stem to the cleaned final pre-imputation files to be merged back into the final files

sample_ids = the full path and name of the list of IDs for the subset to be QC'ed, output from the end of part 1 corresponding to one of the ancestral groups defined using SNPWeights

subset_label = the short label to identify this subset (will be added to the plink file name)

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
while getopts 'o:g:r:l:cpm:h' flag; do
  case "${flag}" in
    o) output_stem="${OPTARG}" ;;
    g) preimputation_geno="${OPTARG}" ;;
    r) sample_ids="${OPTARG}" ;;
    l) subset_label="${OPTARG}" ;;
    c) skip_cleanup='true' ;;
    m) plink_memory_limit="${OPTARG}";;
    p) skip_pccalc='true' ;;
    h) display_usage ; exit ;;
    \?|*) display_usage
       exit 1;;
  esac
done

########################################## Input validation ##################################################

#check to make sure necessary arguments are present
if [ -z "$output_stem" ] || [ -z "$preimputation_geno" ] || [ -z "$sample_ids" ] || [ -z "$subset_label" ];
then
    printf "Error: Necessary arguments not present! Please supply input stem, pre-imputation genotype stem, and list of samples to be QC'ed in.\n\n"
    display_usage
    exit 1
fi

#print out inputs
printf "GWAS QC Post-imputation Script, part 2

Output file path and stem for cleaned imputed files : $output_stem
Pre-imputation genotype file stem: $preimputation_geno
Sample ID list: $sample_ids
Subset label: $subset_label
"

#validate the genotyped files
if [ ! -f "${preimputation_geno}.fam" ];
then
    printf "Cannot see the preimputation genotype files ($preimputation_geno)! Please check the argument supplied to -g and try again!\n"
    exit 1
fi

if [ "$plink_memory_limit" ];
then 
    printf "Memory limit for plink calls: $plink_memory_limit \n"
    plink_memory_limit=$( echo "--memory $plink_memory_limit" )
fi

# get output folder
output_folder=${output_stem%/*}

################################## Step 1: Merge in pre-imputation genotypes #################################

printf "\nStep 1: Subset imputed and genotyped data to the samples in the current subset, filter for variant missingness, and merge imputed and genotyped data.\n"

printf "Subsetting imputed data and filtering variants for missingness...\n"
# Subset the imputed file
output=${output_stem}_${subset_label}_geno01
plink --bfile $output_stem --keep $sample_ids --geno 0.01 --make-bed --out $output $plink_memory_limit > /dev/null
samples=$( grep "pass filters and QC" ${output}.log | awk '{ print $4 }' )
variants=$( grep "pass filters and QC" ${output}.log | awk '{ print $1 }' )
printf "$samples and $variants remain in the imputed data.\n"

printf "Subsetting genotyped data and filtering variants for missingness...\n"
# Subset the genotyped file
geno_output=${preimputation_geno}_${subset_label}_geno01
plink --bfile $preimputation_geno --keep $sample_ids --geno 0.01 --make-bed --out $geno_output $plink_memory_limit > /dev/null
samples=$( grep "pass filters and QC" ${output}.log | awk '{ print $4 }' )
variants=$( grep "pass filters and QC" ${output}.log | awk '{ print $1 }' )
printf "$samples and $variants remain in the genotyped data.\n"

printf "Removing genotyped variants from imputed file...\n"
# make file with all the genotyped variant ids
awk '{ print $2 }' ${geno_output}.bim > ${output_folder}/${subset_label}_genotyped_variants.txt

# remove variants for which there are genotypes from the bim file
output_last=$output
output=${output}_nogeno
plink --bfile ${output_last} --exclude ${output_folder}/${subset_label}_genotyped_variants.txt --make-bed --out $output $plink_memory_limit > /dev/null
printf "$(grep -e "variants remaining" ${output}.log | awk '{ print $2 }') variants remaining after removing genotyped variants from imputation results.\n"

# merge the genotyped and imputed data
plink --bfile ${output} --bmerge ${geno_output} --make-bed --out ${output}_merged $plink_memory_limit > /dev/null

#check for complete sample overlap in the log and throw an error if there is not complete overlap
new_samples=$( grep "base dataset" ${output}_merged.log | awk 'NR==1{ print $3 }' )
if [ "$new_samples" != 0 ];
then 
    printf "There is incomplete overlap between the genotype and imputed fam files! Please ensure you input the correct genotypes and resolve any issues. (Hint: Check how many people whose ids were updated above and compare with the total. If fewer, then supply a race/sex file with all samples in the imputed files.)\n"
    exit 1
fi

#check for same position warnings
sameposwarnings=$( grep "Warning: Variants" ${output}_merged.log | head -n1 )
if [ ! -z "$sameposwarnings" ] ;
then 
    printf "Getting same position warnings. Removing those variants from the imputed dataset and re-attempting merge.\n"
    grep "Warning: Variants" ${output}_merged.log | awk  '{ print $3"\n"$5 }' | sed -e "s/'//g" >  ${output_folder}/${subset_label}_genotyped_variants_sameposwarnings.txt
    plink --bfile ${output} --exclude ${output_folder}/${subset_label}_genotyped_variants_sameposwarnings.txt --make-bed --out ${output}2 $plink_memory_limit > /dev/null
    printf "$(grep -e "variants remaining" ${output}.log | awk '{ print $2 }' ) variants at the same position in genotyped and imputed data. These variants were removed from the imputation results.\n"
    plink --bfile ${output}2 --bmerge ${geno_output} --make-bed --out ${output}_merged $plink_memory_limit > /dev/null
    
    #check for more warnings
    sameposwarnings=$( grep "Warning: Variants" ${output}_merged.log | head -n1 )
    if [ ! -z "$sameposwarnings" ] ;
    then
	printf "\nGetting more same position warnings! Please check and handle manually!\n"
    fi
fi
#update output variable
output=${output}_merged
printf "After merging in genotypes, there are $( grep "pass filters and QC" ${output}.log | sed 's/pass filters and QC.//' ).\nOutput file: $output\n"

########################################## Step 2: apply SNP filters ##################################################

printf "\nStep 2: Applying Hardy-Weinberg equilibrium and MAF filters.\n"

# SNP filters
output=${output_stem}_hwe6_maf01
plink --bfile ${output_stem} --hwe 0.000001 --maf 0.01 --make-bed --out ${output} $plink_memory_limit > /dev/null
grep -e "removed due to minor allele threshold" -e "hwe: " -e 'pass filters and QC' ${output}.log
printf "Output file: $output"

########################################## Step 3: heterozygosity check ##################################################

# prune for heterozygosity check
printf "\n\nStep 3: Pruning and running heterozygosity check\n"
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

########################################## Step 4: PC calculation (optional) ##################################################

if [ "$skip_pccalc" = 'false' ];
then 
    printf "\nStep 4: Calculating post-imputation PCs\n\n"

    # edit the memory limit variable for the purposes of supplying to the PC calculation script
    if [ "$plink_memory_limit" ]; then plink_memory_limit=$( echo $plink_memory_limit | sed 's/-memory/m/' ); fi 

    # Calculate PCs 
    sh calc_plot_PCs.sh -i $output $plink_memory_limit
    
    printf "\nPlease check PC plots for outliers, remove if present, and recalculate PCs. If there are no outliers to remove, GWAS QC for this dataset is (probably) complete!\n"

else
    printf "\nSkipping calculation of post-imputation PCs because the -p flag was supplied.\n\n"
fi

########################################## Cleanup (optional) ##################################################

if [ "$skip_cleanup" = 'false' ];
then 
    #do some clean-up
    printf "Doing some cleanup...\n"
    files_to_remove=$( find ${output_stem}*.bed | grep -v "${output}.bed" )
    rm $files_to_remove
else
    #tell the user to do the cleanup
    printf "The -c flag was specified to skip clean-up. Please remove any unnecessary *.bed files once complete to conserve space!\n"
fi
