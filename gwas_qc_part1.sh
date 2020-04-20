#!/bin/bash
# Author: Vaibhav Janve 2020-04-07
# Author 2: Emily Mahoney 2020-04-09

stem=$1
input_fileset=$2
output_dir=$3
race_sex_file=$4

#fail on error
set -e

#define usage
display_usage() {
    printf "GWAS QC Part 1 script
Completes the first stage in standard GWAS QC, including initial variant and person filters, relatedness and sex checks, restriction to autosomes, and PC calculation.

Usage:
SCRIPTNAME.sh [stem] [input_fileset] [output_dir] [race_sex_file]

stem = the beginning part of all QC'ed files
input_fileset = the full path and file stem for the raw plink set '*[bed,bim,fam]'
output_dir = where the intermediate QC files will be created
race_sex_file = a file with FID and IID (corresponding to the fam file), 1 column indicating both race and ethnicity for PC plots, and another indicating sex for the sex check (1 for males, 2 for females, 0 if unknown), with NO header. Non-hispanic whites need to be indicated with 'White.' No other values in the race column must be fixed.

Note: assumes PLINK 1.9 is available and that this is being run from the scripts folder.
"
        }
#check if there are at least 3 arguments
if [  $# -lt 4 ]
        then
                display_usage
                exit 1
fi


#print out inputs
echo -e "GWAS QC Part 1 Script\n"
echo "stem : "$stem
echo "raw input data : "$input_fileset
echo "QC output directory : "$output_dir
echo "file with race/ethnicity and sex information : "$race_sex_file

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


##### initial SNP filters #####
echo -e "\n\nStep 1: Remove SNPs with >5% missingness or with MAF <0.01\n"
output=$( echo ${output_dir}/${stem}_geno05_maf01 )
plink --bfile $input_fileset --geno 0.05 --maf 0.01 --make-bed --out $output > /dev/null

grep 'loaded from' $output.log  | sed  's/loaded from.*//' | head -n2
grep -e 'missing genotype data' -e 'minor allele threshold' -e ' people pass filters and QC' $output.log
echo -e "Output file: $output \n"

##### person missingness #####
echo -e "\nStep 2: remove subjects w/ >1% missingness\n"
output_last=$output
output=${output}_mind01
plink --bfile $output_last --mind 0.01 --make-bed --out $output > /dev/null

grep -e 'removed due to missing genotype data' -e ' people pass filters and QC' $output.log
echo -e "Output file: $output \n"

##### Relatedness #####
echo -e "\nStep 3: Remove related individuals\n"
echo -e "   Identify relateds"
plink --bfile $output --genome full unbounded nudge --min 0.20 --out ${output}_relatedness > /dev/null

#print out # or related individuals at each threshold
for pi_hat in 0.9 0.5 0.25; do
    echo "Pairs above pi-hat of $pi_hat: "$(awk -v val="$pi_hat" '{ if (NR!=1 && $10 > val ) print }' ${output}_relatedness.genome | wc -l)
done

#print out all basically identical pairs which will be removed entirely to allow for quick scanning for intended sample duplicates (re-genotyped samples, the unlikely event of actual identical twins, ect)
if [ "$( awk '{ if(NR==1 || $10 > 0.9 ) print }' ${output}_relatedness.genome | wc -l )" -gt 1 ];
then
    echo -e "Sample pairs with pi-hat above 0.9 which will be entirely removed:"
    awk '{ if(NR==1 || $10 > 0.9 ) print }' ${output}_relatedness.genome
fi

# Rscript for decisions
Rscript get_related_ids.R ${output}_relatedness

# remove selected ids from the last generated genotype file set
echo -e "   Remove relateds"
echo -e "Removing $( wc -l ${output}_relatedness_related_ids.txt ) individuals for relatedness."
output_last=$output
output=${output}_norelated
plink --bfile $output_last --remove ${output_last}_relatedness_related_ids.txt --make-bed --out $output > /dev/null
grep ' people pass filters and QC' $output.log
echo -e "Output file: $output \n"


##### sex check #####
echo -e "\nStep 4: sex check\n"

#only update sex if there is something more than missing values for sex in the provided file
if [ "$( awk '{ print $4 }' $race_sex_file | sort -u )" != 0 ];
then 
    output_last=$output
    output=${output}_sex
    plink --bfile $output_last --update-sex $race_sex_file 2 --make-bed --out $output > /dev/null
    echo -e "Updated sex from the provided file: $race_sex_file"
else
    echo -e "No non-missing sex information in the provided file (${race_sex_file}). Using sex from fam file."
fi
#do the check
plink --bfile $output --check-sex --out ${output}_checking_sex > /dev/null
grep -e 'check-sex: ' ${output}_checking_sex.log 

#write mismatched sex iids in text file
awk '{ if($5=="PROBLEM" && $4 != 0 && $3 != 0) print $1" "$2 }' ${output}_checking_sex.sexcheck > ${output}_mismatched_sex_ids.txt
echo -e "$( wc -l < ${output}_mismatched_sex_ids.txt ) real sex mismatches (e.g. not ambiguous in the fam or indeterminate based on SNPs)"
echo -e "$( awk '{ if($3 == 0) print }' ${output}_checking_sex.sexcheck | wc -l ) out of $( wc -l < ${output}.fam ) samples are missing sex."

#remove individuals in text file
if [ $( wc -l < ${output}_mismatched_sex_ids.txt ) -gt 0 ];
then
    sex_mismatch_file=${output}_mismatched_sex_ids.txt
    output_last=$output
    output=${output}_nomismatchedsex
    plink --bfile $output_last --remove ${sex_mismatch_file} --make-bed --out $output > /dev/null
    grep -e ' people pass filters and QC' ${output}.log
    echo -e "Output file: $output \n"
    
fi

##### restrict to autosomes #####
echo -e "\nStep 5 : restrict to autosomes\n"
output_last=$output
output=${output}_keep_autosomes
plink --bfile $output_last --chr 1-22 --make-bed --out $output > /dev/null

#Get numbers of removed variants
# total
non_autosomal_var=$(( $(wc -l <${output_last}.bim) - $(wc -l <${output}.bim) ))
# unmapped (chr=0)
unmapped_var=$( awk '{ if($1 == 0) print }' $output_last.bim | wc -l)
# sex variants (chr=23,24,25)
sex_var=$( awk '{ if($1 == 23 || $1 == 24 || $1 == 25) print }' $output_last.bim | wc -l)
# mitochondrial (chr=26)
mito_var=$(awk '{ if($1 == 26) print }' $output_last.bim | wc -l)

echo -e "Removed $non_autosomal_var variants total ($unmapped_var unmapped, $sex_var sex chr, $mito_var mitochondrial)"
grep -e ' people pass filters and QC' ${output}.log
echo -e "Output file: $output \n"

##### prep for PC calculations #####

echo -e "\nStep 6: Checking ids and pruning dataset for heterozygosity check and PC calculation\n"

#do quick check of length of ids
Rscript check_id_length.R $output

#get names of new person/SNP ids files (will only have been created if there need to be some updates)
fam_ids=${output}_dummy_famids.txt
bim_ids=${output}_shorter_bimids.txt

if [ -f "$fam_ids" ];
then
    output_last=$output
    output=${output}_dummy_famids
    plink --bfile $output_last --update-ids $fam_ids --make-bed --out ${output} > /dev/null
fi

if [ -f "$bim_ids" ];
then
    output_last=$output
    output=${output}_ids
    plink --bfile $output_last --update-name $bim_ids --make-bed --out ${output} > /dev/null
fi

#prune and set phenotypes to 1
plink --bfile ${output} --indep-pairwise 200 100 0.2 --allow-no-sex --out ${output}_prune > /dev/null
plink --bfile ${output} --output-missing-phenotype 1 --extract ${output}_prune.prune.in --make-bed --out ${output}_pruned > /dev/null
echo -e "$( wc -l < ${output}_pruned.bim ) variants out of $( wc -l < ${output}.bim ) left after pruning."

echo -e "\nStep 7: Running heterozygosity check\n"
##### heterozygosity check #####
plink --bfile ${output}_pruned --het --out ${output}_pruned_hetcheck > /dev/null
#plot and check for outliers
Rscript plot_het_check_outliers.R ${output}_pruned_hetcheck

#if there are outliers, remove them
if [ -f "${output}_pruned_hetcheck_outliers.txt" ];
then
    output_last=${output}
    output=${output}_nohetout
    plink --bfile $output_last --remove ${output}_pruned_hetcheck_outliers.txt --make-bed --out ${output} > /dev/null
    grep -e ' people pass filters and QC' ${output}.log
    echo -e "Output file: $output \n"
    echo -e "Redoing pruning in the new fileset in preparation for PC calculation."

    #redo pruning for PC calculation
    plink --bfile ${output} --indep-pairwise 200 100 0.2 --allow-no-sex --out ${new_file_name}_prune > /dev/null
    plink --bfile ${output} --extract ${new_file_name}_prune.prune.in --make-bed --out ${output}_pruned > /dev/null
    echo -e "$( wc -l < ${output}_pruned.bim ) variants out of $( wc -l < ${output}.bim ) left after pruning."
fi


##### Calculate PCs #####
echo -e "\nStep 8: Running PC calculation with smartpca\n"

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

#step6.2: create plots
Rscript plot_PCs_generate_ids_to_keep.R ${output}_pruned.pca.evec $race_sex_file

printf "PCs calculated and plots are saved here: ${output}_pruned.pdf. A file with ids for NHW who were not outliers on PCs1-3 is written out for your convenience if all outliers should be removed: ${output}_pruned_withoutoutliers.txt Please check PC plots and decide what individuals to remove before proceeding to imputation preparation. \n"
