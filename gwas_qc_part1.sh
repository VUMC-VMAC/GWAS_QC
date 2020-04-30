#!/bin/bash
# Author: Vaibhav Janve 2020-04-07
# Author 2: Emily Mahoney 2020-04-09

#fail on error
set -e

#define usage
display_usage() {
    printf "GWAS QC Part 1 script

Completes the first stage in standard GWAS QC, including initial variant and person filters, relatedness and sex checks, restriction to autosomes, and PC calculation.

Usage:
SCRIPTNAME.sh -o [output_stem] -i [input_fileset] -r [race_sex_file] -G [stem_1000G]


output_stem = the beginning part of all QC'ed files, including the full file path to the directory where the files are to be saved

input_fileset = the full path and file stem for the raw plink set '*[bed,bim,fam]'

race_sex_file = a file with FID and IID (corresponding to the fam file), 1 column indicating both race and ethnicity for PC plots, and another indicating sex for the sex check (1 for males, 2 for females, 0 if unknown), with NO header. Non-hispanic whites need to be indicated with 'White.' No other values in the race column must be fixed; however, the race column must not include spaces.

stem_1000G (optional) = the full path and stem to the 1000G genotype files in plink format. There must also be a file with FID, IID, race with the same stem and _race.txt as the suffix (ie for a plink file set like this: all_1000G.bed, all_1000G.bim, all_1000G.fam the race file would be like this all_1000G_race.txt)

-h will show this usage
"
        }

while getopts 'o:i:r:G:h' flag; do
  case "${flag}" in
    o) output_stem="${OPTARG}" ;;
    i) input_fileset="${OPTARG}" ;;
    r) race_sex_file="${OPTARG}" ;;
    G) stem_1000G="${OPTARG}" ;;
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
printf "GWAS QC Part 1 Script

Input data : $input_fileset
Output stem : $output_stem 
File with race/ethnicity and sex information : $race_sex_file
"
#print message if 1000G dataset is not specified
if [ -z "$stem_1000G" ];
then
    printf "No location was specified for the 1000G data, so no PCs will be calculated including them!\n\n"
else
    echo $stem_1000G
    printf "1000G data for PC calculation : ${stem_1000G}\n\n"
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


##### initial SNP filters #####
printf '%s\n\n' "Step 1: Remove SNPs with >5% missingness or with MAF <0.01"
output=$( printf ${output_stem}_geno05_maf01 )
plink --bfile $input_fileset --geno 0.05 --maf 0.01 --make-bed --out $output > /dev/null

grep 'loaded from' $output.log  | sed  's/loaded from.*//' | head -n2
grep -e 'missing genotype data' -e 'minor allele threshold' -e ' people pass filters and QC' $output.log
printf "Output file: $output \n"

##### person missingness #####
printf '%s\n' "Step 2: remove subjects w/ >1% missingness"
output_last=$output
output=${output}_mind01
plink --bfile $output_last --mind 0.01 --make-bed --out $output > /dev/null

grep -e 'removed due to missing genotype data' -e ' people pass filters and QC' $output.log
printf "Output file: $output \n"

##### Relatedness #####
printf "\nStep 3: Calculate relatedness and remove related individuals\n"
plink --bfile $output --genome full unbounded nudge --min 0.20 --out ${output}_relatedness > /dev/null

#print out # or related individuals at each threshold
for pi_hat in 0.9 0.5 0.25; do
    n_rel=$(awk -v val="$pi_hat" '{ if (NR!=1 && $10 > val ) print }' ${output}_relatedness.genome | wc -l)
    printf "Pairs above pi-hat of $pi_hat: $n_rel \n"
done

#print out all basically identical pairs which will be removed entirely to allow for quick scanning for intended sample duplicates (re-genotyped samples, the unlikely event of actual identical twins, ect)
if [ "$( awk '{ if(NR==1 || $10 > 0.9 ) print }' ${output}_relatedness.genome | wc -l )" -gt 1 ];
then
    printf "Sample pairs with pi-hat above 0.9 which will be entirely removed:\n"
    awk '{ if(NR==1 || $10 > 0.9 ) print $1" "$2" "$3" "$4" "$10 }' ${output}_relatedness.genome
fi

# Rscript for decisions
Rscript get_related_ids.R ${output}_relatedness

# remove selected ids from the last generated genotype file set
printf "Removing $( wc -l ${output}_relatedness_related_ids.txt ) individuals for relatedness.\n"
output_last=$output
output=${output}_norelated
plink --bfile $output_last --remove ${output_last}_relatedness_related_ids.txt --make-bed --out $output > /dev/null
grep ' people pass filters and QC' $output.log
printf "Output file: $output \n"


##### sex check #####
printf "\nStep 4: sex check\n"

#only update sex if there is something more than missing values for sex in the provided file
if [ "$( awk '{ print $4 }' $race_sex_file | sort -u | tr -d '[:space:]' )" != 0 ];
then 
    output_last=$output
    output=${output}_sex
    plink --bfile $output_last --update-sex $race_sex_file 2 --make-bed --out $output > /dev/null
    printf "Updated sex from the provided file: $race_sex_file \n"
else
    printf "No non-missing sex information in the provided file (${race_sex_file}). Using sex from fam file.\n"
fi
#do the check
plink --bfile $output --check-sex --out ${output}_checking_sex > /dev/null
grep -e 'check-sex: ' ${output}_checking_sex.log 

#write mismatched sex iids in text file
awk '{ if($5=="PROBLEM" && $4 != 0 && $3 != 0) print $1" "$2 }' ${output}_checking_sex.sexcheck > ${output}_mismatched_sex_ids.txt
printf "$( wc -l < ${output}_mismatched_sex_ids.txt ) real sex mismatches (e.g. not ambiguous in the fam or indeterminate based on SNPs)
$( awk '{ if($3 == 0) print }' ${output}_checking_sex.sexcheck | wc -l ) out of $( wc -l < ${output}.fam ) samples are missing sex.\n"

#remove individuals in text file
if [ $( wc -l < ${output}_mismatched_sex_ids.txt ) -gt 0 ];
then
    sex_mismatch_file=${output}_mismatched_sex_ids.txt
    output_last=$output
    output=${output}_nomismatchedsex
    plink --bfile $output_last --remove ${sex_mismatch_file} --make-bed --out $output > /dev/null
    grep -e ' people pass filters and QC' ${output}.log
    printf "Output file: $output \n"
    
fi

##### restrict to autosomes #####
printf "\nStep 5 : restrict to autosomes\n"
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

printf "Removed $non_autosomal_var variants total ($unmapped_var unmapped, $sex_var sex chr, $mito_var mitochondrial)\n"
grep -e ' people pass filters and QC' ${output}.log
printf "Output file: $output \n"

##### prep for PC calculations #####

printf "\nStep 6: Checking ids and pruning dataset for heterozygosity check and PC calculation\n"

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
printf "$( wc -l < ${output}_pruned.bim ) variants out of $( wc -l < ${output}.bim ) left after pruning.\n"

##### heterozygosity check #####
printf "\nStep 7: Running heterozygosity check\n"
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
    printf "Output file: $output \n\n"
    printf "Redoing pruning in the new fileset in preparation for PC calculation within this dataset.\n"

    #redo pruning for PC calculation
    plink --bfile ${output} --indep-pairwise 200 100 0.2 --allow-no-sex --out ${new_file_name}_prune > /dev/null
    plink --bfile ${output} --output-missing-phenotype 1 --extract ${new_file_name}_prune.prune.in --make-bed --out ${output}_pruned > /dev/null
    printf "$( wc -l < ${output}_pruned.bim ) variants out of $( wc -l < ${output}.bim ) left after pruning.\n"
fi

##### Calculate PCs within this dataset #####

printf "\nStep 8: Running PC calculation with smartpca\n"

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

#create plots
#arguments: PCA file, race/sex file, NA or location of 1000G race file, yes/no to indicate whether to output a file for outlier exclusion
if [ ! -z "$stem_1000G" ]
then
    #if 1000G is set, then don't output a file for exclusion and print message about technical issues
    Rscript plot_PCs_generate_ids_to_keep.R ${output}_pruned.pca.evec $race_sex_file NA no
    printf "PCs calculated within the current dataset and plots are saved here: ${output}_pruned_PCplots.pdf. 
Please scan PC plots for any potential technical issues. \n"
else
    #if 1000f is not set, then output a file for default exclusions, print out message about it, and exit the script
    Rscript plot_PCs_generate_ids_to_keep.R ${output}_pruned.pca.evec $race_sex_file NA yes
    printf "\nPCs calculated within the current dataset and plots are saved here: ${output}_pruned_PCplots.pdf. 
A file with ids for NHW who were not >5 SD from the mean in European samples on PCs1-3 is written out for your convenience if all outliers should be removed: ${output}_nooutliers.txt 
Please check PC plots and decide what individuals to remove before proceeding to imputation preparation. \n"
    exit
fi


#### prep for PCs with 1000G data ###
printf "\nStep 9: Merging with 1000G for PC calculation to make ancestry filtering decisions\n"

#get variants in current dataset
awk '{ print $2 }' ${output}.bim > ${output}_snps.txt

#subset 1000G dataset to only variants present here
plink --bfile ${stem_1000G}  --output-missing-phenotype 1 --extract ${output}_snps.txt --make-bed --out ${output}_1000G_overlapping > /dev/null
num_1000G_snps=$( wc -l < ${output}_1000G_overlapping.bim )
printf "$num_1000G_snps overlap between 1000G and the current dataset\n"

#attempt merge (keep the script from failing when this command fails)
plink --bfile ${output} --bmerge ${output}_1000G_overlapping --allow-no-sex --make-bed --out ${output_stem}_1000G || true &> /dev/null
with_1000G=${output_stem}_1000G

if [ -f "${output_stem}_1000G-merge.missnp" ];
then
    #remove the mismatching variants
    plink --bfile ${output} --exclude ${with_1000G}-merge.missnp --allow-no-sex --make-bed --out ${output_stem}_formerge > /dev/null
    plink --bfile ${output}_1000G_overlapping --exclude ${with_1000G}-merge.missnp --allow-no-sex --make-bed --out ${output_stem}_1000G_formerge  > /dev/null
    plink --bfile ${output_stem}_formerge --bmerge ${output_stem}_1000G_formerge --allow-no-sex --make-bed --out ${with_1000G}  > /dev/null

    #filter to overlapping
    plink --bfile ${output_stem}_1000G --geno 0.05 --allow-no-sex --make-bed --out ${output_stem}_1000G_geno05  > /dev/null
    printf "$( wc -l < ${output_stem}_1000G_geno05.bim ) variants in dataset merged with 1000G\n"
    with_1000G=${output_stem}_1000G_geno05
fi

#prune
plink --bfile ${with_1000G} --indep-pairwise 200 100 0.2 --allow-no-sex --out ${with_1000G}_forprune > /dev/null
plink --bfile ${with_1000G} --output-missing-phenotype 1 --extract ${with_1000G}_forprune.prune.in --allow-no-sex --make-bed --out ${with_1000G}_pruned > /dev/null
printf "$( wc -l < ${with_1000G}_pruned.bim ) variants out of $( wc -l < ${with_1000G}.bim ) left after pruning.\n"

##### Calculate PCs #####
printf "\nStep 10: Running PC calculation including 1000G samples with smartpca\n"

#Run smartpca
printf "genotypename: ${with_1000G}_pruned.bed
snpname: ${with_1000G}_pruned.bim
indivname: ${with_1000G}_pruned.fam
evecoutname: ${with_1000G}_pruned.pca.evec
evaloutname: ${with_1000G}_pruned.eigenvalues
altnormstyle: NO
numoutevec: 10
numoutlieriter: 0
numoutlierevec: 10
outliersigmathresh: 6
qtmode: 0" > ${with_1000G}_pruned.par
smartpca -p ${with_1000G}_pruned.par > ${with_1000G}_pruned_pccalc.log

#create plots
#arguments: PCA file, race/sex file, yes/no to indicate whether to output a file for outlier exclusion
Rscript plot_PCs_generate_ids_to_keep.R ${with_1000G}_pruned.pca.evec $race_sex_file ${stem_1000g}_race.txt yes

printf "PCs calculated including 1000G samples and plots are saved here: ${with_1000G}_pruned_PCplots.pdf. 
A file with ids for NHW who were not >5 SD from the mean in European samples on PCs1-3 is written out for your convenience if all outliers should be removed: ${with_1000G}_nooutliers.txt 
Please check PC plots and decide what individuals to remove before proceeding to imputation preparation. \n"
