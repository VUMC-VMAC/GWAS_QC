#!/bin/bash
# Author: Vaibhav Janve 2020-04-07
# Author 2: Emily Mahoney 2020-04-09

#fail on error
set -e

#define usage
display_usage() {
    printf "GWAS QC Pre-Imputation Script (Relateds)

Completes the first stage in standard GWAS QC, including initial variant and person filters, relatedness and sex checks, restriction to autosomes, HWE filtering, and preparation for upload to the imputation server. This is a modified version of the main script to retain related individuals.

Usage:
SCRIPTNAME.sh -i [input_fileset] -o [output_stem] -f [sex_file] -R [ref_file_stem] -b [input genome build] -m [plink_memory_limit] -n -c 

output_stem = the beginning part of all QC'ed files, including the full file path to the directory where the files are to be saved

input_fileset = the full path and file stem for the raw plink set '*[bed,bim,fam]'

sex_file = a file with FID and IID (corresponding to the fam file), 1 column indicating sex for the sex check (1 for males, 2 for females, 0 if unknown), with NO header.

ref_file_stem = the full file path and stem for the reference panel. Assumes it is split by chromosome and named STEM_chr*.txt.gz

input_genome_build (optional) = argument indicating the build of the input dataset; default assumption is b37 but can supply b34, b35, b36, or b38. This will impact the liftOver step.

plink_memory_limit (optional) = argument indicating the memory limit for plink to use rather than the default of half the RAM. This is useful for running this step of QC locally.

-n = optional argument set to indicate not to exclude variants for not being in or not matching the reference panel; default is to exclude

-c = optional argument to skip clean-up

-h will show this usage
"
        }

# set default build 
build=b37
while getopts 'o:i:f:R:b:m:nch' flag; do
  case "${flag}" in
    o) output_stem="${OPTARG}" ;;
    i) input_fileset="${OPTARG}" ;;
    f) sex_file="${OPTARG}" ;;
    R) ref_file_stem="${OPTARG}" ;;
    b) build="${OPTARG}" ;;
    m) plink_memory_limit="${OPTARG}";;
    n) noexclude='true' ;;
    c) skip_cleanup='true' ;;
    h) display_usage ; exit ;;
    \?|*) display_usage
       exit 1;;
  esac
done

#check to make sure necessary arguments are present
if [ -z "$output_stem" ] || [ -z "$input_fileset" ] || [ -z "$sex_file" ] || [ -z "$ref_file_stem" ];
then
    printf "Error: Necessary arguments not present!\n\n"
    display_usage
    exit 1
fi

printf "GWAS QC Pre-imputation (Including Related Individuals)\n"

# Print out singularity information for reproducability
[ ! -z "$SINGULARITY_CONTAINER" ] && printf "\nSingularity image: $SINGULARITY_CONTAINER"
[ ! -z "$SINGULARITY_COMMAND" ] && printf "\nSingularity shell: $SINGULARITY_COMMAND"
[ ! -z "$SINGULARITY_BIND" ] && printf "\nMapped directories:\n $SINGULARITY_BIND\n" | sed 's/,/\n /g'

#print out inputs
printf "\n\nInput data : $input_fileset
Output stem : $output_stem 
File with sex information : $sex_file
Reference panel SNP file : $ref_file_stem
"

# check the build argument
if [ "$build" != "b34" ] && [ "$build" != "b35" ] && [ "$build" != "b36" ] && [ "$build" != "b37" ] && [ "$build" != "b38" ]; 
then 
    printf "Error: Invalid build argument supplied ($build). Please supply one of the following: b34, b35, b36, b37, b38.\n"
    exit 1
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

if [ "$plink_memory_limit" ];
then 
    printf "Memory limit for plink calls: $plink_memory_limit \n"
    plink_memory_limit=$( echo "--memory $plink_memory_limit" )
fi

#get file path for input
input_path=${input_fileset%/*}
output_path=${output_stem%/*}

#get input file name to be the stem for the output
input_stem=${input_fileset##*/}

###################################### Begin QC process #########################################

###### format plinkset: update varIDs #####
printf "\nStep 0: update varIDs to standard CHR_POS_REF_ALT format or CHR_POS_I_D for indels and remove duplicates.\n"
output=${output_stem}_stdnames
## generate an ID with just variant info rather than rsID and shorten IDs that are too long
awk '{if( length($5)>25 || length($6)>25) { print $2,$1"_"$4"_I_D"} else {print $2,$1"_"$4"_"$5"_"$6} }' ${input_fileset}.bim > ${output}.txt
## add "DUPS" with the number of times that value has been seen to the IDs of duplicates
awk 'seen[$2]++{$2=$2"_DUPS_"seen[$2]} 1' ${output}.txt > ${output}_temp && mv ${output}_temp ${output}.txt
## update using plink
plink --bfile ${input_fileset} --update-name ${output}.txt $plink_memory_limit --make-bed --out ${output} > /dev/null

########## initial SNP filters ##########
printf '\n%s\n' "Step 1: Remove SNPs with >5% missingness or with MAF <0.01"
output_last=${output}
output=${output}_geno05_maf01
plink --bfile $output_last --geno 0.05 --maf 0.01 --make-bed --out $output $plink_memory_limit > /dev/null

# document
## get the number of variants dropped for initial variant filters
vargeno=$(  grep "removed due to missing genotype data" ${output}.log | awk '{print $1;}' )
varmaf=$(  grep "removed due to minor allele threshold(s)" ${output}.log | awk '{print $1;}' )
printf '%s\n' "Removed $vargeno SNPs for >5% missingness and $varmaf SNPs for MAF <0.01"
## get the resulting number of samples and variants
variants=$( grep "pass filters and QC" ${output}.log | awk '{print $1;}' )
samples=$( grep "pass filters and QC" ${output}.log | awk '{print $4;}' )
printf "$samples samples
$variants variants\n"

printf "Output file: $output \n\n"

############# initial sample filters #################

###### person missingness ######
printf '%s\n' "Step 2: remove subjects w/ >1% missingness"
output_last=$output
output=${output}_mind01
plink --bfile $output_last --mind 0.01 --make-bed --out $output $plink_memory_limit > /dev/null

# document
## get the number of samples removed for missingness
samplesgeno=$( grep "removed due to missing genotype data" ${output}.log | awk '{print $1;}' )
printf '%s\n' "Removed $samplesgeno people for >1% missingness."
## get the resulting number of samples and variants
variants=$( grep "pass filters and QC" ${output}.log | awk '{print $1;}' )
samples=$( grep "pass filters and QC" ${output}.log | awk '{print $4;}' )
printf "$samples samples
$variants variants\n"
printf "Output file: $output \n"

########## Relatedness ##########
printf "\nStep 3: Calculate relatedness and remove identical individuals (leaving in related individuals)\n"
plink --bfile $output --genome full unbounded nudge --min 0.20 --out ${output}_relatedness $plink_memory_limit > /dev/null

#print out # or related individuals at each threshold
for pi_hat in 0.9 0.5 0.25; do
    n_rel=$(awk -v val="$pi_hat" '{ if (NR!=1 && $10 > val ) print }' ${output}_relatedness.genome | wc -l)
    printf "Pairs above pi-hat of $pi_hat: $n_rel \n"
done

#print out all basically identical pairs which will be removed entirely to allow for quick scanning for intended sample duplicates (re-genotyped samples, the unlikely event of actual identical twins, ect)
if [ "$( awk '{ if(NR==1 || $10 > 0.9 ) print }' ${output}_relatedness.genome | wc -l )" -gt 1 ];
then
    printf "Sample pairs with pi-hat above 0.9 which will be entirely removed:\n"
    awk '{ if(NR==1 || $10 > 0.9 ) print $1" "$2" "$3" "$4" "$10 }' ${output}_relatedness.genome | tee ${output}_relatedness_identical.txt

    # remove selected ids from the last generated genotype file set
    output_last=$output
    output=${output}_noidentical
    plink --bfile $output_last --remove ${output_last}_relatedness_identical.txt --make-bed --out $output $plink_memory_limit > /dev/null

    # document
    ## get the number of samples dropped for being identical
    samplesidentical=$(($samples - $( grep "people remaining" ${output}.log | awk '{print $2;}' )))
    printf "Removed $samplesidentical related individuals.\n"
    ## get the resulting number of samples and variants
    variants=$( grep "pass filters and QC" ${output}.log | awk '{print $1;}' )
    samples=$( grep "pass filters and QC" ${output}.log | awk '{print $4;}' )
    printf "$samples samples
$variants variants\n"
    printf "Output file: $output \n"

else
    printf "No identical samples to remove.\n"
fi

########## sex check ##########
printf "\nStep 4: sex check\n"

num_x_chr=$( awk '{ if($1 == 23) print }' ${output}.bim | wc -l )
if [ "$num_x_chr" -gt 0 ];
then 
    printf "$num_x_chr present to check reported sex against genotypic sex.\n"
    #only update sex if there is something more than missing values for sex in the provided file
    if [ "$( awk '{ print $3 }' $sex_file | sort -u | tr -d '[:space:]' )" != 0 ];
    then 
	output_last=$output
	output=${output}_sex
	plink --bfile $output_last --update-sex $sex_file --make-bed --out $output $plink_memory_limit > /dev/null
	printf "Updated sex from the provided file: $sex_file \n"
    else
	printf "No non-missing sex information in the provided file (${sex_file}). Using sex from fam file.\n"
    fi
    #do the check
    plink --bfile $output --check-sex --out ${output}_checking_sex $plink_memory_limit > /dev/null
    grep -e 'check-sex: ' ${output}_checking_sex.log 

    #write mismatched sex iids in text file
    awk '{ if($5=="PROBLEM" && $4 != 0 && $3 != 0) print $1" "$2 }' ${output}_checking_sex.sexcheck > ${output}_mismatched_sex_ids.txt
    printf "$( wc -l < ${output}_mismatched_sex_ids.txt ) real sex mismatches (e.g. not ambiguous in the fam or indeterminate based on SNPs)\n$( awk '{ if($3 == 0) print }' ${output}_checking_sex.sexcheck | wc -l ) out of $( wc -l < ${output}.fam ) samples are missing sex.\n"

    #remove individuals in text file
    if [ $( wc -l < ${output}_mismatched_sex_ids.txt ) -gt 0 ];
    then
	sex_mismatch_file=${output}_mismatched_sex_ids.txt
	output_last=$output
	output=${output}_nomismatchedsex
	plink --bfile $output_last --remove ${sex_mismatch_file} --make-bed --out $output $plink_memory_limit > /dev/null

	# document
	## get the number of samples which fail sex check
	samplessex=$(($samples - $( grep "people remaining" ${output}.log | awk '{print $2;}' )))
	printf "Removed $sampelssex individuals with mismatched sex.\n"
	## get the resulting number of samples and variants
	variants=$( grep "pass filters and QC" ${output}.log | awk '{print $1;}' )
	samples=$( grep "pass filters and QC" ${output}.log | awk '{print $4;}' )
	printf "$samples samples\n$variants variants\n"
	printf "Output file: $output \n"
    else
	printf "Removed 0 individuals with mismatched sex.\n"
	printf "$samples samples\n$variants variants\n"
    fi
else
    printf "No X chromosomes present to compare against self-report sex. Skipping sex check.\n"
fi

########## restrict to autosomes ##########
printf "\nStep 5 : restrict to autosomes\n"
output_last=$output
output=${output}_keep_autosomes
plink --bfile $output_last --chr 1-22 --make-bed --out $output $plink_memory_limit > /dev/null

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
# document
## get the resulting number of samples and variants
variants=$( grep "pass filters and QC" ${output}.log | awk '{print $1;}' )
samples=$( grep "pass filters and QC" ${output}.log | awk '{print $4;}' )
printf "$samples samples
$variants variants\n"
printf "Output file: $output \n"

########## Heterozygosity check ##########

printf "\nStep 6: Pruning genotypes and running heterozygosity check.\n"

# prune and set phenotypes to 1
plink --bfile ${output} --indep-pairwise 200 100 0.2 --allow-no-sex --out ${output}_prune $plink_memory_limit > /dev/null
plink --bfile ${output} --output-missing-phenotype 1 --extract ${output}_prune.prune.in --make-bed --out ${output}_pruned $plink_memory_limit > /dev/null
rm ${output}_prune.*
printf "$( wc -l < ${output}_pruned.bim ) variants out of $( wc -l < ${output}.bim ) left after pruning.\n"

# heterozygosity check 
plink --bfile ${output}_pruned --het --out ${output}_pruned_hetcheck $plink_memory_limit > /dev/null
# plot and check for outliers
Rscript plot_het_check_outliers.R ${output}_pruned_hetcheck

#if there are outliers >6 sd from the F stat mean, remove them
if [ -f "${output}_pruned_hetcheck_outliers.txt" ];
then
    output_last=${output}
    output=${output}_nohetout
    plink --bfile $output_last --remove ${output}_pruned_hetcheck_outliers.txt --make-bed --out ${output} $plink_memory_limit > /dev/null

    # document
    ## get the number of samples which fail sex check
    samplessex=$(($samples - $( grep "people remaining" ${output}.log | awk '{print $2;}' )))
    printf "Removed $sampelssex individuals with mismatched sex.\n"
    ## get the resulting number of samples and variants
    variants=$( grep "pass filters and QC" ${output}.log | awk '{print $1;}' )
    samples=$( grep "pass filters and QC" ${output}.log | awk '{print $4;}' )
    printf "$samples samples\n$variants variants\n"
    printf "Output file: $output \n\n"
else
    printf "$samples samples\n$variants variants\n"
    printf "Output file: $output \n\n"
fi

########## HWE filter ##########

printf "\nStep 7: Running the Hardy-Weinberg Equilibrium SNP filter \n"  
plink --bfile $output --hwe 0.000001 --make-bed --out ${output}_hwe6 $plink_memory_limit > /dev/null
output=${output}_hwe6

# document
## get the number of variants dropped for initial variant filters
varhwe=$(  grep "removed due to Hardy-Weinberg exact test" ${output}.log | awk '{print $2;}' )
printf "Removed $varhwe SNPs for HWE p<1e-6.\n"
## get the resulting number of samples and variants
variants=$( grep "pass filters and QC" ${output}.log | awk '{print $1;}' )
samples=$( grep "pass filters and QC" ${output}.log | awk '{print $4;}' )
printf "$samples samples
$variants variants\n"
printf "Output file: $output \n"

########## remove palindromic ##########

printf "\nStep 8: Removing palindromic variants \n"

#get palindromic variants
awk '{ if( ($5 == "A" && $6 == "T") || ($5 == "T" && $6 == "A") || ($5 == "C"  && $6 == "G") || ($5 == "G" && $6 == "C")) print }'  ${output}.bim > ${output_path}/palindromic_snps.txt

#remove them
output_last=$output
output=${output}_nopal
plink --bfile $output_last --exclude ${output_path}/palindromic_snps.txt --make-bed --out $output $plink_memory_limit > /dev/null
printf "Removed $( wc -l < ${output_path}/palindromic_snps.txt ) palindromic variants.\n"
varspalindromic=$( wc -l < ${output_path}/palindromic_snps.txt )
printf "Output file: $output \n"

########## liftOver ##########

#check the build, decide whether to do the liftOver and, if so, which chain file to use
if [ "$build" = "b37" ];
then
    printf "\nStep 9: Lifting over the genotype file from build 37 to build 38 to match the TOPMed reference panel\n"
    chain_file="hg19ToHg38.over.chain.gz"
elif [ "$build" = "b36" ];
then 
    printf "\nStep 9: Lifting over the genotype file from build 36 to build 38 to match the TOPMed reference panel\n"
    chain_file="hg18ToHg38.over.chain.gz"
elif [ "$build" = "b35" ];
then 
    printf "\nStep 9: Lifting over the genotype file from build 35 to build 38 to match the TOPMed reference panel\n"
    chain_file="hg17ToHg38.over.chain.gz"
elif [ "$build" = "b34" ];
then 
    printf "\nStep 9: Lifting over the genotype file from build 34 to build 38 to match the TOPMed reference panel\n"
    chain_file="hg16ToHg38.over.chain.gz"
elif [ "$build" = "b38" ];
then
    printf "\nStep 9: Input data was specified as already on build 38, so no lift-over is necessary. Proceeding to the next step.\n"
fi

#if the chain file variable is set, then run the liftOver
if [ ! -z "$chain_file" ];
then
    #create bed file to start the liftover
    awk '{ print "chr"$1" "$4 -1" "$4" "$2 }' ${output}.bim > ${output}_toliftover.txt
    #lift over
    ./liftOver ${output}_toliftover.txt $chain_file ${output}_lifted.txt ${output}_unlifted.txt
    #remove "chr" from chromosome column
    sed -i 's/^chr//g' ${output}_lifted.txt

    #get list of variants to remove (those with non-standard chr codes or which cannot be lifted)
    awk 'length($1)>2 { print $4 }' ${output}_lifted.txt > ${output}_lift_issue_snps_todrop.txt
    varsnonstd=$( wc -l < ${output}_lift_issue_snps_todrop.txt )
    awk '{ print $4 }' ${output}_unlifted.txt | grep -v '^$' | sort -u >> ${output}_lift_issue_snps_todrop.txt
    varsnotpresent=$( awk '{ print $4 }' ${output}_unlifted.txt | grep -v '^$' | sort -u | wc -l )
    varsfaillift=$( wc -l < ${output}_lift_issue_snps_todrop.txt )
    printf "Removing $varsfaillift variants which fail liftOver ($varsnonstd 'random', 'alt', 'fix', or other, and $varsnotpresent not present in b38).\n"

    #drop 
    plink --bfile $output --exclude ${output}_lift_issue_snps_todrop.txt --make-bed --out ${output}_noprobSNPs $plink_memory_limit > /dev/null
    #update the chr and BP positions for remaining
    plink --bfile ${output}_noprobSNPs --update-chr ${output}_lifted.txt 1 4 --make-bed --out ${output}_noprobSNPs_chr $plink_memory_limit > /dev/null
    plink --bfile ${output}_noprobSNPs_chr --update-map ${output}_lifted.txt 3 4 --make-bed --out ${output}_noprobSNPs_chr_bplifted $plink_memory_limit > /dev/null
    output=${output}_noprobSNPs_chr_bplifted
fi

########## remove same-position variants ##########

printf "\nStep 10: Removing multi-allelic and duplicated variants.\n"

## first remove actual duplicated variants using plink
## the --list-duplicate-vars with the suppress-first modifier checks for duplicates based on position and allele (not ID)
output_last=$output
output=${output}_nodups
plink --bfile $output_last --list-duplicate-vars ids-only suppress-first --out ${output_last}_dups $plink_memory_limit > /dev/null
dupvars=$( wc -l < ${output_last}_dups.dupvar )
printf "Removing $dupvars duplicate variants.\n"
plink --bfile $output_last --exclude ${output_last}_dups.dupvar --make-bed --out $output $plink_memory_limit > /dev/null
printf "Output file: $output \n"

## now get variants that are at the same position. This will be only multi-allelic variants now. 
awk '{ print $2" "$1"_"$4 }' ${output}.bim | sort -T ${output_path}/ -k2 | uniq -f1 -D | awk '{ print $1 }' > ${output}_samepos_vars.txt
## if there are any present, remove them
if [ "$( wc -l < ${output}_samepos_vars.txt )" -gt 0 ];
then
    varssamepos=$( wc -l < ${output}_samepos_vars.txt )
    printf "Removing $varssamepos same-position variants.\n"
    plink --bfile ${output} --exclude ${output}_samepos_vars.txt --make-bed --out ${output}_nosamepos $plink_memory_limit > /dev/null
    output=${output}_nosamepos
    printf "Output file: $output\n"
    # add up for documentation
    varssamepos=$(($dupvars+$varssamepos))
else 
    printf "No multi-allelic variants to remove.\n"
    varssamepos=$dupvars
fi

########## Imputation prep ##########

printf "\nStep 11: Comparing with the reference panel and preparing files for imputation for each chromosome.\n"

#make set with short name and freq file
output_last=${output}
output=${output_stem}_forimputation
plink --bfile ${output_last} --allow-no-sex --freq --make-bed --out ${output} $plink_memory_limit > /dev/null

#run the imputation checking script
perl HRC-1000G-check-bim.pl -b ${output}.bim  -f ${output}.frq -r ${ref_file_stem}.txt.gz -h -n > /dev/null

#if noexclude, then skip the exclusion step in the plink file
if [ "$noexclude" = true ];
then
    sed -i -e '1s/.*/#&/' -e "s|${output_path}/TEMP1|${output}|g" ${output_path}/Run-plink.sh
fi

#run created script with all the plink commands
# muting all output since this will throw quite a few warnings if there are any funky alleles
sed -i "s|rm TEMP|rm ${output_path}/TEMP|" ${output_path}/Run-plink.sh
sh ${output_path}/Run-plink.sh > /dev/null 2>&1

#run imputation prep for each chromosome
for i in $( seq 1 22 ) ;
do
    printf "Beginning chr ${i}...\n"
    #update chr code to have chr# rather than just #, sort and output vcf
    mkdir ${output_path}/tmp${i}/
    bcftools annotate --rename-chrs update_chr_names_b38.txt ${output}-updated-chr${i}.vcf | bcftools sort --temp-dir ${output_path}/tmp${i}/ -O z -o ${output}-updated-chr${i}.vcf.gz > /dev/null
    printf "\nChr ${i} complete...\n"
done

# document (format for the workflow)
## calculate variants removed bc of mismatch with reference
varsmismatchref=$(($( wc -l < ${output}.bim )-$( wc -l < ${output}-updated.bim )))
printf "\nSummary:\nRemoved $varspalindromic palindromic, $varsfaillift failing liftOver to b38 , $varssamepos same position SNPs, and $varsmismatchref SNPs not matching or not present in reference panel.\n"
## get the resulting number of samples and variants
variants=$( grep "pass filters and QC" ${output}*-updated.log | awk '{print $1;}' )
samples=$( grep "pass filters and QC" ${output}*-updated.log | awk '{print $4;}' )
printf "$samples samples
$variants variants\n"

## commenting this out for now to prevent accidental deletion of files (since the pipeline has been updated)
##remove the intermediate .vcf and .bed files
###make sure that it only removes the files from running this script, not the initial files or files for other datsets
#if [ "$skip_cleanup" = false ];
#then
#    rm ${output_path}/*.vcf
#    files_to_remove=$( find ${output_path}/${input_stem}_hwe6*.bed | grep -v "${output}.bed" | grep -v "${output}-updated.bed" )
#    rm $files_to_remove
#else
#    printf "Skipping cleanup. Please manually remove unnecessary files.\n"
#fi

printf "\nConversion complete! Upload the files (${output}-updated-chr*.vcf.gz) to the imputation server.\n"

# # get rid of all the bim files except the last one
# ## make sure that it's only removing files from this set in case others are being run in the same folder
# printf "PC plots complete. Doing some clean-up...\n"
# files_to_delete=$( find ${output_stem}*.bed | grep -v "${output}.bed" )
# rm $files_to_delete
