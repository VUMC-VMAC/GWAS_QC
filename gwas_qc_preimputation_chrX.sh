#!/bin/bash
# Author: Vaibhav Janve 2020-04-07
# Author 2: Emily Mahoney 2020-04-09
# Aug 2021, Vaibhav Janve: script update to process X-chr
# July 2024, ERM: SNPWeights update

#fail on error
set -e

#define usage
display_usage() {
    printf "GWAS QC for X Chromosome Preimputation 

Completes the first stage (including X-chromosome) in standard GWAS QC, including initial variant and person filters, relatedness and sex checks, restriction to autosomes, and PC calculation. The phenotype column in .fam is updated with sex in process.

Usage:
gwas_qc_part1x.sh -o [output_stem] -i [input_fileset] -f [sex_file] -b [b37] -m [plink_memory_limit] -R [ref_file]

output_stem = the beginning part of all QC'ed files, including the full file path to the directory where the files are to be saved

input_fileset = the full path and file stem for the raw plink set '*[bed,bim,fam]'

sex_file = a file with FID and IID (corresponding to the fam file), 1 column indicating sex for the sex check (1 for males, 2 for females, 0 if unknown), with NO header.

ref_file = the full file path and name for the reference panel

-b = optional argument indicating the build of the input dataset; default assumption is b37 but can supply b36 or b38. This will impact the x-chromosome split step.

plink_memory_limit (optional) = argument indicating the memory limit for plink to use rather than the default of half the RAM. This is useful for running this step of QC locally.

-h will show this usage
"
        }

build=b37
while getopts 'o:i:f:b:R:m:h' flag; do
  case "${flag}" in
    o) output_stem="${OPTARG}" ;;
    i) input_fileset="${OPTARG}" ;;
    f) sex_file="${OPTARG}" ;;
    b) build="${OPTARG}" ;;
    m) plink_memory_limit="${OPTARG}" ;;
    R) ref_file="${OPTARG}" ;;
    h) display_usage ; exit ;;
    \?|*) display_usage
       exit 1;;
  esac
done

#check to make sure necessary arguments are present
if [ -z "$output_stem" ] || [ -z "$input_fileset" ] || [ -z "$sex_file" ] || [ -z "$ref_file" ] ;
then
    printf "Error: Necessary arguments not present!\n\n"
    display_usage
    exit 1
fi

#print out inputs
printf "GWAS QC for X Chromosome Preimputation\n" 

# Print out singularity information for reproducability
[ ! -z "$SINGULARITY_CONTAINER" ] && printf "\nSingularity image: $SINGULARITY_CONTAINER"
[ ! -z "$SINGULARITY_COMMAND" ] && printf "\nSingularity shell: $SINGULARITY_COMMAND"
[ ! -z "$SINGULARITY_BIND" ] && printf "\nMapped directories:\n $SINGULARITY_BIND\n" | sed 's/,/\n /g'

printf "\nInput data : $input_fileset
Output stem : $output_stem
File with sex information : $sex_file
Reference panel SNP file : $ref_file
"

# check the build argument
if [ "$build" != "b36" ] && [ "$build" != "b37" ] && [ "$build" != "b38" ]; 
then 
    printf "Error: Invalid build argument supplied ($build). Please supply one of the following: b36, b37, b38.\n"
    exit 1
fi

# get X-chromosome SNP count
n_x=$( awk '$1==23{print}' ${input_fileset}.bim | wc -l )

if [ ${n_x} -lt 300 ];
then
    # updated this to simply fail and exit since there is no point in attempting QC with fewer than this number of variants
    printf "%s\n" "ERROR: too few X-chromosomes to perform QC";
    exit 1
else
    # print out the number of X chromosome variants for QC
    printf "Number of X-chromosome SNPs: ${n_x}\n"
fi

# print message for input build
if [ -z ${build} ]
then
   printf "Input data build not specified assuming 'b37' \n"
   build="b37"
else
   printf "Input data build: ${build} \n"
fi

# report plink memory limit
if [ "$plink_memory_limit" ];
then 
    printf "Memory limit for plink calls: $plink_memory_limit \n"
    plink_memory_limit=$( echo "--memory $plink_memory_limit" )
fi


############### add check to make sure there is self-report sex present. If this is not present we cannot validate the correspondence of genotypes to phenotypes, so it makes sense to make this a requirement (although this is not a requirement for autosomal QC). 

#check to make sure this is being run in the scripts folder (checking if necessary script is present)
if test ! -f get_related_ids.R ; then
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


#get output folder
output_folder=${output_stem%/*}
input_plinkset=${input_fileset##*/}

###################################### Initial variant filters ##################################################3

###### format plinkset: update varIDs #####
printf "\nStep 0: update varIDs to standard CHR_POS_REF_ALT format or CHR_POS_I_D for indels and remove duplicates.\n"
output=${output_folder}/${input_plinkset}_stdnames
## generate an ID with just variant info rather than rsID and shorten IDs that are too long
awk '{if( length($5)>25 || length($6)>25) { print $2,$1"_"$4"_I_D"} else {print $2,$1"_"$4"_"$5"_"$6} }' ${input_fileset}.bim > ${output}.txt
## add "DUPS" with the number of times that value has been seen to the IDs of duplicates
awk 'seen[$2]++{$2=$2"_DUPS_"seen[$2]} 1' ${output}.txt > ${output}_temp && mv ${output}_temp ${output}.txt
## update using plink
plink --bfile ${input_fileset} --update-name ${output}.txt $plink_memory_limit --make-bed --out ${output}  > /dev/null

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

###################################### Initial sample filters ##################################################3

##### person missingness #####
printf '%s\n' "Step 2: remove subjects w/ >1% missingness"
output_last=$output
output=${output}_mind01
plink --bfile $output_last --mind 0.01 --make-bed --out $output $plink_memory_limit > /dev/null

# document                                                                                                                                         
## get the number of samples removed for missingness                                                                                               
samplesgeno=$( grep "removed due to missing genotype data" ${output}.log | awk '{print $1;}' )
printf '%s\n' "Removed $samplesgeno people for >1% missingness.\n"
## get the resulting number of samples and variants                                                                                                
variants=$( grep "pass filters and QC" ${output}.log | awk '{print $1;}' )
samples=$( grep "pass filters and QC" ${output}.log | awk '{print $4;}' )
printf "$samples samples                                                                                                                           
$variants variants\n"
printf "Output file: $output \n"


##### Relatedness #####
printf "\nStep 3: Calculate relatedness and remove related individuals\n"
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
    awk '{ if(NR==1 || $10 > 0.9 ) print $1" "$2" "$3" "$4" "$10 }' ${output}_relatedness.genome
fi

# Rscript for decisions
Rscript get_related_ids.R ${output}_relatedness

# remove selected ids from the last generated genotype file set
output_last=$output
output=${output}_norelated
plink --bfile $output_last --remove ${output_last}_relatedness_related_ids.txt --make-bed --out $output $plink_memory_limit > /dev/null

# document                                                                                                                                         
## get the number of samples dropped for relatedness                                                                                               
samplesrelated=$(($samples - $( grep "people remaining" ${output}.log | awk '{print $2;}' )))
printf "Removed $samplesrelated related individuals.\n"
## get the resulting number of samples and variants                                                                                                
variants=$( grep "pass filters and QC" ${output}.log | awk '{print $1;}' )
samples=$( grep "pass filters and QC" ${output}.log | awk '{print $4;}' )
printf "$samples samples                                                                                                                           
$variants variants\n"
printf "Output file: $output \n"

##### sex check #####
printf "\nStep 4: sex check\n"

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

# check to confirm that sex information is present and fail if not
if [[ $( grep "ambiguous" ${output}.fam ) -ne 0 ]];
then
    printf "Error: Missing sex for $( awk '{ if($5==0) print }' /data/h_vmac/waltes2/X_QC/WHICAP/1099_GWAS/PrePart1/1099GWAS_part1X_nopal.fam | wc -l ) samples. Please check your input files! Self-report sex is required to confirm data quality and to perform sex-specific variant filters on X chromosome variants.\n"
    exit 1
fi

#do the check
plink --bfile $output --check-sex --out ${output}_checking_sex $plink_memory_limit > /dev/null
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


plinkset_x_out=${output}
# confirm sex is updated
n_sex=$(awk '{print $5}' ${plinkset_x_out}.fam | sort | uniq -dc | wc -l)
sex=$(awk '{print $5}' ${plinkset_x_out}.fam | sort | uniq -dc | awk '{print $2}')

if [ $n_sex -ne 2 ];
then
    # printing this out for information in case there are sets which have only males/females
    printf "\nWarning: more or less than 2 sexes found in fam:
  ${plinkset_x_out}.fam , 1=Male, 2=Female;
  $(awk '{print $5}' ${plinkset_x_out}.fam | sort | uniq -dc )\n"
fi

###################################### Begin X chr specific processing ##################################################3

printf "\nBeginning X-chromosome processing...\n"

##### save plinkset for x-chromosome processing #####
plinkset_x_in=${plinkset_x_out}
plinkset_x_first=${plinkset_x_in}

#### STEP 5 Subset to X chromosome ####

# merge and then split chrX with position range based on build
printf "\nStep 5: Subset to X chromosome, removign autosomal and PAR variants.\n"
printf "\nFirst, assign pseudoautosomal regions (PARS) based on build ${build} coordinates.\n"
plinkset_x_in=${plinkset_x_out}
plinkset_x_out=${plinkset_x_in}_mergeX
plink --bfile ${plinkset_x_in} --merge-x no-fail --make-bed --out ${plinkset_x_out} $plink_memory_limit > /dev/null
grep -e 'merge-x: ' ${plinkset_x_out}.log ||
printf " --merge-x: 0 chromosome codes changed.\n no chr 25 found. Please confirm %s\n" ${plinkset_x_out}.log  >&2

printf "Now, split PAR region from X chr.\n"
plinkset_x_in=${plinkset_x_out}
plinkset_x_out=${plinkset_x_in}_splitX
plink --bfile ${plinkset_x_in} --split-x ${build} no-fail --make-bed --out ${plinkset_x_out} $plink_memory_limit > /dev/null
grep -e 'split-x: ' ${plinkset_x_out}.log ||
printf " --split-x: 0 chromosome codes changed.\n no chr 25 found. Please confirm %s\n" ${plinkset_x_out}.log  >&2

printf '%s\n' "Now, setting het. haploid genotypes missing."
# set hh missing to deal with het haplotype warning (since after sex check)
# check no hh warnings on females and set hh to missing
if [  $n_sex -eq 1 ] && [ ${sex} -eq 1 ];
then
    printf "\nOnly males or other sex present \n"
else
    printf "\nConfirm no hh warning for females only selection below: 'Warning: ... het. haploid genotypes present' \n"
    plink --bfile ${plinkset_x_out} --filter-females --freq $plink_memory_limit --out ${plinkset_x_out} > /dev/null
fi

# set hh missing
plinkset_x_in=${plinkset_x_out}
plinkset_x_out=${plinkset_x_out}_nohh
printf "\nSet het. haploid genotypes missing in ${plinkset_x_out}\n"
plink --bfile ${plinkset_x_in} --set-hh-missing --make-bed --out ${plinkset_x_out} $plink_memory_limit > /dev/null

## now that all other variants have been split or removed, subset to chr 23 variants only
printf "Finally, subsetting to X-chromosome variants only.\n"
plinkset_x_in=${plinkset_x_out}
plinkset_x_out=${plinkset_x_first}_X
#subset to X-chr SNPS
plink --bfile ${plinkset_x_in} --chr 23 --make-bed --out ${plinkset_x_out} $plink_memory_limit > /dev/null    
## log the number of variants removed 
n_x2=$( wc -l < ${plinkset_x_out}.bim )
printf "Removed $(( ${n_x} - ${n_x2})) PARS SNPs.\n"
## get the resulting number of samples and variants                                                                                                
samples=$( grep "pass filters and QC" ${plinkset_x_out}.log | awk '{print $4;}' )
printf "$samples samples
$n_x2 variants\n"
printf "Output file: $plinkset_x_out \n"

#### STEP 6 Performing differential missingness filter ####

if [  $n_sex -eq 1 ];
then
    printf "\nSkipping Step 6: Diffential missingness calculations because only one sex is present.\n"
else
    printf "\nStep 6: Performing differential missingness calculations.\n"
    
    # update pheno with sex
    printf "First, updating phenotype in plink files to correspond to sex for ease of missingness calculations...\n"
    plinkset_x_in=${plinkset_x_out}
    plinkset_x_out=${plinkset_x_out}_sexpheno
    awk -F ' ' '{print $1,$2,$5}' ${plinkset_x_in}.fam > ${plinkset_x_out}.txt
    plink --bfile ${plinkset_x_in} --make-pheno ${plinkset_x_out}.txt 2 --make-bed --out ${plinkset_x_out} $plink_memory_limit > /dev/null

    # generate male-female diff miss stats
    printf "Now, calculating missingness statistics by sex.\n"
    plink --bfile ${plinkset_x_out} --test-missing --out ${plinkset_x_out}_diffmiss $plink_memory_limit > /dev/null
    # drop SNPs with P < 0.0000001
    ## generate list of variants
    awk '$5 <= 0.0000001{print $2}' ${plinkset_x_out}_diffmiss.missing > ${plinkset_x_out}_diffmiss_to_drop.txt

    ## remove using plink
    printf "Finally, removing using plink.\n"
    plinkset_x_in=${plinkset_x_out}
    plinkset_x_out=${plinkset_x_in}_e7diffmis
    plink --bfile ${plinkset_x_in} --exclude ${plinkset_x_in}_diffmiss_to_drop.txt --make-bed --out ${plinkset_x_out} $plink_memory_limit > /dev/null

    # document
    printf "\nRemoved $(wc -l < ${plinkset_x_in}_diffmiss_to_drop.txt) SNPs with male/female differential missingness.\n"
    ## get the resulting number of samples and variants                             
    variants=$( grep "pass filters and QC" ${plinkset_x_out}.log | awk '{print $1;}' )
    samples=$( grep "pass filters and QC" ${plinkset_x_out}.log | awk '{print $4;}' )
    printf "$samples samples\n$variants variants\n"
    printf "Output file: $plinkset_x_out \n"

fi


##### Step 7: Hardy Weinberg equilibrium 1e6 (based on females) #####
if [  $n_sex -eq 1 ] && [ ${sex} -eq 1 ];
then
    printf "\nSkipping step 7: Hardy Weinberg equilibrium 1e6 (based on females) because only makes are present.\n"
else
    printf "\nStep 7: Hardy Weinberg equilibrium 1e6 (based on females)\n"
    plinkset_x_in=${plinkset_x_out}
    plinkset_x_out_fem=${plinkset_x_in}_1e6femHWE
    plinkset_x_out=${plinkset_x_in}_femHWE6
    # get SNPS out of HWE in females
    plink --bfile ${plinkset_x_in} --filter-females --hwe 0.000001 --make-just-bim --out ${plinkset_x_out_fem} $plink_memory_limit > /dev/null
    awk '{print $2}' ${plinkset_x_out_fem}.bim > ${plinkset_x_out_fem}_SNPs.txt
    plink --bfile ${plinkset_x_in} --extract ${plinkset_x_out_fem}_SNPs.txt --make-bed --out ${plinkset_x_out} $plink_memory_limit > /dev/null

    # document
    printf "\nRemoved $(wc -l < ${plinkset_x_out_fem}_SNPs.txt) SNPs for HWE p<1e-6 in females.\n"
    ## get the resulting number of samples and variants                             
    variants=$( grep "pass filters and QC" ${plinkset_x_out}.log | awk '{print $1;}' )
    samples=$( grep "pass filters and QC" ${plinkset_x_out}.log | awk '{print $4;}' )
    printf "$samples samples\n$variants variants\n"
    printf "Output file: $plinkset_x_out \n"

fi


##################################### Begin preparation for imputation #####################################3

printf "\nStep 8: Removing palindromic variants \n"

#get palindromic variants
awk '{ if( ($5 == "A" && $6 == "T") || ($5 == "T" && $6 == "A") || ($5 == "C"  && $6 == "G") || ($5 == "G" && $6 == "C")) print }'  ${plinkset_x_out}.bim > ${output_folder}/palindromic_snps.txt
#remove them
output_last=$plinkset_x_out
output=${plinkset_x_out}_nopal
plink --bfile $output_last --exclude ${output_folder}/palindromic_snps.txt --output-chr M --make-bed --out $output > /dev/null
varspalindromic=$( wc -l < ${output_folder}/palindromic_snps.txt )
printf "Removed $varspalindromic palindromic variants.\n"
printf "Output file: $output \n"

#check the build, decide whether to do the liftOver and, if so, which chain file to use
if [ "$build" = "b37" ];
then
    printf "\nStep 9: Lifting over the genotype file from build 37 to build 38 to match the TOPMed reference panel\n"
    chain_file="hg19ToHg38.over.chain.gz"
elif [ "$build" = "b36" ];
then 
    printf "\nStep 9: Lifting over the genotype file from build 36 to build 38 to match the TOPMed reference panel\n"
    chain_file="hg18ToHg38.over.chain.gz"
elif [ "$build" = "b38" ];
then
    printf "\nSkipping step 9: Input data was specified as already on build 38, so no lift-over is necessary. Proceeding to the next step.\n"
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

    #get list of variants to remove (those with non-standard allele codes or which cannot be lifted)
    grep -e "random" -e "alt" -e "fix" ${output}_lifted.txt | awk '{ print $4 }' > ${output}_lift_issue_snps_todrop.txt
    awk '{ print $4 }' ${output}_unlifted.txt | sort -u >> ${output}_lift_issue_snps_todrop.txt
    varsfaillift=$( wc -l < ${output}_lift_issue_snps_todrop.txt )
    printf "Removing $varsfaillift variants which fail liftOver ($( grep "random" ${output}_lifted.txt | wc -l ) 'random', $( grep "alt" ${output}_lifted.txt | wc -l ) 'alt', $( grep "fix" ${output}_lifted.txt | wc -l ) 'fix', and $( awk '{ print $4 }' ${output}_unlifted.txt | sort -u | wc -l ) not present in b38).\n"

    #drop 
    plink --bfile $output --exclude ${output}_lift_issue_snps_todrop.txt --output-chr M --make-bed --out ${output}_noprobSNPs > /dev/null
    #update the chr and BP positions for remaining
    plink --bfile ${output}_noprobSNPs --update-chr ${output}_lifted.txt 1 4 --output-chr M --make-bed --out ${output}_noprobSNPs_chr > /dev/null
    plink --bfile ${output}_noprobSNPs_chr --update-map ${output}_lifted.txt 3 4 --make-bed --out ${output}_noprobSNPs_chr_bplifted > /dev/null
    output=${output}_noprobSNPs_chr_bplifted
fi

#### remove same-position variants ####

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
awk '{ print $2" "$1"_"$4 }' ${output}.bim | sort -T ${output_folder}/ -k2 | uniq -f1 -D | awk '{ print $1 }' > ${output}_samepos_vars.txt
## if there are any present, remove them
varssamepos=$( wc -l < ${output}_samepos_vars.txt )
if [ "$varssamepos" -gt 0 ];
then
    printf "Removing $varssamepos same-position variants.\n"
    plink --bfile ${output} --exclude ${output}_samepos_vars.txt --make-bed --out ${output}_nosamepos > /dev/null
    output=${output}_nosamepos
    printf "Output file: $output\n"
    # add up for documentation                                                                                                                     
    varssamepos=$(($dupvars+$varssamepos))
else 
    printf "No multi-allelic variants to remove.\n"
    varssamepos=$dupvars
fi

#### Imputation prep ####

printf "\nStep 11: Comparing with the reference panel and preparing files for imputation for each chromosome.\n"

#make set with short name and freq file
output_last=${output}
output=${output_stem}
plink --bfile ${output_last} --freq --make-bed --out ${output} > /dev/null

#run the imputation checking script
perl HRC-1000G-check-bim.pl -b ${output}.bim  -f ${output}.frq -r ${ref_file} -h -n > /dev/null

#if noexclude, then skip the exclusion step in the plink file
if [ "$noexclude" = true ];
then
    sed -i -e '1s/.*/#&/' -e "s|${output_folder}/TEMP1|${output}|g" ${output_folder}/Run-plink.sh
fi

#run created script with all the plink commands
#muting all output since this will throw quite a few warnings if there are any funky alleles
sed -i "s|rm TEMP|rm ${output_folder}/TEMP|" ${output_folder}/Run-plink.sh
sh ${output_folder}/Run-plink.sh > /dev/null 2>&1

## calculate variants removed bc of mismatch with reference
varsmismatchref=$(($( wc -l < ${output}.bim )-$( wc -l < ${output}-updated-chr23.bim )))

printf "\n update variant IDs to TOPMED imputation server format chrX:POS:REF:ALT \n"
# the genotyped variant ids in TOPMED imputation server format chrX:POS:REF:ALT
output_last=${output}-updated-chr23
output=${output_last}_TOPMED_varID
awk '{ print "chrX:"$4":"$6":"$5,$2}' ${output_last}.bim > ${output_last}_TOPMED_varID.txt

plink --bfile ${output_last} --update-name ${output_last}_TOPMED_varID.txt 1 2 --make-bed $plink_memory_limit --out ${output}  > /dev/null

# create VCF
## there was another filter to just chr 23 here as well as the --real-ref-alleles flag which has already been applied above and nothing would have changed the ref alleles
## since both were superfluous they were removed
plink --bfile ${output} --recode-vcf --out ${output} $plink_memory_limit > /dev/null

printf "Output: ${output} \n"
grep -e ' people pass filters and QC' ${output}.log

#run imputation prep
#update chr code to have chr# rather than just #, sort and output vcf
printf "Sorting vcf and updating chromosome code for imputation server...\n"
mkdir ${output_folder}/tmpX/
bcftools annotate --rename-chrs update_chr_names_b38.txt ${output}.vcf | bcftools sort --temp-dir  ${output_folder}/tmpX/ -O z -o ${output}_chrX.vcf.gz > /dev/null
printf "Sorting complete.\n\n"

##print out number of variants actually excluded or which would have been excluded
#if [ "$noexclude" = true ];
#then 
#    exclude_var=$( cat ${output_folder}/Exclude-* | wc -l )
#    current_var=$( cat ${output}.bim | wc -l )
#    hypothetical_var=$(( $excluded_var - $current_var ))
#    printf "Would have removed ${excluded_var} variants for mismatch with the reference panel, being palindromic with MAF > 0.4, or being absent from the reference panel leaving ${hypothetical_var} for imputation, but the no-exclude option was specified.\n"
#else
#    printf "Removed $( cat ${output_folder}/Exclude-* | wc -l ) variants for mismatch with the reference panel, being palindromic with MAF > 0.4, or being absent from the reference panel leaving $( cat ${output}.bim | wc -l ) for imputation.\n"
#fi

# document (format for the workflow)
printf "\nSummary:\nRemoved $varspalindromic palindromic, $varsfaillift failing liftOver to b38 , $varssamepos same position SNPs, and $varsmismatchref SNPs not matching or not present in reference panel.\n"
## get the resulting number of samples and variants                                                                                                
variants=$( grep "pass filters and QC" ${output}*.log | awk '{print $1;}' )
samples=$( grep "pass filters and QC" ${output}*.log | awk '{print $4;}' )
printf "$samples samples                                                                                                                           
$variants variants\n"

# commenting this out to prevent removal of files that ought to be kept
#remove the intermediate .vcf and .bed files
##make sure that it only removes the files from running this script, not the initial files or files for other datsets
#if [ "$skip_cleanup" = false ];
#then
#    rm ${output_folder}/*.vcf
#    files_to_remove=$( find ${output_folder}/${input_stem}*.bed | grep -v "${output}.bed" | grep -v "${output_stem}-updated.bed" )
#    rm $files_to_remove
#else
#    printf "Skipping cleanup. Please manually remove unnecessary files.\n"
#fi

printf "\nConversion complete! Upload the file (${output_stem}_chrX.vcf.gz) to the imputation server.\n"
