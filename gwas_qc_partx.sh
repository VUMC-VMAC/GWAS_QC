#!/bin/bash
# Author: Vaibhav Janve 2020-04-07
# Author 2: Emily Mahoney 2020-04-09
# Aug 2021, Vaibhav Janve: script update to process X-chr
# July 2024, ERM: SNPWeights update

# VJ: singularity image and directories bound for reproducibility
[ ! -z "$SINGULARITY_CONTAINER" ] && printf "\nSingularity image: $SINGULARITY_CONTAINER"
[ ! -z "$SINGULARITY_BIND" ] && printf "\nMapped directories:\n $SINGULARITY_BIND\n" | sed 's/,/\n /g'

#fail on error
set -e

#define usage
display_usage() {
    printf "GWAS QC for X Chromosome Preimputation 

Completes the first stage (including X-chromosome) in standard GWAS QC, including initial variant and person filters, relatedness and sex checks, restriction to autosomes, and PC calculation. The phenotype column in .fam is updated with sex in process.

Usage:
gwas_qc_part1x.sh -o [output_stem] -i [input_fileset] -f [sex_file] -b [b37] -p [ds_plinkset] -m [plink_memory_limit] -s -R [ref_file]

output_stem = the beginning part of all QC'ed files, including the full file path to the directory where the files are to be saved

input_fileset = the full path and file stem for the raw plink set '*[bed,bim,fam]'

sex_file = a file with FID and IID (corresponding to the fam file), 1 column indicating sex for the sex check (1 for males, 2 for females, 0 if unknown), with NO header.

ref_file = the full file path and name for the reference panel

-p = autosomal plinkset stem with IID and FID matching the input plinkset and PC outliers removed

-b = optional argument indicating the build of the input dataset; default assumption is b37 but can supply b36 or b38. This will impact the x-chromosome split step.

plink_memory_limit (optional) = argument indicating the memory limit for plink to use rather than the default of half the RAM. This is useful for running this step of QC locally.

-s optional, to skip dup removal step after initial id format
-h will show this usage
"
        }

skip_flag='false'
while getopts 'o:i:f:b:p:R:sh' flag; do
  case "${flag}" in
    o) output_stem="${OPTARG}" ;;
    i) input_fileset="${OPTARG}" ;;
    f) sex_file="${OPTARG}" ;;
    p) ds_plinkset="${OPTARG}" ;;
    b) build="${OPTARG}" ;;
    R) ref_file="${OPTARG}" ;;
    s) skip_flag='true' ;;
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

printf "\n%s\n\n" "Brief summary:
Step0: format plinkset
	update varIDs to standard CHR:POS:REF:ALT format and remove duplicates(keeps unique, only duplicates removed)
Step1: initial SNP filters
        Remove SNPs with >5% missingness or with MAF <0.01 (--geno 0.05 --maf 0.01)
Step2: person missingness
        remove individuals w/ >1% missingness (--mind 0.01)
Step3: Relatedness
        removes 1 of each pair of second degree relatives, 0.9> pi-hat >0.25, and both of each pair with a pi-hat >=0.9
Step4: Sex check
        update sex in .fam with sex_file
        removes individuals with discrepancies
########## X-chromosome Processing ##########
Step1: Assign pseudoautosomal regions (PARS) based on build
        merge and split the genotype on build ${build} coordinates preparing to exclude PARs
Step2: Set het. haploid genotypes missing
        confirm no het. haplotype warnings on female only selection
Step3: Subset to X-chromosome
Step4: *** SKIPPED *** Filter out SNPs with absolute MAF sex difference > 0.02 (2%)
        confirm sex is updated
        new fam with updated pheno with sex
        filter out SNPs with absolute male-female differential MAF: delta > 0.02
Step5:  [if both sexes present] Filter out SNPs on differential missingness
        drop SNPs with P < 0.0000001 (1e-7)
Step6: *** SKIPPED ***(performed in autosomal part1) Heterozygocity + Autosmal PC outliers
Step7: *** SKIPPED ***(we no longer remove samples bc of PCs in autosomal QC at this stage and never have heterozygosity outliers) Remove individuals based on autosomal PCs and heterozygocity outliers (and for restrict to NHW individuals)
Step8: [if females are present] Hardy Weinberg equilibrium(HWE) 1e6 (based on females)
        drop SNPS in both males and females that are out of HWE."


#print out inputs
printf "GWAS QC for X Chromosome Preimputation 

Input data : $input_fileset
Output stem : $output_stem
File with sex information : $sex_file
Reference panel SNP file : $ref_file
"

# get X-chromosome SNP count
n_x=$(awk '$1==23{print}' ${input_fileset}.bim | wc -l)

if [[ ${n_x} -lt 300 ]]; 
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
   printf " Input data build not specified assuming 'b37' \n"
   build="b37"
else
   printf " Input data build: ${build} \n"
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

##### format plinkset: update varIDs #####
printf '%s\n\n' "Step 0: update varIDs to standard CHR:POS:REF:ALT format or CHR:POS:I:D for indels and remove duplicates"
#varID -> CHR:POS:REF:ALT
plinkset_in=${input_fileset}
plinkset_out=${output_folder}/${input_plinkset}_varID
awk '{print $1,$1"_"$4"_"$5"_"$6,$3,$4,$5,$6}' ${input_fileset}.bim > ${plinkset_out}_temp.bim
plink --bed ${plinkset_in}.bed --fam ${plinkset_in}.fam --bim ${plinkset_out}_temp.bim $plink_memory_limit --make-bed --out ${plinkset_out}  > /dev/null

printf " Input: ${plinkset_in}\n"
grep 'loaded from' ${plinkset_out}.log  | sed  's/loaded from.*//' | head -n2
printf " Output: ${plinkset_out} \n\n"

# long indel: smart pca limit 39 characters so keeping limit to 25
plinkset_in=${plinkset_out}
plinkset_out=${plinkset_in}_indel
awk '{if( length($5)>25 || length($6)>25) { print $1,$1"_"$4"_I_D",$3,$4,$5,$6} else {print $0} }' ${plinkset_in}.bim > ${plinkset_out}_temp.bim
plink --bed ${plinkset_in}.bed --fam ${plinkset_in}.fam --bim ${plinkset_out}_temp.bim $plink_memory_limit --make-bed --out ${plinkset_out}  > /dev/null

grep 'loaded from' ${plinkset_out}.log  | sed  's/loaded from.*//' | head -n2
grep -e ' people pass filters and QC' ${plinkset_out}.log
printf "Output file: ${plinkset_out} \n\n"


#dups: only remove the duplicates keep the first variant(TODO keep with least missing)
plinkset_in=${plinkset_out}
plinkset_out=${plinkset_in}_dups
awk 'seen[$2]++{$2=$2"_DUPS_"seen[$2]} 1' ${plinkset_in}.bim > ${plinkset_out}_temp.bim
# dups marked
plink --fam  ${plinkset_in}.fam --bed  ${plinkset_in}.bed --bim ${plinkset_out}_temp.bim $plink_memory_limit --make-bed --out ${plinkset_out} >/dev/null

if [ "$skip_flag" = false ];
then		
	plinkset_in=${plinkset_out}
        plinkset_out=${plinkset_in}Excl
	# plink --bfile  ${plinkset_in} --exclude <( awk '{print $2}' ${plinkset_in}.bim | grep 'DUPS' | awk '{print $2}' ) \
#--make-bed --out ${plinkset_out} >/dev/null
	plink --bfile ${plinkset_in} --list-duplicate-vars suppress-first -out ${plinkset_in} $plink_memory_limit  > /dev/null
	plink --bfile ${plinkset_in} --exclude ${plinkset_in}.dupvar --make-bed --out ${plinkset_out} $plink_memory_limit >/dev/null
	printf " Dupli/Multiplicate variants excluded: $(wc -l <${plinkset_in}.dupvar)\n"
fi
        printf " Input: ${plinkset_in}\n"
        grep 'loaded from' ${plinkset_out}.log  | sed  's/loaded from.*//' | head -n2
        grep -e ' people pass filters and QC' ${plinkset_out}.log
        printf "Output file: ${plinkset_out} \n\n"

	
##### initial SNP filters #####
printf '%s\n\n' "Step 1: Remove SNPs with >5% missingness or with MAF <0.01"
output=$( printf ${output_stem}_geno05_maf01 )
plink --bfile ${plinkset_out} --geno 0.05 --maf 0.01 --make-bed --out $output $plink_memory_limit > /dev/null

grep 'loaded from' $output.log  | sed  's/loaded from.*//' | head -n2
grep -e 'missing genotype data' -e 'minor allele threshold' -e ' people pass filters and QC' $output.log
printf "Output file: $output \n\n"

##### person missingness #####
printf '%s\n' "Step 2: remove subjects w/ >1% missingness"
output_last=$output
output=${output}_mind01
plink --bfile $output_last --mind 0.01 --make-bed --out $output $plink_memory_limit > /dev/null

grep -e 'removed due to missing genotype data' -e ' people pass filters and QC' $output.log
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
printf "Removing $( wc -l ${output}_relatedness_related_ids.txt ) individuals for relatedness.\n"
output_last=$output
output=${output}_norelated
plink --bfile $output_last --remove ${output_last}_relatedness_related_ids.txt --make-bed --out $output $plink_memory_limit > /dev/null
grep ' people pass filters and QC' $output.log
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
    grep -e ' people pass filters and QC' ${output}.log
    printf "Output file: $output \n"
fi

plinkset_x_out=${output}
# confirm sex is updated
n_sex=$(awk '{print $5}' ${plinkset_x_out}.fam | sort | uniq -dc | wc -l)
sex=$(awk '{print $5}' ${plinkset_x_out}.fam | sort | uniq -dc | awk '{print $2}')

if [ $n_sex -ne 2 ];
then
  ####################### make this an error? but earlier in the pipeline
  printf "\n Warning: more or less than 2 sexes found in fam:
  ${plinkset_x_out}.fam , 1=Male, 2=Female;
  $(awk '{print $5}' ${plinkset_x_out}.fam | sort | uniq -dc )\n"
fi

################## X-chromosome Processing ###############################
printf "\n ##### X-chromosome Processing ##### \n"
##### save plinkset for x-chromosome processing #####
plinkset_x_in=${plinkset_x_out}
#plinkset_x_start=${plinkset_x_in}
plinkset_x_first=${plinkset_x_in}

# merge and split-X to exclude PARs regions
# VJ: merge and split combination takes care of any inconsistent assignment (move it before sex check)
# also if chr25(PAR) is already present than the --split-x option results in error

#### STEP 1 ####
# merge and then split chrX with position range based on build
printf "\nStep 1: Assign pseudoautosomal regions (PARS) based on build\n"
printf "\n merge and split the genotype on build ${build} coordinates preparing to exclude PARs\n"
    plinkset_x_in=${plinkset_x_out}
    plinkset_x_out=${plinkset_x_in}_mergeX
    printf " Input: ${plinkset_x_in}\n"
    printf " Output: ${plinkset_x_out}\n"
    plink --bfile ${plinkset_x_in} --merge-x no-fail --make-bed --out ${plinkset_x_out} $plink_memory_limit > /dev/null
    grep -e 'merge-x: ' ${plinkset_x_out}.log ||
    printf " --merge-x: 0 chromosome codes changed.\n no chr 25 found. Please confirm %s\n" ${plinkset_x_out}.log  >&2


    plinkset_x_in=${plinkset_x_out}
    plinkset_x_out=${plinkset_x_in}_splitX
    printf " Input: ${plinkset_x_in} \n"
    printf " Output: ${plinkset_x_out} \n"
    plink --bfile ${plinkset_x_in} --split-x ${build} no-fail --make-bed --out ${plinkset_x_out} $plink_memory_limit > /dev/null
    grep -e 'split-x: ' ${plinkset_x_out}.log ||
    printf " --split-x: 0 chromosome codes changed.\n no chr 25 found. Please confirm %s\n" ${plinkset_x_out}.log  >&2

#### STEP 2 ####
printf '%s\n\n' "Step 2: set het. haploid genotypes missing"
#set hh missing todeal with het haplotype warning (since after sex check)
# check no hh warnings on females and set hh to missing
if [  $n_sex -eq 1 ] && [ ${sex} -eq 1 ];
then
        printf "\n Only males or other sex present \n"
else
    printf "\n Confirm no hh warning for females only selection below: 'Warning: ... het. haploid genotypes present' \n"
    # will have issues with males only cohort. skip this for male only
    plink --bfile ${plinkset_x_out} --filter-females --freq $plink_memory_limit --out ${plinkset_x_out} > /dev/null
fi

# set hh missing
plinkset_x_in=${plinkset_x_out}
plinkset_x_out=${plinkset_x_out}_nohh
printf "\n set het. haploid genotypes missing in ${plinkset_x_out}\n"
plink --bfile ${plinkset_x_in} --set-hh-missing --make-bed --out ${plinkset_x_out} $plink_memory_limit > /dev/null
printf "Input file: ${plinkset_x_in} \n"
printf "Output file: ${plinkset_x_out} \n\n"

#### STEP 3 ####
printf "\nStep 3: Subset to X-chromosome \n"
plinkset_x_in=${plinkset_x_out}
plinkset_x_out=${plinkset_x_first}_X


#subset to X-chr SNPS (X:23, Y: 24, XY:25, mitochondria:26)
# TODO: * need to confirm if chr is not coded as X instead of 23
#       * may need to split-x before subsetting (confirm)
plink --bfile ${plinkset_x_in} --chr 23 --make-bed --out ${plinkset_x_out} $plink_memory_limit > /dev/null
    printf " Input: ${plinkset_x_in} \n"
    printf " Output: ${plinkset_x_out}\n"
grep -e ' people pass filters and QC.' ${plinkset_x_out}.log
n_x2=$(awk '$1==23{print}' ${plinkset_x_out}.bim | wc -l)
printf "Removed PARS SNPs: $(( ${n_x} - ${n_x2}))\n"
printf "number of X-chromosome SNPs remaining: $(awk '$1==23{print}' ${plinkset_x_out}.bim | wc -l)\n"
# exclude pseudoautosomal regions (PARs)
#n_par=$(awk '$1==25{print}' ${plinkset_x_out}.bim | wc -l)
#printf "check: ${n_par} SNPs in PARs left.(should be = 0)\n"

#### STEP 4 ####
############## This step occurs above and erroring was already commented out here. 
## confirm sex is updated
#n_sex=$(awk '{print $5}' ${plinkset_x_out}.fam | sort | uniq -dc | wc -l)
#sex=$(awk '{print $5}' ${plinkset_x_out}.fam | sort | uniq -dc | awk '{print $2}') 
#if [ $n_sex -ne 2 ];
#then
#printf "\n ERROR: sex not proper/updated in fam:
#${plinkset_x_out}.fam , 1=Male, 2=Female;
#$(awk '{print $5}' ${plinkset_x_out}.fam | sort | uniq -dc )\n"
#        #exit 1;
#fi

# update pheno with sex
printf "Updating phenotype in plink files to correspond to sex...\n"
plinkset_x_in=${plinkset_x_out}
plinkset_x_out=${plinkset_x_out}_sexpheno
awk -F ' ' '{print $1,$2,$5}' ${plinkset_x_in}.fam > ${plinkset_x_out}.txt
plink --bfile ${plinkset_x_in} --make-pheno 2 --make-bed --out ${plinkset_x_out} $plink_memory_limit > /dev/null


##### diffential missingness calculations #####
if [  $n_sex -eq 1 ];
then
        printf "\nSkipping Step 5: diffential missingness calculations #####  \n"
        printf "\n Only one sex present \n"
else
        printf "\nStep 5: ##### diffential missingness calculations #####  \n"
        # generate male-female diff miss stats
        plink --bfile ${plinkset_x_out} --test-missing --out ${plinkset_x_out}_diffmiss $plink_memory_limit > /dev/null
        # drop SNPs with P < 0.0000001
	## generate list of variants
        awk '$5 <= 0.0000001{print $2}' ${plinkset_x_out}_diffmiss.missing > ${plinkset_x_out}_diffmiss_to_drop.txt
        printf "\n male-female diff miss SNPs with P < 0.0000001(1e-7) dropped: $(wc -l < ${plinkset_x_out}_diffmiss_to_drop.txt)\n"
	## remove using plink
        plinkset_x_in=${plinkset_x_out}
        plinkset_x_out=${plinkset_x_in}_e7diffmiss
        plink --bfile ${plinkset_x_in} --exclude ${plinkset_x_in}_diffmiss_to_drop.txt --make-bed --out ${plinkset_x_out} $plink_memory_limit > /dev/null
        printf " Input: ${plinkset_x_in} \n"
        printf " Output: ${plinkset_x_out} \n"
        grep -e 'people pass filters and QC.' ${plinkset_x_out}.log

fi


##### Step 7: Remove individuals based on autosomal PCs and heterozygocity outliers (and for restrict to NHW individuals) #####
printf "\nStep 7: #### Remove individuals based on autosomal PCs and heterozygocity outliers (and for restrict to NHW individuals)  #### \n"
# skip for now as need to figure out how to pass this file that requires manual intervention
# (Done: required as input)VJ need to get the final plinkset from autosomal processing

#plinkset_autosomal_part1=${ds}_genotyped_geno05_maf01_mind01_norelated_sex_nomismatchedsex_keep_autosomes_nooutliers
#awk '{print $1,$2}' ${ds_plinkset}.fam > ${plinkset_x_out}_ind_to_keep.txt
printf "\n individuals in final autosomal dataset: $(wc -l < ${ds_plinkset}.fam) \n"

plinkset_x_in=${plinkset_x_out}
plinkset_x_out=${plinkset_x_in}_noout ##################### This should have a better name
plink --bfile ${plinkset_x_in} --keep ${ds_plinkset}.fam --make-bed --out ${plinkset_x_out} $plink_memory_limit > /dev/null
    printf " Input: ${plinkset_x_in} \n"
    printf " Output: ${plinkset_x_out} \n"
#printf "removed individuals based on autosomal PC plots and heterozygocity outliers (and for restrict to NHW individuals), plinkset:\n ${plinkset_x_out}\n"
#### Step 8: Hardy Weinberg equilibrium 1e6 (based on females) only females are diploids XX ####
if [  $n_sex -eq 1 ] && [ ${sex} -eq 1 ];
then
        printf "\nSkipping step 8: #### Hardy Weinberg equilibrium 1e6 (based on females) only females are diploids XX #### \n"
        printf "\n Only males or other sex present \n"
else
        printf "\nStep 8: #### Hardy Weinberg equilibrium 1e6 (based on females) only females are diploids XX #### \n"
        plinkset_x_in=${plinkset_x_out}
        plinkset_x_out_fem=${plinkset_x_in}_1e6femHWE
	plinkset_x_out=${plinkset_x_in}_femHWE6
        # get SNPS out of HWE in females
        plink --bfile ${plinkset_x_in} --filter-females --hwe 0.000001 --make-just-bim --out ${plinkset_x_out_fem} $plink_memory_limit > /dev/null
        awk '{print $2}' ${plinkset_x_out_fem}.bim > ${plinkset_x_out_fem}_SNPs.txt
        plink --bfile ${plinkset_x_in} --extract ${plinkset_x_out_fem}_SNPs.txt --make-bed --out ${plinkset_x_out} $plink_memory_limit > /dev/null
        grep -e 'extract:' ${plinkset_x_out}.log
        grep -e 'variants loaded from .bim file' ${plinkset_x_out}.log
        grep -e 'people pass filters and QC.' ${plinkset_x_out}.log
        printf " Output: ${plinkset_x_out} \n"
fi


############################# insert preparation for imputation

printf "\nStep 2: Removing palindromic variants \n"

#get palindromic variants
awk '{ if( ($5 == "A" && $6 == "T") || ($5 == "T" && $6 == "A") || ($5 == "C"  && $6 == "G") || ($5 == "G" && $6 == "C")) print }'  ${plinkset_x_out}.bim > ${output_path}/palindromic_snps.txt
#remove them
output_last=$plinkset_x_out
output=${plinkset_x_out}_nopal
plink --bfile $output_last --exclude ${output_path}/palindromic_snps.txt --output-chr M --make-bed --out $output > /dev/null
printf "Removed $( wc -l < ${output_path}/palindromic_snps.txt ) palindromic variants.\n"
grep ' people pass filters and QC' ${output}.log
printf "Output file: $output \n"

#check the build, decide whether to do the liftOver and, if so, which chain file to use
if [ "$build" = "b37" ];
then
    printf "\nStep 3: Lifting over the genotype file from build 37 to build 38 to match the TOPMed reference panel\n"
    chain_file="hg19ToHg38.over.chain.gz"
elif [ "$build" = "b36" ];
then 
    printf "\nStep 2: Lifting over the genotype file from build 36 to build 38 to match the TOPMed reference panel\n"
    chain_file="hg18ToHg38.over.chain.gz"
elif [ "$build" = "b38" ];
then
    printf "\nStep 3: Input data was specified as already on build 38, so no lift-over is necessary. Proceeding to the next step.\n"
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
    printf "Removing $( wc -l < ${output}_lift_issue_snps_todrop.txt ) variants which fail liftOver ($( grep "random" ${output}_lifted.txt | wc -l ) 'random', $( grep "alt" ${output}_lifted.txt | wc -l ) 'alt', $( grep "fix" ${output}_lifted.txt | wc -l ) 'fix', and $( awk '{ print $4 }' ${output}_unlifted.txt | sort -u | wc -l ) not present in b38).\n"

    #drop 
    plink --bfile $output --exclude ${output}_lift_issue_snps_todrop.txt --output-chr M --make-bed --out ${output}_noprobSNPs > /dev/null
    #update the chr and BP positions for remaining
    plink --bfile ${output}_noprobSNPs --update-chr ${output}_lifted.txt 1 4 --output-chr M --make-bed --out ${output}_noprobSNPs_chr > /dev/null
    plink --bfile ${output}_noprobSNPs_chr --update-map ${output}_lifted.txt 3 4 --make-bed --out ${output}_noprobSNPs_chr_bplifted > /dev/null
    output=${output}_noprobSNPs_chr_bplifted
fi

#### remove same-position variants ####

printf "\nStep 4: Removing multi-allelic and duplicated variants.\n"
awk '{ print $2" "$1"_"$4 }' ${output}.bim | sort -T ${output_path}/ -k2 | uniq -f1 -D | awk '{ print $1 }' > ${output}_samepos_vars.txt

if [ "$( wc -l < ${output}_samepos_vars.txt )" -gt 0 ];
then
    printf "Removing $( wc -l < ${output}_samepos_vars.txt ) same-position variants.\n"
    plink --bfile ${output} --exclude ${output}_samepos_vars.txt --make-bed --out ${output}_nosamepos > /dev/null
    output=${output}_nosamepos
    grep ' people pass filters and QC' $output.log
    printf "Output file: $output\n"
else 
    printf "No multi-allelic or duplicated variants to remove.\n"
fi

#### Imputation prep ####

printf "\nStep 5: Comparing with the reference panel and preparing files for imputation for each chromosome.\n"

#make set with short name and freq file
plink --bfile ${output} --freq --make-bed --out ${output_stem} > /dev/null

#run the imputation checking script
perl HRC-1000G-check-bim.pl -b ${output_stem}.bim  -f ${output_stem}.frq -r ${ref_file} -h -n > /dev/null

#if noexclude, then skip the exclusion step in the plink file
if [ "$noexclude" = true ];
then
    sed -i -e '1s/.*/#&/' -e "s|${output_path}/TEMP1|${output_stem}|g" ${output_path}/Run-plink.sh
fi

#run created script with all the plink commands
#muting all output since this will throw quite a few warnings if there are any funky alleles
sed -i "s|rm TEMP|rm ${output_path}/TEMP|" ${output_path}/Run-plink.sh
sh ${output_path}/Run-plink.sh > /dev/null 2>&1


##### Subset to Overlapping SNPs between sexes #####
printf "\nStep: Subset to Overlapping SNPs between sexes \n\n"
output=${output_stem}-updated-chr23

# skip overlapping SNP step for single sex dataset
n_sex=$(awk '{print $5}' ${output}.fam | sort | uniq -dc | wc -l)
sex=$(awk '{print $5}' ${output}.fam | sort | uniq -dc | awk '{print $2}') 
if [ $n_sex -ne 2 ];
then
    printf "Warning: single sex dataset skipping subset to overlapping SNPs between sexes \n"
    printf "N sexes: $n_sex, sex code: $sex \n"
else
    plink --bfile ${output} --filter-females --geno 0.01 --write-snplist --out ${output}_females $plink_memory_limit > /dev/null
    plink --bfile ${output} --filter-males --geno 0.01 --write-snplist --out ${output}_males $plink_memory_limit > /dev/null
    # sort those files and get overlapping variant IDs        
    cat ${output}_males.snplist | sort > ${output}_males_sorted.snplist
    cat ${output}_females.snplist | sort > ${output}_females_sorted.snplist
    comm -12 ${output}_males.snplist ${output}_females.snplist  > ${output}_overlapping_SNPs.txt
 
    # Keep the allele order defined in the .bim file,instead of forcing A2 to be the major allele.--real-ref-alleles also removes 'PR' from the INFO values emitted by --recode vcf{-fid/-iid}.
    ## there was a filter (--chr 23) here but it was removed because there already should just be chr 23 variants present here
    output_last=$output
    output=${output}_overlap
    plink --bfile ${output_last} --extract ${output_last}_overlapping_SNPs.txt --real-ref-alleles --make-bed --out ${output} $plink_memory_limit > /dev/null
    printf "Output: ${output} \n"
    printf " female SNPs: $(wc -l < ${output_last}_female_snps.txt) \n"
    printf " male SNPs: $(wc -l < ${output_last}_male_snps.txt) \n"
    printf " Overlapping SNPs in both sexes: $(wc -l < ${output_last}_overlapping_SNPs.txt) \n"
    grep -e ' people pass filters and QC' ${output}.log

fi

printf "\n update variant IDs to TOPMED imputation server format chrX:POS:REF:ALT \n"
# the genotyped variant ids in TOPMED imputation server format chrX:POS:REF:ALT
output_last=$output
output=${output}_TOPMED_varID

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
mkdir ${output_path}/tmpX/
bcftools annotate --rename-chrs update_chr_names_b38.txt ${output}.vcf | bcftools sort --temp-dir  ${output_path}/tmpX/ -O z -o ${output}_chrX.vcf.gz > /dev/null
printf "Sorting complete.\n\n"

#print out number of variants actually excluded or which would have been excluded
if [ "$noexclude" = true ];
then 
    printf "Would have removed $( cat ${output_path}/Exclude-* | wc -l ) variants for mismatch with the reference panel, being palindromic with MAF > 0.4, or being absent from the reference panel leaving $( cat ${output}-updated-chr23_overlap.bim | wc -l ) for imputation, but the no-exclude option was specified.\n"
else
    printf "Removed $( cat ${output_path}/Exclude-* | wc -l ) variants for mismatch with the reference panel, being palindromic with MAF > 0.4, or being absent from the reference panel leaving $( cat ${output}-updated-chr23_overlap.bim | wc -l ) for imputation.\n"
fi

# commenting this out to prevent removal of files that ought to be kept
#remove the intermediate .vcf and .bed files
##make sure that it only removes the files from running this script, not the initial files or files for other datsets
#if [ "$skip_cleanup" = false ];
#then
#    rm ${output_path}/*.vcf
#    files_to_remove=$( find ${output_path}/${input_stem}*.bed | grep -v "${output}.bed" | grep -v "${output_stem}-updated.bed" )
#    rm $files_to_remove
#else
#    printf "Skipping cleanup. Please manually remove unnecessary files.\n"
#fi

printf "\nConversion complete! Upload the file (${output_stem}_chrX.vcf.gz) to the imputation server.\n"
