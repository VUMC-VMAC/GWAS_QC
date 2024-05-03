#!/bin/bash
# Author: Vaibhav Janve 2020-04-07
# Author 2: Emily Mahoney 2020-04-09
# Aug 2021, Vaibhav Janve: script update to process X-chr


# VJ: singularity image and directories bound for reproducibility
[ ! -z "$SINGULARITY_CONTAINER" ] && printf "\nSingularity image: $SINGULARITY_CONTAINER"
[ ! -z "$SINGULARITY_BIND" ] && printf "\nMapped directories:\n $SINGULARITY_BIND\n" | sed 's/,/\n /g'

#fail on error
set -e

#define usage
display_usage() {
    printf "GWAS QC Part 1 (including X-chromosome)

Completes the first stage (including X-chromosome) in standard GWAS QC, including initial variant and person filters, relatedness and sex checks, restriction to autosomes, and PC calculation.
the phenotype column in .fam is updated with sex in process.

Usage:
gwas_qc_part1x.sh -o [output_stem] -i [input_fileset] -r [race_file] -f [sex_file] -G [stem_1000G] -b [b37] -p [ds_plinkset] -s

output_stem = the beginning part of all QC'ed files, including the full file path to the directory where the files are to be saved

input_fileset = the full path and file stem for the raw plink set '*[bed,bim,fam]'

sex_file = a file with FID and IID (corresponding to the fam file), 1 column indicating sex for the sex check (1 for males, 2 for females, 0 if unknown), with NO header.

race_file (optional) = a file with FID and IID (corresponding to the fam file) and 1 column indicating both race and ethnicity for PC plots, with NO header. Non-hispanic whites need to be indicated with 'White' or 'EUR.' No other values in the race column must be fixed; however, the race column must not include spaces. This is only needed if you want to color PC plots based on race (which at this stage will only be self-report). Ancestral categories will be calculated post-imputation using SNPWeights. 


stem_1000G (optional) = the full path and stem to the 1000G genotype files in plink format. There must also be a file with FID, IID, race with the same stem and _race.txt as the suffix (i
e for a plink file set like this: all_1000G.bed, all_1000G.bim, all_1000G.fam the race file would be like this all_1000G_race.txt)
-b = optional argument indicating the build of the input dataset; default assumption is b37 but can supply b36 or b38. This will impact the x-chromosome split step.
-p = autosomal plinkset stem with IID and FID matching the input plinkset and PC outliers removed
-s optional, to skip dup removal step after initial id format
-h will show this usage
"
        }

skip_flag='false'
while getopts 'o:i:r:f:G:b:p:sh' flag; do
  case "${flag}" in
    o) output_stem="${OPTARG}" ;;
    i) input_fileset="${OPTARG}" ;;
    r) race_file="${OPTARG}" ;;
    f) sex_file="${OPTARG}" ;;
    G) stem_1000G="${OPTARG}" ;;
    p) ds_plinkset="${OPTARG}" ;;
    b) build="${OPTARG}" ;;
    s) skip_flag='true' ;;
    h) display_usage ; exit ;;
    \?|*) display_usage
       exit 1;;
  esac
done


# X-chromosome # of SNP's threshold
nx_threshold=300
skip_sexMAF="yes"
skip_heterozygocity_outlier_check="yes"
x_process_flag=1


#check to make sure necessary arguments are present
if [ -z "$output_stem" ] || [ -z "$input_fileset" ] || [ -z "$sex_file" ] ;
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
Step7: Remove individuals based on autosomal PCs and heterozygocity outliers (and for restrict to NHW individuals)
Step8: [if females are present] Hardy Weinberg equilibrium(HWE) 1e6 (based on females)
        drop SNPS in both males and females that are out of HWE."

plinkset_in=${input_fileset}

# get X-chromosome SNP count
n_x=$(awk '$1==23{print}' ${plinkset_in}.bim | wc -l)

#print out inputs
printf "GWAS QC Part 1 (including X-chromosome) Script

Input data : $input_fileset
Output stem : $output_stem
File with sex information : $sex_file
Number of X-chromosome SNPs: ${n_x}
"
if [[ ${n_x} -lt 300 ]]; 
then
  printf "%s\n" "CAUTION: too few X-chromosomes, script might error at sex check";
fi


# print message for input build

if [ -z ${build} ]
then
   printf " Input data build not specified assuming 'b37' \n"
   build="b37"
else
   printf " Input data build: ${build} \n"
fi

#print message if 1000G dataset is not specified
if [ -z "$stem_1000G" ];
then
    printf 'No location was specified for the 1000G data, so no PCs will be calculated including them!\n\n'
else
    echo $stem_1000G
    printf "1000G data for PC calculation : ${stem_1000G}\n\n"
fi

#print message if the race file is not specified
if [ -z "$race_file" ];
then
    printf "No file was supplied with self-report race information, so PCs will not be colored based on these values.\n\n"
else
    printf "Self-report race values will be drawn from ${race_file} in order to color PCs.\n\n"
fi


#check to make sure this is being run in the scripts folder (checking if necessary script is present)
# VJ: these checks may not be necessary as we have moved scripts to container
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


if [ ${n_x} -le ${nx_threshold} ];
then
        printf "SKIP X-chromosome processing: too few X-chromosome SNPs (${n_x}) to process (threshold ${nx_threshold})\n";
        xprocess_flag=0;
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
	plink --bed ${plinkset_in}.bed --fam ${plinkset_in}.fam --bim ${plinkset_out}_temp.bim --make-bed --out ${plinkset_out}  > /dev/null
	
	printf " Input: ${plinkset_in}\n"
	grep 'loaded from' ${plinkset_out}.log  | sed  's/loaded from.*//' | head -n2
	printf " Output: ${plinkset_out} \n\n"


        # long indel: smart pca limit 39 characters so keeping limit to 25
        plinkset_in=${plinkset_out}
        plinkset_out=${plinkset_in}_indel
        awk '{if( length($5)>25 || length($6)>25) { print $1,$1"_"$4"_I_D",$3,$4,$5,$6} else {print $0} }' ${plinkset_in}.bim > ${plinkset_out}_temp.bim
        plink --bed ${plinkset_in}.bed --fam ${plinkset_in}.fam --bim ${plinkset_out}_temp.bim --make-bed --out ${plinkset_out}  > /dev/null

        grep 'loaded from' ${plinkset_out}.log  | sed  's/loaded from.*//' | head -n2
        grep -e ' people pass filters and QC' ${plinkset_out}.log
        printf "Output file: ${plinkset_out} \n\n"


        #dups: only remove the duplicates keep the first variant(TODO keep with least missing)
        plinkset_in=${plinkset_out}
        plinkset_out=${plinkset_in}_dups
                awk 'seen[$2]++{$2=$2"_DUPS_"seen[$2]} 1' ${plinkset_in}.bim > ${plinkset_out}_temp.bim
        # dups marked
        plink --fam  ${plinkset_in}.fam --bed  ${plinkset_in}.bed --bim ${plinkset_out}_temp.bim --make-bed --out ${plinkset_out} >/dev/null
if [ "$skip_flag" = false ];
then		
	plinkset_in=${plinkset_out}
        plinkset_out=${plinkset_in}Excl
	# plink --bfile  ${plinkset_in} --exclude <( awk '{print $2}' ${plinkset_in}.bim | grep 'DUPS' | awk '{print $2}' ) \
#--make-bed --out ${plinkset_out} >/dev/null
	plink --bfile ${plinkset_in} --list-duplicate-vars suppress-first -out ${plinkset_in}  > /dev/null
	plink --bfile ${plinkset_in} --exclude ${plinkset_in}.dupvar --make-bed --out ${plinkset_out} >/dev/null
	printf " Dupli/Multiplicate variants excluded: $(wc -l <${plinkset_in}.dupvar)\n"
fi
        printf " Input: ${plinkset_in}\n"
        grep 'loaded from' ${plinkset_out}.log  | sed  's/loaded from.*//' | head -n2
        grep -e ' people pass filters and QC' ${plinkset_out}.log
        printf "Output file: ${plinkset_out} \n\n"

	
##### initial SNP filters #####
printf '%s\n\n' "Step 1: Remove SNPs with >5% missingness or with MAF <0.01"
output=$( printf ${output_stem}_geno05_maf01 )
plink --bfile ${plinkset_out} --geno 0.05 --maf 0.01 --make-bed --out $output > /dev/null

grep 'loaded from' $output.log  | sed  's/loaded from.*//' | head -n2
grep -e 'missing genotype data' -e 'minor allele threshold' -e ' people pass filters and QC' $output.log
printf "Output file: $output \n\n"

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
if [ "$( awk '{ print $3 }' $sex_file | sort -u | tr -d '[:space:]' )" != 0 ];
then
    output_last=$output
    output=${output}_sex
    plink --bfile $output_last --update-sex $sex_file --make-bed --out $output > /dev/null
    printf "Updated sex from the provided file: $sex_file \n"
else
    printf "No non-missing sex information in the provided file (${sex_file}). Using sex from fam file.\n"
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

plinkset_x_out=${output}
# confirm sex is updated
n_sex=$(awk '{print $5}' ${plinkset_x_out}.fam | sort | uniq -dc | wc -l)
sex=($(awk '{print $5}' ${plinkset_x_out}.fam | sort | uniq -dc | awk '{print $2}') )

if [ $n_sex -ne 2 ];
then
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
    plink --bfile ${plinkset_x_in} --merge-x no-fail --make-bed --out ${plinkset_x_out} --memory 15000 > /dev/null
    grep -e 'merge-x: ' ${plinkset_x_out}.log ||
    printf " --merge-x: 0 chromosome codes changed.\n no chr 25 found. Please confirm %s\n" ${plinkset_x_out}.log  >&2


    plinkset_x_in=${plinkset_x_out}
    plinkset_x_out=${plinkset_x_in}_splitX
    printf " Input: ${plinkset_x_in} \n"
    printf " Output: ${plinkset_x_out} \n"
    plink --bfile ${plinkset_x_in} --split-x ${build} no-fail --make-bed --out ${plinkset_x_out} --memory 15000 > /dev/null
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
    plink --bfile ${plinkset_x_out} --filter-females --freq --memory 15000 --out ${plinkset_x_out} > /dev/null
fi

# set hh missing
plinkset_x_in=${plinkset_x_out}
plinkset_x_out=${plinkset_x_out}_nohh
printf "\n set het. haploid genotypes missing in ${plinkset_x_out}\n"
plink --bfile ${plinkset_x_in} --set-hh-missing --make-bed --out ${plinkset_x_out} --memory 15000 > /dev/null
printf "Input file: ${plinkset_x_in} \n"
printf "Output file: ${plinkset_x_out} \n\n"

#### STEP 3 ####
printf "\nStep 3: Subset to X-chromosome \n"
plinkset_x_in=${plinkset_x_out}
plinkset_x_out=${plinkset_x_first}_X


#subset to X-chr SNPS (X:23, Y: 24, XY:25, mitochondria:26)
# TODO: * need to confirm if chr is not coded as X instead of 23
#       * may need to split-x before subsetting (confirm)
plink --bfile ${plinkset_x_in} --chr 23 --make-bed --out ${plinkset_x_out} --memory 15000 > /dev/null
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
# confirm sex is updated
n_sex=$(awk '{print $5}' ${plinkset_x_out}.fam | sort | uniq -dc | wc -l)
sex=($(awk '{print $5}' ${plinkset_x_out}.fam | sort | uniq -dc | awk '{print $2}') )
if [ $n_sex -ne 2 ];
then
printf "\n ERROR: sex not proper/updated in fam:
${plinkset_x_out}.fam , 1=Male, 2=Female;
$(awk '{print $5}' ${plinkset_x_out}.fam | sort | uniq -dc )\n"
        #exit 1;
fi

# update pheno with sex
# create new fam file to use
awk -F ' ' '{print $1,$2,$3,$4,$5,$5}' ${plinkset_x_out}.fam > ${plinkset_x_out}_pheno.fam
printf "\n modified fam: updated phenotype in ${plinkset_x_out}_pheno.fam with sex \n"

#### male-female differential MAF : delta > 0.02 ####
if [ $skip_sexMAF = "yes" ]
then
            printf "\nSkipping Step 4: #### male-female differential MAF : delta > 0.02 #### \n"
else
    if [  $n_sex -eq 1 | ];
    then
            printf "\nSkipping Step 4: #### male-female differential MAF : delta > 0.02 #### \n"
            printf "\n Only one sex present \n"
    else
            printf "\nStep 4: #### male-female differential MAF : delta > 0.02 #### \n"
            #VJ: since we have the sex as pheno type we can generate case-control frequency
                    #we get: ${plinkset_x_out}.frq
            printf "Generate case-control phenotype-stratified allele frequency report with sex as phenotype: ${plinkset_x_out}.frq.cc \n"
            plink --bed ${plinkset_x_out}.bed --bim ${plinkset_x_out}.bim --fam ${plinkset_x_out}_pheno.fam --freq case-control --out ${plinkset_x_out} --memory 15000 > /dev/null
                    #we get: ${plinkset_x_out}.frq.cc
            # generating genotype count report
            printf " Generating genotype count report: ${plinkset_x_out}.frqx \n"
            plink --bfile ${plinkset_x_out} --freqx --out ${plinkset_x_out} --memory 15000 > /dev/null
            printf " Generating basic allele frequency report: ${plinkset_x_out}.frq \n"
            plink --bfile ${plinkset_x_out} --freq --out ${plinkset_x_out} --memory 15000 > /dev/null
            awk 'function abs(x){return ((x < 0.0) ? -x : x)} {print $0, abs($5-$6)}' OFS='\t' ${plinkset_x_out}.frq.cc | \
            awk '$9>0.02 {print $2}' > ${plinkset_x_out}_MAFdiff02_to_exclude.txt

            printf "total x-chr SNPs: $( awk '$1==23{print}' ${plinkset_x_out}.frq.cc | wc -l) \n"
            printf "x-chr SNPs to exclude for male-female MAF difference of greater than 0.02: \
            $( wc -l < ${plinkset_x_out}_MAFdiff02_to_exclude.txt)\n"

            plinkset_x_in=${plinkset_x_out}
            plinkset_x_out=${plinkset_x_in}_MAFdiff02

            plink --bfile ${plinkset_x_in} --exclude ${plinkset_x_in}_MAFdiff02_to_exclude.txt --make-bed --out ${plinkset_x_out} --memory 15000 > /dev/null
                printf " Input: ${plinkset_x_in} \n"
                printf " Output: ${plinkset_x_out} \n"
            grep -e "people pass filters and QC." ${plinkset_x_out}.log
    fi
fi
##### diffential missingness calculations #####
if [  $n_sex -eq 1 ];
then
        printf "\nSkipping Step 5: ##### diffential missingness calculations #####  \n"
        printf "\n Only one sex present \n"
else
        printf "\nStep 5: ##### diffential missingness calculations #####  \n"
        awk -F ' ' '{print $1,$2,$3,$4,$5,$5}' ${plinkset_x_out}.fam > ${plinkset_x_out}_pheno.fam
        printf "\n modified fam: updated phenotype in ${plinkset_x_out}_pheno.fam with sex \n"
        # generate male-female diff miss stats
        plink --bed ${plinkset_x_out}.bed --bim ${plinkset_x_out}.bim --fam ${plinkset_x_out}_pheno.fam --test-missing --out ${plinkset_x_out}_diffmiss --memory 15000 > /dev/null
        grep -B4 "... done" ${plinkset_x_out}_diffmiss.log
        # drop SNPs with P < 0.0000001
        awk '$5 <= 0.0000001{print $2}' ${plinkset_x_out}_diffmiss.missing > ${plinkset_x_out}_diffmiss_to_drop.txt
        printf "\n male-female diff miss SNPs with P < 0.0000001(1e-7) dropped: $(wc -l < ${plinkset_x_out}_diffmiss_to_drop.txt)\n"

        plinkset_x_in=${plinkset_x_out}
        plinkset_x_out=${plinkset_x_in}_e7diffmiss
        plink --bfile ${plinkset_x_in} --exclude ${plinkset_x_in}_diffmiss_to_drop.txt --make-bed --out ${plinkset_x_out} --memory 15000 > /dev/null
            printf " Input: ${plinkset_x_in} \n"
            printf " Output: ${plinkset_x_out} \n"
        grep -e 'people pass filters and QC.' ${plinkset_x_out}.log

fi
# skipping as autosomal QC is taking care of heterozygocity outliers
# printf "\nStep6: #### skipping heterozygocity check (males only?,skipping as autosomal QC is taking care of heterozygocity outliers and\n --het requires autosomal ####\n"


##### heterozygocity check (males only?, using all for now) #####
if [ ${skip_heterozygocity_outlier_check}="yes" ]
then 
    printf "\nSkipping Step6: Heterozygocity outlier check \n"
else
  if [  $n_sex -eq 1 ] && [ ${sex} -eq 2 ];
  then
        printf "\nSkipping Step6: #### Heterozygocity check (males only) #### \n"
        printf "\n Only females or other sex present \n"
  else


        #printf "\nStep6: #### Heterozygocity check (males only?, using all for now) ####\n"
        printf "\nStep6: #### Heterozygocity check (males only) #### \n"
        printf " PS: plink requires autosomal chromosome for --het calculations,  using dog as species to treat X-chromosome as autosomal for this step.
         and circumvent 'Error: --het requires at least one polymorphic autosomal marker.' \n"
        ## filter males and prune set (removed --filter-males flag)
        #plink --bfile ${plinkset_x_out} --indep-pairwise 200 100 0.2 --allow-no-sex --out ${plinkset_x_out}_prune --memory 15000 > /dev/null
        plink --bfile ${plinkset_x_out} --filter-males --indep-pairwise 200 100 0.2 --allow-no-sex --out ${plinkset_x_out}_prune --memory 15000 > /dev/null
        plink --bfile ${plinkset_x_out} --output-missing-phenotype 1 --extract ${plinkset_x_out}_prune.prune.in --make-bed --out ${plinkset_x_out}_pruned > /dev/null
        rm ${plinkset_x_out}_prune.*
        printf "$( wc -l < ${plinkset_x_out}_pruned.bim ) variants(not males only) out of $( wc -l < ${plinkset_x_out}.bim ) left after pruning.\n"


        ###### heterozygosity check #####
        # using "dog" as species to treat X as autosomal and overcome "Error: --het requires at least one polymorphic autosomal marker."
        plink --bfile ${plinkset_x_out}_pruned --het --dog --out ${plinkset_x_out}_pruned_hetcheck > /dev/null
        ##plot and check for outliers
        Rscript plot_het_check_outliers.R ${plinkset_x_out}_pruned_hetcheck
        printf "Heterozygocity outliers(X-chr, males only) review figure: ${plinkset_x_out}_pruned_hetcheck.png \n"


        #if there are outliers >6 sd from the F stat mean, remove them
        if [ -f "${plinkset_x_out}_pruned_hetcheck_outliers.txt" ];
        then
            plinkset_x_in=${plinkset_x_out}
            plinkset_x_out=${plinkset_x_out}_nohetout
            plink --bfile ${plinkset_x_in} --remove ${plinkset_x_in}_pruned_hetcheck_outliers.txt --make-bed --out ${plinkset_x_out} --memory 15000 > /dev/null
            printf " Input: ${plinkset_x_in} \n"
            grep -e ' people pass filters and QC' ${plinkset_x_out}.log
            printf "Output: ${plinkset_x_out} \n"
        fi
  fi
fi
##### Step 7: Remove individuals based on autosomal PCs and heterozygocity outliers (and for restrict to NHW individuals) #####
printf "\nStep 7: #### Remove individuals based on autosomal PCs and heterozygocity outliers (and for restrict to NHW individuals)  #### \n"
# skip for now as need to figure out how to pass this file that requires manual intervention
# (Done: required as input)VJ need to get the final plinkset from autosomal processing


#plinkset_autosomal_part1=${ds}_genotyped_geno05_maf01_mind01_norelated_sex_nomismatchedsex_keep_autosomes_nooutliers
awk '{print $1,$2}' ${ds_plinkset}.fam > ${plinkset_x_out}_ind_to_keep.txt
printf "\n individuals in final autosomal dataset(individual list: ${plinkset_x_out}_ind_to_keep.txt): $(wc -l < ${plinkset_x_out}_ind_to_keep.txt) \n"

plinkset_x_in=${plinkset_x_out}
plinkset_x_out=${plinkset_x_in}_noout
plink --bfile ${plinkset_x_in} --keep ${plinkset_x_in}_ind_to_keep.txt --make-bed --out ${plinkset_x_out} --memory 15000 > /dev/null
    printf " Input: ${plinkset_x_in} \n"
    printf " Output: ${plinkset_x_out} \n"
printf "removed individuals based on autosomal PC plots and heterozygocity outliers (and for restrict to NHW individuals), plinkset:\n ${plinkset_x_out}\n"
#### Step 8: Hardy Weinberg equilibrium 1e6 (based on females) only females are diploids XX ####
if [  $n_sex -eq 1 ] && [ ${sex} -eq 1 ];
then
        printf "\nSkipping step 8: #### Hardy Weinberg equilibrium 1e6 (based on females) only females are diploids XX #### \n"
        printf "\n Only males or other sex present \n"
else
        printf "\nStep 8: #### Hardy Weinberg equilibrium 1e6 (based on females) only females are diploids XX #### \n"
        plinkset_x_in=${plinkset_x_out}
        plinkset_x_out=${plinkset_x_in}_1e6femHWE
        # get SNPS out of HWE in females
        plink --bfile ${plinkset_x_in} --filter-females --hwe 0.000001 --make-bed --out ${plinkset_x_out} --memory 15000 > /dev/null
        awk '{print $2}' ${plinkset_x_out}.bim > ${plinkset_x_out}_SNPs.txt
                printf " Input: ${plinkset_x_in} \n"
        grep -e 'variants loaded from .bim file' ${plinkset_x_out}.log
        grep -e ' removed due to Hardy-Weinberg exact test.' ${plinkset_x_out}.log
        #grep -e 'people pass filters and QC.' ${plinkset_x_out}.log
        plink --bfile ${plinkset_x_in} --extract ${plinkset_x_out}_SNPs.txt --make-bed --out ${plinkset_x_out} --memory 15000 > /dev/null
        grep -e 'extract:' ${plinkset_x_out}.log
        grep -e 'variants loaded from .bim file' ${plinkset_x_out}.log
        grep -e 'people pass filters and QC.' ${plinkset_x_out}.log
        printf " Output: ${plinkset_x_out} \n"
fi
plinkset_x_in=${plinkset_x_out}
plinkset_x_out=${output_stem}
plink --bfile ${plinkset_x_in} --make-bed --out ${plinkset_x_out} --memory 15000 > /dev/null
printf " Input: ${plinkset_x_in} \n"
grep -e 'people pass filters and QC.' ${plinkset_x_out}.log
    printf " Output: ${plinkset_x_out} \n"
printf "\n Final X-chromosome plinkset:\n ${plinkset_x_out} \n Part 1 ... Done! \n"
