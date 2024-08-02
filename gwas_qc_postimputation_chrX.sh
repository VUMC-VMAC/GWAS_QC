#!/bin/bash
# VJ: 20200709 modified post imputation script for single chromosome run
# ERM: 20240705 modified for SNPWeights update

#fail on error
set -e

#define usage
display_usage() {
    printf "GWAS QC Post-imputation Script, X Chromosome

This script will unzip the imputation results for single chromosome specified (assuming password is saved in pass.txt in the same folder as the imputation results files and will perform standard post-imputation QC for our common variant pipeline. This includes filtering for R2, removing multi-allelic variants and filtering out variants for low MAF or HWE disequilibrium. Finally, PCs will be calculated on the final file-set.

Usage:
gwas_qc_postimputation_chrX.sh -o [output_stem] -i [imputation_results_folder] -f [sex_file] -s [snp_names_file] -g [preimputation_geno_X] -c [subset] -l [lab] -p [autosomal_stem] -z -x -d -m [plink_memory_limit]

output_stem = the beginning part of all QC'ed files including the full path to the folder in which they should be created

imputation_results_folder = the folder to which the imputation results have been downloaded which also contains a file called pass.txt with the password to unzip the imputation results files

sex_file = a file with FID and IID (corresponding to the fam file), 1 column indicating sex for the sex check (1 for males, 2 for females, 0 if unknown), with NO header.

snp_names_file = the file stem for converting the SNP names from imputation results to rs numbers. There should be one for each chromosome and each must have 2 columns: imputation result SNP ids and rs numbers. Can have header but it will be ignored.

preimputation_geno_X = the full path and stem to the cleaned final pre-imputation files for X chromosome to be merged back into the final files

subset = the full path and name of the file with IDs for the current subset, typically output from the assign_ancestral_categories.R script

lab = label for the current subset, typically one of EUR, AA, LatHisp, CaribHisp, etc. 

autosomal_stem = the full path and stem to the fam file of the final autosomal plinkset (with PC outliers removed)

plink_memory_limit (optional) = argument indicating the memory limit for plink to use rather than the default of half the RAM. This is useful for running this step of QC locally.

-z indicates that the imputation results will need to be unzipped. All *.zip files in imputation_results_folder will be unzipped with password provided in pass.txt
-x indicates to skip the first filtering of the individual chr files. This comes in handy if there were an issue with the next step(s) because this first step is the longest. If you supply this flag, please give the full stem of the file that this run should start with to the output_stem command.
-d will skip the clean-up at the end of the script which removes intermediate *.bed files.
-h will display this message
\n\n"
        }

#parse options
do_unzip='false'
skip_first_filters='false'
skip_cleanup='false'
while getopts 'o:i:f:s:g:p:l:m:c:zxdh' flag; do
  case "${flag}" in
    o) output_stem="${OPTARG}" ;;
    i) imputation_results_folder="${OPTARG}" ;;
    f) sex_file="${OPTARG}" ;;
    s) snp_names_file="${OPTARG}" ;;\
    z) do_unzip='true' ;;
    g) preimputation_geno_X="${OPTARG}" ;;
    c) subset="${OPTARG}" ;;
    p) autosomal_stem="${OPTARG}" ;;
    l) lab="${OPTARG}" ;;
    m) plink_memory_limit="${OPTARG}";;
    x) skip_first_filters='true' ;;
    d) skip_cleanup='true' ;;
    h) display_usage ; exit ;;
    \?|*) display_usage
       exit 1;;
  esac
done

printf "GWAS QC Post-imputation, X Chromosome\n\n"

# Print out singularity information for reproducability
[ ! -z "$SINGULARITY_CONTAINER" ] && printf "\nSingularity image: $SINGULARITY_CONTAINER"
[ ! -z "$SINGULARITY_COMMAND" ] && printf "\nSingularity shell: $SINGULARITY_COMMAND"
[ ! -z "$SINGULARITY_BIND" ] && printf "\nMapped directories:\n $SINGULARITY_BIND\n" | sed 's/,/\n /g'

#check to make sure necessary arguments are present                                                                                                    
if [ -z "$output_stem" ] || [ -z "$imputation_results_folder" ] || [ -z "$sex_file" ] || [ -z "${snp_names_file}" ] ;
then
    printf "Error: Necessary arguments not present! \n\n"
    display_usage
    exit 1
fi

#print out required inputs
printf "\nOutput file path and stem for cleaned imputed files : $output_stem
Imputation results folder : ${imputation_results_folder}
Sex information file : ${sex_file}
Stem for files for SNP name conversion : ${snp_names_file}_chrX
File for subsetting to the current set: ${subset}
Label for the current set: ${lab}
Autosomal plink files for removing heterozygosity/PC outliers: ${autosomal_stem}\n"

# indicate whether results will be unzipped
if [ "$do_unzip" = 'false' ];
then 
    printf "The -z was not specified so imputation results must already be unzipped.\n"
fi

#validate the genotyped files
if [ ! -f "${preimputation_geno_X}.fam" ];
then
    printf "Cannot see the preimputation genotype files ($preimputation_geno_X)! Please check the argument supplied to -g and try again! \n"
    exit 1
fi

## process arguments for the current subsets
if [ "$( echo $subset | grep ',' )" ];
then
    printf "Multiple subsets supplied. Validating these inputs...\n"

    # make arrays from the comma separated arguments
    IFS=',' read -r -a subset <<< "$subset"
    IFS=',' read -r -a lab <<< "$lab"
    IFS=',' read -r -a autosomal_stem <<< "$autosomal_stem"
    
    #test array lengths to make sure they are the same
    if [ "${#subset[@]}" = "${#lab[@]}" ] && [ "${#lab[@]}" = "${#autosomal_stem[@]}" ]; 
    then
	for index in "${!subset[@]}"
	do
	    #validate the file for the current subset
	    if [ ! -f "${subset[index]}" ];
	    then
		printf "Cannot see the file to subset to the subset (${subset[index]})! Please check the argument supplied to -c and try again! \n"
		exit 1
	    fi

	    #validate the final plinkset files(only .fam needed)
	    if [ ! -f "${autosomal_stem[index]}.fam" ];
	    then
		printf "Cannot see the final autosomal genotype files (${autosomal_stem[index]}.fam)! Please check the argument supplied to -p and try again! \n"
		exit 1
	    fi
	done
    else
	printf "Please supply the same number of subset files, labels, and final autosomal plink stems. Below are the supplied inputs. Please check, fix, and reattempt.\n"
	printf "Subset files: ${subset}\nSubset labels: ${lab}\nFinal autosomal files: ${autosomal_stem}\n"
	exit 1
    fi
fi

#validate the sex file
if [ ! -f "$sex_file" ]; 
then
    printf "File with sex information ($sex_file) does not exist! Please try again, specifying the correct input file.\n"
    exit 1
fi

if [ "$plink_memory_limit" ];
then 
    printf "Memory limit for plink calls: $plink_memory_limit \n"
    plink_memory_limit=$( echo "--memory $plink_memory_limit" )
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

#get the output folder
output_folder=${output_stem%/*}/
# get the stem of genotype files in case it is in another folder
geno_stem=${preimputation_geno_X##*/}

################# Start the post-imputation QC ######################


if [ $skip_first_filters = 'false' ];
then

    ####################### Step 1: Unzipping imputation results #######################
    if [ "$do_unzip" = 'true' ];
    then
	#unzip imputation results using password
	printf "\nStep 1 : Unzipping imputation results for chr X\n\n"
	unzip -P $( cat ${imputation_results_folder}/pass.txt ) ${imputation_results_folder}/chr_X.zip -d $imputation_results_folder
    else
	#make sure the imputation results are actually unzipped
	if test ! -f ${imputation_results_folder}/chrX.dose.vcf.gz ;
	then 
	    printf "\nThe -z was not specified but cannot find the imputation results (ie ${imputation_results_folder}/chrX.dose.vcf.gz)! Please check whether imputation results have been unzipped. If they haven't, make sure the decryption password is saved in ${imputation_results_folder}/pass.txt and specify the -z flag. \n"
	    exit 1
	else
	    printf "\nStep 1 (unzipping the imputation results) is already complete. Proceeding the step 2...\n\n"
	fi
    fi

    ####################### Step 2: Filter for R2 and multi-allelic variants, update SNP names, and merge chr together #######################
    printf "\nStep 2 : Filtering imputation results for R2<0.8 and multi-allelic variants\n"

    #restrict to variants with R2>=0.80    
    plink2 --vcf ${imputation_results_folder}/chrX.dose.vcf.gz --const-fid 0 --exclude-if-info "R2<0.8" $plink_memory_limit --make-bed --out ${output_stem}_chrX_temp > /dev/null;
    # exclude multi-allelic variants
    awk '{ print $2,$4 }' ${output_stem}_chrX_temp.bim | uniq -f1 -D | awk '{print $1}' >> ${output_stem}_chrX.dups
    plink2 --bfile ${output_stem}_chrX_temp --exclude ${output_stem}_chrX.dups $plink_memory_limit --make-bed --out ${output_stem}_chrX_temp_nodups > /dev/null;

    # if there are missing rs numbers, then update with the snp names file provided
    if [ -z "$( grep -m 1 'rs' ${output_stem}_chrX_temp_nodups.bim )" ];
    then
	# update name with rs numbers
	plink2 --bfile ${output_stem}_chrX_temp_nodups --update-name ${snp_names_file}_chrX.txt $plink_memory_limit --make-bed --out ${output_stem}_chrX_temp_nodups_names > /dev/null; 
        output=${output_stem}_chrX_temp_nodups_names
    else
	# if there are rs numbers are already present, then set this as the output stem and move on
        output=${output_stem}_chrX_temp_nodups
    fi

    #print out numbers of variants
    total_var=$(( $(grep "out of" ${output_stem}_chrX_temp.log | awk 'BEGIN { ORS="+" } { print $4 }' | sed 's/\(.*\)+/\1 /' ) ))
    afterR2=$( cat ${output_stem}_chrX_temp.bim | wc -l )
    nomulti=$( cat ${output_stem}_chrX_temp_nodups.bim | wc -l )
    printf "$total_var variants after imputation, $afterR2 variants with R2>0.8, and $nomulti variants with unique positions\n\n"
    

    samples=$( grep 'samples' ${output}.log | awk '{ print $1 }' )
    variants=$( grep 'variants' ${output}.log | awk '{ print $1 }' )
    printf "$samples samples and $variants variants\n\n"
else
    output=${output_stem}
    printf "\nSkipping steps 1 and 2 (the unzipping, conversion, and filtering of the individual chromosome) because the -x flag was supplied. Picking up at updating sample IDs and sex. Assuming the X chr file to start with is $output\n"
fi

############################ Step 3: Update SNP names and add sex to fam ######################

#update SNP names, sex, and perform standard SNP filtering
printf "\nStep 3: Updating IDs and adding sex back to the fam file.\n"
#update person ids using the race and sex file
output_last=$output
output=${output}_IDs
awk '{ print "0 "$1"_"$2" "$1" "$2 }' $sex_file > ${output_last}_update_ids.txt
plink --bfile $output_last --update-ids ${output_last}_update_ids.txt $plink_memory_limit --make-bed --out $output > /dev/null

#update SNP names and add sex back into fam file
output_last=$output
output=${output}_sex
plink --bfile ${output_last} --update-sex ${sex_file} $plink_memory_limit --make-bed --out ${output} > /dev/null

# log as a sanity check. Everyone should be updated. 
tmp=$( grep "people updated" ${output}.log | awk '{ print $2 }' )
printf "$tmp people whose ids and sex were able to be updated using the sex file.\n"


################################## subset to the current dataset ######################################3

printf "\nAll steps from this point on are specific to each subset. There will be ${#subset[@]} subsets.\n"

presubset_output=$output

for index in "${!subset[@]}" ; 
do 
    current_subset="${subset[index]}"
    current_lab="${lab[index]}"
    current_autosomal_stem="${autosomal_stem[index]}"

    printf "\nBeginning the rest of post-imputation QC for ${current_lab}.\nFile for initial subsetting: ${current_subset}\nPlink stem for final subsetting: ${current_autosomal_stem}\n"

    printf "\nStep 4 (${current_lab}): Subset to samples present in ${current_lab}\n"

    ### Imputed set
    output_last=${presubset_output}
    output=${presubset_output}_${current_lab}_geno05
    plink --bfile ${output_last} --keep ${current_subset} --geno 0.05 --make-bed --out ${output} $plink_memory_limit > /dev/null
    ## print notice if there are variants removed at this step
    if [[ $( grep "variants removed" ${output}.log | awk '{ print $1 }' ) -gt 0 ]]; then printf "Warning: There were variants removed for missingness from the imputed set when subsetting to ${current_lab}. This is unusual! Please check the imputed files!\n" ; fi
    ## report the number of samples and variants at this step
    samples=$( grep "pass filters and QC" ${output}.log | awk '{ print $4 }' )
    variants=$( grep "pass filters and QC" ${output}.log | awk '{ print $1 }' )
    printf "Subsetted to samples present in ${current_lab} in imputed set, leaving ${samples} samples and ${variants} variants.\n"

    ### Genotyped set
    output_geno=${output_folder}/${geno_stem}_${current_lab}_geno05
    plink --bfile ${preimputation_geno_X} --keep ${current_subset} --geno 0.05 --make-bed --out ${output_geno} $plink_memory_limit > /dev/null
    ## get the number of variants removed here for logging later. 
    geno=$( grep "variants removed" ${output_geno}.log | awk '{ print $1 }' )
    ## report the number of samples and variants at this step
    samples=$( grep "pass filters and QC" ${output_geno}.log | awk '{ print $4 }' )
    variants=$( grep "pass filters and QC" ${output_geno}.log | awk '{ print $1 }' )
    printf "Subsetted to samples present in ${current_lab} in genotyped set, leaving ${samples} samples and ${variants} variants.\n"


    ###### merge back in genotypes #####
    printf "\nStep 5 (${current_lab}): Merging back in the original genotypes.\n"

    # the genotyped variant ids in TOPMED imputation server format chrX:POS:REF:ALT
    # this step was performed in part2, however, keeping it for robustness
    printf "Updating pre-imputation variant IDs to TOPMED imputation server format chrX:POS:REF:ALT...\n"
    output_last_geno=$output_geno
    output_geno=${output_geno}_varID
    awk '{ print "chrX:"$4":"$6":"$5,$2}' ${output_last_geno}.bim > ${output_geno}.txt
    plink --bfile ${output_last_geno} --update-name ${output_geno}.txt 1 2 --make-bed $plink_memory_limit --out ${output_geno} > /dev/null

    # update to RSIDs
    printf "Now updating pre-imputation variant IDs to RSID...\n"
    output_last_geno=$output_geno
    output_geno=${output_geno}_rsid
    plink --bfile ${output_last_geno} --update-name ${snp_names_file}_chrX.txt --make-bed $plink_memory_limit --out ${output_geno} > /dev/null

    # the genotyped variant ids to exclude
    printf "Now removing those variants from the imputed fileset...\n"
    awk '{ print $2 }' ${output_geno}.bim > ${output_geno}_genotyped_variants.txt
    # remove variants in the imputation file for which there are genotypes
    output_last=$output
    output=${output}_nogeno
    plink --bfile ${output_last} --exclude ${output_geno}_genotyped_variants.txt --make-bed $plink_memory_limit --out $output > /dev/null
    
    # document the number removed                                                                              
    startvar=$( grep "variants loaded from" ${output}.log | awk '{ print $1 }')
    remaining_var=$( grep -e "variants remaining" ${output}.log | awk '{ print $2 }' )
    imputed_geno=$(( $startvar - $remaining_var ))
    printf "$remaining_var variants remaining after removing genotyped variants from imputation results.\n"

    # merge the genotyped and imputed data
    plink --bfile ${output} --bmerge ${output_geno} --make-bed $plink_memory_limit --out ${output}_merged > /dev/null
    ## calculate the number of genotyped (but not imputed) variants added                                      
    merged_var=$( grep "pass filters" ${output}_merged.log | awk '{ print $1 }' )
    geno_only=$(( $merged_var - $startvar ))

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
	printf "Getting same position warnings. Removing those variants from the imputed data and re-attempting merge.\n"
	# make a file with the variants at the same position
	grep "Warning: Variants" ${output}_merged.log | awk  '{ print $3"\n"$5 }' | sed -e "s/'//g" >  ${output_folder}/genotyped_variants_sameposwarnings.txt
	# remove those variants from imputed file
	plink --bfile ${output} --exclude ${output_folder}/genotyped_variants_sameposwarnings.txt --make-bed --out ${output}2 > /dev/null
	# get the number of variants removed
	startvar=$( grep "variants loaded from" ${output}2.log | awk '{ print $1 }')
	samepos=$( grep "variants remaining" ${output}2.log | awk '{ print $2 }' )
	removed_var=$(( $startvar - $samepos ))
	printf "${removed_var} variants at the same position in genotyped and imputed data. These variants were removed from the imputation results.\n"
	## re-attempting merge
	plink --bfile ${output}2 --bmerge ${output_geno} --make-bed --out ${output}_merged > /dev/null

	#check for more warnings
	sameposwarnings=$( grep "Warning: Variants" ${output}_merged.log | head -n1 )
	if [ ! -z "$sameposwarnings" ] ;
	then
	    printf "\nGetting more same position warnings! Please check and handle manually!\n"
	fi
    fi
    #update output variable
    output=${output}_merged
    printf "After merging in genotypes, there are $( grep "pass filters and QC" ${output}.log | sed 's/pass filters and QC.//' ).\n"
    printf "Output file: ${output} \n"


    ################################## SNP filters ###########################

    # SNP filters
    printf "\nStep 6 (${current_lab}): Applying Hardy-Weinberg equilibrium and MAF filters.\nHardy Weinberg equilibrium 1e6 based on females as only females are diploids XX\n"

    if [  "$n_sex" = 1 ] && [ "${sex}" = 1 ];
    then
	    printf "\nSkipping Hardy Weinberg equilibrium 1e6 based on females: only males present \n"
    else

	    output_last=${output}
	    output_fem=${output}_1e6femHWE
	    output=${output}_femHWE6

	    # get list of SNPS to subset to based on HWE test in females
	    plink --bfile ${output_last} --filter-females --hwe 0.000001 $plink_memory_limit --make-bed --out ${output_fem}  > /dev/null
	    awk '{print $2}' ${output_fem}.bim > ${output_fem}_SNPs.txt
	    # extract those variants in the whole sample
	    plink --bfile ${output_last} --extract ${output_fem}_SNPs.txt $plink_memory_limit --make-bed --out ${output} > /dev/null
	    # document, but don't print yet
	    hwe=$(( "$( wc -l < ${output_last}.bim )" - "$( wc -l < ${output_fem}.bim )" ))
	    variants=$( grep "pass filters and QC" ${output}.log | awk '{print $1;}' )
	    samples=$( grep "pass filters and QC" ${output}.log | awk '{print $4;}' )
    fi

    
    output_last=$output
    output=${output}_maf01
    plink --bfile ${output_last} --maf 0.01 $plink_memory_limit --make-bed --out ${output} > /dev/null
    # document
    maf=$( grep "removed due to minor allele threshold(s)" ${output}.log | awk '{print $1;}' )
    variants=$( grep "pass filters and QC" ${output}.log | awk '{print $1;}' )
    samples=$( grep "pass filters and QC" ${output}.log | awk '{print $4;}' )

    ## now print out numbers since variant exclusions are done
    printf '%s\n' "Removed $geno genotyped SNPs for >5% missingness, merged in genotyped SNPs (-$removed_var same position, +$geno_only genotyped-only), removed $hwe SNPs for HWE p<1e-6 and $maf SNPs for MAF <0.01."
    printf "$samples samples\n$variants variants\n"

    # check no hh warnings on females and set hh to missing
    printf "\nConfirm no hh warning for females only below: 'Warning: ... het. haploid genotypes present' \n"
    # will have issues with males only cohort. skip this for male only
    plink --bfile ${output} --filter-females --freq $plink_memory_limit --out $output > /dev/null

    # set hh missing
    output_last=$output
    output=${output}_nohh
    printf "\nSetting het. haploid genotypes missing...\n"
    plink --bfile ${output_last} --set-hh-missing $plink_memory_limit --make-bed --out ${output} > /dev/null

    
    ##### Step 7: Remove individuals based on final autosomal files #####
    ## first check whether anyone will be removed
    if [[ $( wc -l < ${output}.fam ) -eq $( wc -l < ${current_autosomal_stem}.fam ) ]]; 
    then
	printf "\nSkipping Step 7 (${current_lab}) because no samples appear to have been dropped in the final autosomal file for this set.\n"
	#log the final numbers and file name
	final_samples=$( grep "pass filters and QC" ${output}.log | awk '{ print $4 }' )
	variants=$( grep "pass filters and QC" ${output}.log | awk '{ print $1 }' )
	printf "$samples samples\n$variants variants\n"
	printf "\nFinal X-chromosome plinkset for ${current_lab[index]}: ${output}\n"


    else
	## if there are a different number of samples, run the exclusion
	printf "\nStep 7 (${current_lab}): Remove individuals based on autosomal PCs outlier, restricting to final plinkset provided  \n"

	output_last=${output}
	output=${output}_noout
	plink --bfile ${output_last} --keep ${current_autosomal_stem}.fam --make-bed --out ${output} $plink_memory_limit > /dev/null
	final_samples=$( grep "pass filters and QC" ${output}.log | awk '{ print $4 }' )
	variants=$( grep "pass filters and QC" ${output}.log | awk '{ print $1 }' )

	# calculate number dropped based on final autosomal files
	dropped_samples=$(( $samples - $final_samples ))
	printf "Removed $dropped_samples samples based on autosomal heterozygosity and PC values.\n"

	printf "$final_samples samples\n$variants variants\n"
	printf "\nFinal X-chromosome plinkset for ${current_lab[index]}: ${output}\n"
    fi

done

printf "\nX chromosome post-imputation QC complete!\n"
