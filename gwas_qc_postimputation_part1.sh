
#!/bin/bash

#fail on error
set -e

#define usage
display_usage() {
    printf "GWAS QC post-imputation script, part 1

This script will unzip the imputation results (assuming password is saved in pass.txt in the same folder as the imputation results files and will perform the first steps of standard post-imputation QC for our common variant pipeline. This includes filtering for R2, removing multi-allelic variants and filtering out variants for low MAF. Then, ancestry and PCs will be estimated using pre-calculated weights. The output of this script will include filtered genotypes, plots of predicted PCs, and a file including ancestral group estimations. Please check these outputs and proceed with the part 2 post-imputation script.

Usage:
SCRIPTNAME.sh -o [output_stem] -i [imputation_results_folder] -f [sex_file] -s [snp_names_file] -w [snpweights_file] -z -x -m [plink_memory_flag] -c

output_stem = the beginning part of all QC'ed files including the full path to the folder in which they should be created

imputation_results_folder = the folder to which the imputation results have been downloaded which also contains a file called pass.txt with the password to unzip the imputation results files

sex_file = a file with FID and IID (corresponding to the fam file), 1 column indicating sex for the sex check (1 for males, 2 for females, 0 if unknown), with NO header.

snp_names_file = the file stem for converting the SNP names from imputation results to rs numbers. There should be one for each chromosome and each must have 2 columns: imputation result SNP ids and rs numbers. Can have header but it will be ignored.

snpweights_file = the full path and name of the file containing pre-calculated weights from which to calculate ancestral estimates; note that there must be a corresponding file with the same stem ending in *_snps.txt which contains a list of variants to subset the current set to before calculating weights

plink_memory_limit (optional) = argument indicating the memory limit for plink to use rather than the default of half the RAM. This is useful for running this step of QC locally.

-z indicates that the imputation results will need to be unzipped.
-x indicates to skip the first R2<0.8 and multi-allelic variant filtering of the individual chr files. This comes in handy if there were an issue with the next step(s) because this first step is the longest.
-c will skip the clean-up at the end of the script which removes intermediate *.bed files.
-h will display this message
"
        }

#parse options
do_unzip='false'
skip_first_filters='false'
skip_cleanup='false'
while getopts 'o:i:f:s:w:zxcm:h' flag; do
  case "${flag}" in
    o) output_stem="${OPTARG}" ;;
    i) imputation_results_folder="${OPTARG}" ;;
    f) sex_file="${OPTARG}" ;;
    s) snp_names_file="${OPTARG}" ;;
    z) do_unzip='true' ;;
    w) snpweights_file="${OPTARG}" ;;
    x) skip_first_filters='true' ;;
    c) skip_cleanup='true' ;;
    m) plink_memory_limit="${OPTARG}";;
    h) display_usage ; exit ;;
    \?|*) display_usage
       exit 1;;
  esac
done

########################################## Input validation ##################################################

#check to make sure necessary arguments are present
if [ -z "$output_stem" ] || [ -z "$imputation_results_folder" ] || [ -z "$sex_file" ] || [ -z "$snp_names_file" ] || [ -z "$snpweights_file" ];
then
    printf "Error: Necessary arguments not present!\n\n"
    display_usage
    exit 1
fi

#print out inputs
printf "GWAS QC Post-imputation Script Part 1

Output file path and stem for cleaned imputed files : $output_stem
Imputation results folder : $imputation_results_folder
Sex information file : $sex_file
Stem for files for SNP name conversion : $snp_names_file
"
if [ "$do_unzip" = 'false' ];
then 
    printf "The -z was not specified so imputation results must already be unzipped.\n"
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
output_folder=${output_stem%/*}


############################## Step 1: Unzipping, filtering, and merging imputation results (R2 and multi-allelic) #############################


if [ $skip_first_filters = 'false' ];
then
    
    # first, unzip the imputation results if needed. 
    if [ "$do_unzip" = 'true' ];
    then

	####################### Step 1: Unzipping imputation results #######################
	printf "\nStep 1 : Unzipping imputation results\n\n"

	#check to make sure password file is present
	if [ ! -f "${imputation_results_folder}/pass.txt" ];
	then
	    printf "Imputation password file not found! Please save the password to unzip imputation results (found in the email notification from the imputation server of the job's completion) to ${imputation_results_folder}/pass.txt\n"
	    exit 1
	fi

	#unzip imputation results using password
	for i in $( ls ${imputation_results_folder}/*.zip );
	do
	    unzip -P $( cat ${imputation_results_folder}/pass.txt ) $i -d $imputation_results_folder
	done
    else
	#make sure the imputation results are actually unzipped
	if test ! -f ${imputation_results_folder}/chr1.dose.vcf.gz ;
	then 
	    printf "The -z was not specified but cannot find the imputation results (ie ${imputation_results_folder}/chr1.dose.vcf.gz for each chromosome )! Please check whether imputation results have been unzipped. If they haven't, make sure the decryption password is saved in ${imputation_results_folder}/pass.txt and specify the -z flag.\n"
	    exit 1
	else
	    printf "\nStep 1 (unzipping the imputation results) is already complete. Proceeding the step 2...\n"
	fi
    fi

    ####################### Step 2: Filter for R2 and multi-allelic variants, update SNP names, and merge chr together #######################
    # before the loop, initialize the file for merging each chr together
    echo "" | tee ${output_stem}_merge_list.txt

    #filter imputation results for R2<0.8 and remove multi-allelic variants (multiple rows in vcf->bim)
    printf "\nStep 2 : Filtering imputation results for R2<0.8 and multi-allelic variants, updating variant names, and merging individual chr files together\n"
    for i in $(seq 1 22); do
	#restrict to variants with R2>=0.80
	plink2 --vcf ${imputation_results_folder}/chr${i}.dose.vcf.gz --const-fid 0 --exclude-if-info "R2<0.8" --make-bed --out ${output_stem}_chr${i}_temp $plink_memory_limit > /dev/null
	
	#get list of variants which have the same position (by base pair since this is just 1 chromosome) and remove
	dup_pos=$( awk '{ print $4 }' ${output_stem}_chr${i}_temp.bim | uniq -d )
	for j in $dup_pos ; do grep -w $j ${output_stem}_chr${i}_temp.bim | awk '{ print $2 }' >> ${output_stem}_chr${i}.dups ; done
	plink2 --bfile ${output_stem}_chr${i}_temp --exclude ${output_stem}_chr${i}.dups --make-bed --out ${output_stem}_chr${i}_temp_nodups $plink_memory_limit > /dev/null
	
	#update variant names with rs numbers
	plink2 --bfile ${output_stem}_chr${i}_temp_nodups --update-name ${snp_names_file}_chr${i}.txt --make-bed --out ${output_stem}_chr${i}_temp_nodups_names $plink_memory_limit > /dev/null
	
	#check the length of variant IDs
	awk 'length($2)>30{ print $2" "$1"_"$4 }' ${output_stem}_chr${i}_temp_nodups_names.bim > ${output_stem}_chr${i}_temp_nodups_names_shortids.txt
	# if there are IDs to update, then update them
	# either way, add the appropriate file name to the merge-list
	if [ $(ls ${output_stem}_chr${i}_temp_nodups_names_shortids.txt) ] ;
	then
	    plink --bfile ${output_stem}_chr${i}_temp_nodups_names --update-name ${output_stem}_chr${i}_temp_nodups_names_shortids.txt --make-bed --out ${output_stem}_chr${i}_temp_nodups_names_shorterids $plink_memory_limit > /dev/null ; 
	    printf "Updating $( grep 'updated' ${output_stem}_chr${i}_temp_nodups_names_shorterids.log | awk '{ print $2 }') too-long IDs from chr ${i}...\n" ; 
	    printf "${output_stem}_chr${i}_temp_nodups_names_shorterids\n" >> ${output_stem}_merge_list.txt
	else
	    printf "IDs in chr ${i} ok.\n"
	    printf "${output_stem}_chr${i}_temp_nodups_names\n" >> ${output_stem}_merge_list.txt ;
	fi

    done

    #print out numbers of variants
    total_var=$(( $(grep "out of" ${output_stem}_chr*_temp.log | awk 'BEGIN { ORS="+" } { print $4 }' | sed 's/\(.*\)+/\1 /' ) ))
    afterR2=$(( $total_var - $( cat ${output_stem}_chr*_temp.bim | wc -l ) ))
    nomulti=$(( $total_var - $afterR2 - $( cat ${output_stem}_chr*_temp_nodups.bim | wc -l ) ))
    printf "$total_var variants after imputation, $afterR2 variants removed for R2<0.8, and $nomulti duplicated variants removed.\n\n"

    #Merge all chromosomes into one file
    printf "Merging all chromosomes into one file\n"

    #merge individual chromosome files to be one large file
    output=${output_stem}
    plink --merge-list ${output_stem}_merge_list.txt --make-bed --out $output $plink_memory_limit > /dev/null
    grep 'pass filters and QC' ${output}.log
else
    output=$output_stem
    # check for the presence of the merged file with the supplied name and fail if not present
    if [ ! -f ${output_stem}.bed ]; then exit 1 ; fi
    # if present, log skipping the first step
    printf "\nSkipping steps 1 and 2 (the unzipping, conversion, and filtering of the individual chromosomes) because the -x flag was supplied! Picking up at updating sample IDs and sex in the merged file. Assuming the merged file stem is $output\n"
fi

################################################# Step 3: Update SNP names and add sex to fam ##########################################

printf "\nStep 3: Updating sample and variant names and adding sex back to the fam file.\n"

#### Updating person ids and sex ####
#update person ids using the sex file
output_last=$output
output=${output}_IDs
awk '{ print "0 "$1"_"$2" "$1" "$2 }' $sex_file > ${output_last}_update_ids.txt
plink --bfile $output_last --update-ids ${output_last}_update_ids.txt --make-bed --out $output $plink_memory_limit > /dev/null

#update SNP names and add sex back into fam file
output_last=$output
output=${output}_sex
plink --bfile ${output_last} --update-sex ${sex_file} --make-bed --out ${output} $plink_memory_limit > /dev/null

tmp=$( grep "people updated" ${output}.log | awk '{ print $2 }' )
printf "$tmp people whose ids and sex were able to be updated using the sex file.\n"

################################################# Step 4: Apply SNPWeights to determine genetic ancestry ##########################################

# Calculate SNPWeights
printf "\nStep 4: Calculate inferred ancestry.\n\nFirst, subsetting to overlap with pre-calculated weights and removing non-ATGC variants...\n"
## Filter down to the variants in 1000G
plink --bfile ${output} --extract ${snpweights_file}_snp.txt --make-bed --out ${output}_overlap $plink_memory_limit > /dev/null 
grep -e 'pass filters and QC' ${output}.log

printf "\nNow, remove non-ATGC variants...\n"
## remove non-ATGC variants
plink --bfile ${output}_overlap --snps-only 'just-acgt' --make-bed --out ${output}_overlap_acgt $plink_memory_limit > /dev/null
grep -e 'pass filters and QC' ${output}.log

## Prune 
printf "Pruning...\n"
plink --bfile ${output}_overlap_acgt --indep-pairwise 200 100 0.2 --allow-no-sex --out ${output}_overlap_acgt_prune $plink_memory_limit > /dev/null
plink --bfile ${output}_overlap_acgt --output-missing-phenotype 1 --extract ${output}_overlap_acgt_prune.prune.in --make-bed --out ${output}_overlap_acgt_pruned $plink_memory_limit > /dev/null
rm ${output}_overlap_acgt_prune.*
printf "$( wc -l < ${output}_overlap_acgt_pruned.bim ) variants out of $( wc -l < ${output}_overlap_acgt.bim ) left after pruning.\n"

## Convert into EIG format
printf "Now, converting to EIG format for SNPWeights calculation...\n"
printf "genotypename: ${output}_overlap_acgt_pruned.bed
snpname: ${output}_overlap_acgt_pruned.bim
indivname: ${output}_overlap_acgt_pruned.fam
outputformat: EIGENSTRAT
genotypeoutname: ${output}_overlap_acgt_pruned.eigenstratgeno
snpoutname: ${output}_overlap_acgt_pruned.snp
indivoutname: ${output}_overlap_acgt_pruned.ind
familynames: YES
" > ${output}_overlap_acgt_pruned_convert_to_EIG.par

convertf -p ${output}_overlap_acgt_pruned_convert_to_EIG.par > /dev/null

## run ancestry inference
printf "Now, inferring genetic ancestry using SNPWeights...\n"
printf "geno: ${output}_overlap_acgt_pruned.eigenstratgeno
snp: ${output}_overlap_acgt_pruned.snp
ind: ${output}_overlap_acgt_pruned.ind
snpwt: ${snpweights_file}
predpcoutput: ${output}_overlap_acgt_pruned.out
" > ${output}_overlap_acgt_pruned_infer_ancestry.par 

python /SNPweights2.1/inferancestry.py --par ${output}_overlap_acgt_pruned_infer_ancestry.par > /dev/null

# Separate samples into ancestral categories based on standard thresholds
printf "Finally, using SNPWeights output to assign ancestral groupings...\n"
Rscript assign_ancestral_categories.R ${output}_overlap_acgt_pruned.out

# generate hash file listing sets to be carried forward (only those which have more than 10 samples)
wc -l ${output}_overlap_acgt_pruned*_keep.txt | awk '{ if($1>10 && $2!="total") print }' | tee ${output}_ancestral_groups.txt

################################################# Cleanup ##########################################

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
