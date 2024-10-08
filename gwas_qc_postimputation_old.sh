#/bin/bash

#fail on error
set -e

#define usage
display_usage() {
    printf "GWAS QC post-imputation script

This script will unzip the imputation results (assuming password is saved in pass.txt in the same folder as the imputation results files and will perform standard post-imputation QC for our common variant pipeline. This includes filtering for R2, removing multi-allelic variants and filtering out variants for low MAF or HWE disequilibrium. Finally, PCs will be calculated on the final file-set.

Usage:
SCRIPTNAME.sh -o [output_stem] -i [imputation_results_folder] -r [race_file] -f [sex_file] -s [snp_names_file] -g [preimputation_geno] -z -x -c

output_stem = the beginning part of all QC'ed files including the full path to the folder in which they should be created

imputation_results_folder = the folder to which the imputation results have been downloaded which also contains a file called pass.txt with the password to unzip the imputation results files

sex_file = a file with FID and IID (corresponding to the fam file), 1 column indicating sex for the sex check (1 for males, 2 for females, 0 if unknown), with NO header.

race_file (optional) = a file with FID and IID (corresponding to the fam file) and 1 column indicating both race and ethnicity for PC plots, with NO header. Non-hispanic whites need to be indicated with 'White' or 'EUR.' No other values in the race column must be fixed; however, the race column must not include spaces. This is only needed if you want to color PC plots based on race (which at this stage will only be self-report). Ancestral categories will be calculated post-imputation using SNPWeights. 

snp_names_file = the file stem for converting the SNP names from imputation results to rs numbers. There should be one for each chromosome and each must have 2 columns: imputation result SNP ids and rs numbers. Can have header but it will be ignored.

preimputation_geno = the full path and stem to the cleaned final pre-imputation files to be merged back into the final files

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
while getopts 'o:i:r:f:g:zxch' flag; do
  case "${flag}" in
    o) output_stem="${OPTARG}" ;;
    i) imputation_results_folder="${OPTARG}" ;;
    r) race_file="${OPTARG}" ;;
    f) sex_file="${OPTARG}" ;;
    s) snp_names_file="${OPTARG}" ;;
    z) do_unzip='true' ;;
    g) preimputation_geno="${OPTARG}" ;;
    x) skip_first_filters='true' ;;
    c) skip_cleanup='true' ;;
    h) display_usage ; exit ;;
    \?|*) display_usage
       exit 1;;
  esac
done

#check to make sure necessary arguments are present                                                                                              
if [ -z "$output_stem" ] || [ -z "$imputation_results_folder" ] || [ -z "$race_file" ] || [ -z "$snp_names_file" ] || [ -z "$preimputation_geno" ];
then
    printf "Error: Necessary arguments not present!\n\n"
    display_usage
    exit 1
fi

#print out inputs
printf "GWAS QC Post-imputation Script

Output file path and stem for cleaned imputed files : $output_stem
Imputation results folder : $imputation_results_folder
Sex information file : $sex_file
Race information file : $race_file
Stem for files for SNP name conversion : $snp_names_file
"
if [ "$do_unzip" = 'false' ];
then 
    printf "The -z was not specified so imputation results must already be unzipped.\n\n"
fi

#validate the genotyped files
if [ ! -f "${preimputation_geno}.fam" ];
then
    printf "Cannot see the preimputation genotype files ($preimputation_geno)! Please check the argument supplied to -g and try again!\n"
    exit 1
fi

#validate the sex file
if [ ! -f "$sex_file" ]; 
then
    printf "File with sex information ($sex_file) does not exist! Please try again, specifying the correct input file.\n"
    exit 1
fi

#validate the race file
if [ ! -f "$race_file" ]; 
then
    printf "File with race information ($race_file) does not exist! Please try again, specifying the correct input file.\n"
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

#get the output folder
output_folder=${output_stem%/*}


################# Start the post-imputation QC ######################

if [ "$do_unzip" = 'true' ];
then
    printf "Step 1 : Unzipping imputation results\n\n"

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
	printf "Step 1 (unzipping the imputation results) is already complete. Proceeding the step 2...\n\n"
    fi
fi

if [ $skip_first_filters = 'false' ];
then
    #filter imputation results for R2<0.8 and remove multi-allelic variants (multiple rows in vcf->bim)
    printf "Step 2 : Filtering imputation results for R2<0.8 and multi-allelic variants\n"
    for i in $(seq 1 22); do
	#restrict to variants with R2>=0.80
	plink2 --vcf ${imputation_results_folder}/chr${i}.dose.vcf.gz --const-fid 0 --exclude-if-info "R2<0.8" --make-bed --out ${output_stem}_chr${i}_temp > /dev/null
	
	#get list of variants which have the same position (by base pair since this is just 1 chromosome) and remove
	dup_pos=$( awk '{ print $4 }' ${output_stem}_chr${i}_temp.bim | uniq -d )
	for j in $dup_pos ; do grep -w $j ${output_stem}_chr${i}_temp.bim | awk '{ print $2 }' >> ${output_stem}_chr${i}.dups ; done
	plink2 --bfile ${output_stem}_chr${i}_temp --exclude ${output_stem}_chr${i}.dups --make-bed --out ${output_stem}_chr${i}_temp_nodups > /dev/null
	
	#update variant names with rs numbers
	plink2 --bfile ${output_stem}_chr${i}_temp_nodups --update-name ${snp_names_file}_chr${i}.txt --make-bed --out ${output_stem}_chr${i}_temp_nodups_names > /dev/null
    done

    #print out numbers of variants
    total_var=$(( $(grep "out of" ${output_stem}_chr*_temp.log | awk 'BEGIN { ORS="+" } { print $4 }' | sed 's/\(.*\)+/\1 /' ) ))
    afterR2=$(( $total_var - $( cat ${output_stem}_chr*_temp.bim | wc -l ) ))
    nomulti=$(( $total_var - $afterR2 - $( cat ${output_stem}_chr*_temp_nodups.bim | wc -l ) ))
    printf "$total_var variants after imputation, $afterR2 variants removed for R2<0.8, and $nomulti duplicated variants removed.\n\n"

    #Merge all chromosomes into one file
    printf "Step 3 : Merging all chromosomes into one file\n"
    #create merge file (emptying old version if present)
    echo "" | tee ${output_stem}_merge_list.txt
    for i in $( seq 1 22 );
    do
	printf "${output_stem}_chr${i}_temp_nodups_names\n" >> ${output_stem}_merge_list.txt ;
    done
    
    #merge individual chromosome files to be one large file
    output=${output_stem}
    plink --merge-list ${output_stem}_merge_list.txt --make-bed --out $output > /dev/null
    grep 'pass filters and QC' ${output}.log
else
    output=$output_stem
    printf "Skipping the conversion and filtering of the individual chromosome because the -x flag was supplied! Picking up at updating sample IDs and sex in the merged file. Assuming the merged file stem is $output\n"
fi

#### Updating person ids and sex ####
#update person ids using the sex file
output_last=$output
output=${output}_IDs
awk '{ print "0 "$1"_"$2" "$1" "$2 }' $sex_file > ${output_last}_update_ids.txt
plink --bfile $output_last --update-ids ${output_last}_update_ids.txt --make-bed --out $output > /dev/null
tmp=$( grep "people updated" ${output}.log | awk '{ print $2 }' )
printf "$tmp people whose ids were able to be updated using the sex file.\n"

#update SNP names and add sex back into fam file
output_last=$output
output=${output}_sex
plink --bfile ${output_last} --update-sex ${sex_file} --make-bed --out ${output} > /dev/null

###### merge back in genotypes #####

printf "\nStep 4: Merging back in the original genotypes.\n"

# make file with all the genotyped variant ids -- this will be problematic if we don't have separate folders for NHW/all-races
awk '{ print $2 }' ${preimputation_geno}.bim > ${output_folder}/genotyped_variants.txt

# remove variants for which there are genotypes from the bim file
output_last=$output
output=${output}_nogeno
plink --bfile ${output_last} --exclude ${output_folder}/genotyped_variants.txt --make-bed --out $output > /dev/null
printf "$(grep -e "variants remaining" ${output}.log | awk '{ print $2 }') variants remaining after removing genotyped variants from imputation results.\n"

# merge the genotyped and imputed data
plink --bfile ${output} --bmerge ${preimputation_geno} --make-bed --out ${output}_merged > /dev/null

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
    grep "Warning: Variants" ${output}_merged.log | awk  '{ print $3"\n"$5 }' | sed -e "s/'//g" >  ${output_folder}/genotyped_variants_sameposwarnings.txt
    plink --bfile ${output} --exclude ${output_folder}/genotyped_variants_sameposwarnings.txt --make-bed --out ${output}2 > /dev/null
    printf "$(grep -e "variants remaining" ${output}.log | awk '{ print $2 }' ) variants remaining after removing genotyped variants from imputation results based on position.\n"
    plink --bfile ${output}2 --bmerge ${preimputation_geno} --make-bed --out ${output}_merged > /dev/null
    
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
printf "Output file: $output"

printf "\nStep 5: Applying standard variant filters\n"

# SNP filters
output_last=$output
output=${output}_maf01_hwe6
plink --bfile ${output_last} --maf 0.01 --hwe 0.000001 --make-bed --out ${output} > /dev/null
grep -e "hwe: " -e "removed due to minor allele threshold" -e 'pass filters and QC' ${output}.log
printf "Output file: $output"

# prune for heterozygosity check
printf "\n\nStep 6: Pruning and running heterozygosity check\n"
plink --bfile ${output} --indep-pairwise 200 100 0.2 --allow-no-sex --out ${output}_prune > /dev/null
plink --bfile ${output} --output-missing-phenotype 1 --extract ${output}_prune.prune.in --make-bed --out ${output}_pruned > /dev/null
rm ${output}_prune.*
printf "$( wc -l < ${output}_pruned.bim ) variants out of $( wc -l < ${output}.bim ) left after pruning.\n"

##### heterozygosity check #####
plink --bfile ${output}_pruned --het --out ${output}_pruned_hetcheck > /dev/null
#plot and check for outliers
Rscript plot_het_check_outliers.R ${output}_pruned_hetcheck

#if there are outliers >6 sd from the F stat mean, remove them
if [ -f "${output}_pruned_hetcheck_outliers.txt" ];
then
    output_last=${output}
    output=${output}_nohetout
    plink --bfile $output_last --remove ${output}_pruned_hetcheck_outliers.txt --make-bed --out ${output} > /dev/null
    grep -e ' people pass filters and QC' ${output}.log
    printf "Output file: $output \n"
fi

##### PC calculation ####
printf "\nStep 7: Calculating post-imputation PCs\n\n"

# Calculate PCs
sh calc_plot_PCs.sh -i $output -r $race_file

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

printf "\nPlease check PC plots for outliers, remove if present, and recalculate PCs. If there are no outliers to remove, GWAS QC for this dataset is (probably) complete!\n"


