#!/bin/bash
# VJ: 20200709 modified post imputation script for single chromosome run
# gwas_qc_postimputation_chr.sh
# added -c option to specify single chromosome number[int 1-22],[X/23]
#fail on error
set -e

# VJ: singularity image and directories bound for reproducibility
[ ! -z "$SINGULARITY_CONTAINER" ] && printf "\nSingularity image: $SINGULARITY_CONTAINER"
[ ! -z "$SINGULARITY_COMMAND" ] && printf "\nSingularity shell: $SINGULARITY_COMMAND"
[ ! -z "$SINGULARITY_BIND" ] && printf "\nMapped directories:\n $SINGULARITY_BIND\n" | sed 's/,/\n /g'

#define usage
display_usage() {
    printf "GWAS QC post-imputation script

This script will unzip the imputation results for single chromosome specified(assuming password is saved in pass.txt in the same folder as the imputation results files and will perform standard post-imputation QC for our common variant pipeline. This includes filtering for R2, removing multi-allelic variants and filtering out variants for low MAF or HWE disequilibrium. Finally, PCs will be calculated on the final file-set.

Usage:
SCRIPTNAME.sh -c [CHR] -o [output_stem] -i [imputation_results_folder] -r [race_file] -f [sex_file] -s [snp_names_file] -g [preimputation_geno] -p [final autosomal genotype plinkset] -z -x -d

CHR = X; Chromosome number [int 1-26] to be processed

output_stem = the beginning part of all QC'ed files including the full path to the folder in which they should be created

imputation_results_folder = the folder to which the imputation results have been downloaded which also contains a file called pass.txt with the password to unzip the imputation results files

sex_file = a file with FID and IID (corresponding to the fam file), 1 column indicating sex for the sex check (1 for males, 2 for females, 0 if unknown), with NO header.

race_file (optional) = a file with FID and IID (corresponding to the fam file) and 1 column indicating both race and ethnicity for PC plots, with NO header. Non-hispanic whites need to be indicated with 'White' or 'EUR.' No other values in the race column must be fixed; however, the race column must not include spaces. This is only needed if you want to color PC plots based on race (which at this stage will only be self-report). Ancestral categories will be calculated post-imputation using SNPWeights. 

snp_names_file = the file stem for converting the SNP names from imputation results to rs numbers. There should be one for each chromosome and each must have 2 columns: imputation result SNP ids and rs numbers. Can have header but it will be ignored.

preimputation_geno = the full path and stem to the cleaned final pre-imputation files to be merged back into the final files

-p = final autosomal plinkset stem with IID and FID matching the input plinkset and PC outliers removed, used in Xchromosome GWAS QC to match autosomal individuals in final plinkset.

-z indicates that the imputation results will need to be unzipped. All *.zip files in imputation_results_folder will be unzipped with password provided in pass.txt
-x indicates to skip the first filtering of the individual chr files. This comes in handy if there were an issue with the next step(s) because this first step is the longest.
-d will skip the clean-up at the end of the script which removes intermediate *.bed files.
-h will display this message


Brief summary:
Step1: unzipping the imputation results 
	-z: all *.zip files in imputation_results_folder will be unzipped
Step2: Filtering imputation results for R2<0.8 and multi-allelic variants
	filter vcf got R2>0.8
Step3: Merging all chromosomes(if more than 1) into one plinkset
	Merging back in the original genotypes.
        remove multi-position variants(assumed multi-allelic)
        update VarID with RSID from ref file provided
Step4: update IDs,SNP names and sex
        IDs: 0 FID_IID -> FID IID
        sex: from race and sex file
        SNP names: from ref file
Step5:  merge back in genotyped variants
        the genotyped variant ids in TOPMED imputation server format chrX:POS:REF:ALT
        update genotyped variants to RSIDs
        merge the genotyped and imputed data
Step6: SNP filters
        HWE(based on females): 1e-6
        MAF: 0.01 (1 percent)
Step7: (skipped performed in part1) Heterozygocity + Autosmal PC outliers 
Step8: keep SNPs overlapping in males and females (skip for cohorts with single sex) 
Step9: keep individuals in final autosomal plinkset \n\n"
        }

#parse options
do_unzip='false'
skip_first_filters='false'
skip_cleanup='false'
while getopts 'c:o:i:r:f:s:g:p:zxdh' flag; do
  case "${flag}" in
    c) CHR="${OPTARG}";; # VJ: TODO include test for CHR -lt 26 and -gt 0
    o) output_stem="${OPTARG}" ;;
    i) imputation_results_folder="${OPTARG}" ;;
    r) race_file="${OPTARG}" ;;
    f) sex_file="${OPTARG}" ;;
    s) snp_names_file="${OPTARG}" ;;\
    z) do_unzip='true' ;;
    g) preimputation_geno="${OPTARG}" ;;
    p) ds_plinkset="${OPTARG}" ;;
    x) skip_first_filters='true' ;;
    d) skip_cleanup='true' ;;
    h) display_usage ; exit ;;
    \?|*) display_usage
       exit 1;;
  esac
done

printf "\n Brief summary:
Step1: unzipping the imputation results
        -z: all *.zip files in imputation_results_folder will be unzipped
Step2: Filtering imputation results for R2<0.8 and multi-allelic variants
        filter vcf got R2>0.8
        remove multi-position variants(assumed multi-allelic)
        update VarID with RSID from ref file provided
Step3: Merging all chromosomes(if more than 1) into one plinkset
Step4: update IDs,SNP names and sex
        IDs: 0 FID_IID -> FID IID
        sex: from race and sex file
        SNP names: from ref file
Step5:  merge back in genotyped variants
	the genotyped variant ids in TOPMED imputation server format chrX:POS:REF:ALT
	update genotyped variants to RSIDs
	merge the genotyped and imputed data
Step6: SNP filters
        HWE(based on females): 1e-6
        MAF: 0.01 (1 percent)
Step7: (skipped performed in part1) Heterozygocity + Autosmal PC outliers 
Step8: keep SNPs overlapping in males and females (skip for single sex cohorts)
Step9: keep individuals in final autosomal plinkset \n\n"

#check to make sure necessary arguments are present                                                                                                    
if [ -z "$output_stem" ] || [ -z "$imputation_results_folder" ] || [ -z "$sex_file" ] || [ -z "$race_file" ] || [ -z "${snp_names_file}" ] ;
then
    printf "Error: Necessary arguments not present! \n\n"
    display_usage
    exit 1
fi

#print out inputs
printf "GWAS QC Post-imputation Script

Output file path and stem for cleaned imputed files : $output_stem
Imputation results folder : ${imputation_results_folder}
Sex information file : ${sex_file}
Race information file : ${race_file}
Stem for files for SNP name conversion : ${snp_names_file}_chr${CHR}
"
if [ "$do_unzip" = 'false' ];
then 
    printf "The -z was not specified so imputation results must already be unzipped.\n\n"
fi

#validate the genotyped files
if [ ! -f "${preimputation_geno}.fam" ];
then
    printf "Cannot see the preimputation genotype files ($preimputation_geno)! Please check the argument supplied to -g and try again! \n"
    exit 1
fi

#validate the final plinkset files(only .fam needed)
if [ ! -f "${ds_plinkset}.fam" ];
then
    printf "Cannot see the preimputation genotype files (${ds_plinkset}.fam)! Please check the argument supplied to -p and try again! \n"
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
output_folder=${output_stem%/*}/


################# Start the post-imputation QC ######################

#VJ: still unzipping all the files present
if [ "$do_unzip" = 'true' ];
then
    #unzip imputation results using password
    printf "Step 1 : Unzipping imputation results\n\n"
    unzip -P $( cat ${imputation_results_folder}/pass.txt ) ${imputation_results_folder}/chr_${CHR}.zip -d $imputation_results_folder
else
    #make sure the imputation results are actually unzipped
    if test ! -f ${imputation_results_folder}/chr${CHR}.dose.vcf.gz ;
    then 
	printf "The -z was not specified but cannot find the imputation results (ie ${imputation_results_folder}/chr${CHR}.dose.vcf.gz)! Please check whether imputation results have been unzipped. If they haven't, make sure the decryption password is saved in ${imputation_results_folder}/pass.txt and specify the -z flag. \n"
	exit 1
    else
	printf "Step 1 (unzipping the imputation results) is already complete. Proceeding the step 2...\n\n"
    fi
fi

if [ $skip_first_filters = 'false' ];
then
	#filter imputation results for R2<0.8 and remove multi-allelic variants (multiple rows in vcf->bim)
	printf "Step 2 : Filtering imputation results for R2<0.8 and multi-allelic variants\n"
	#restrict to variants with R2>=0.80    
	plink2 --vcf ${imputation_results_folder}/chr${CHR}.dose.vcf.gz --const-fid 0 --exclude-if-info "R2<0.8" --memory 15000 --make-bed --out ${output_stem}_chr${CHR}_temp > /dev/null;
	# exclude duplicate variants
	awk '{ print $2,$4 }' ${output_stem}_chr${CHR}_temp.bim | uniq -f1 -D | awk '{print $1}' >> ${output_stem}_chr${CHR}.dups     
	plink2 --bfile ${output_stem}_chr${i}_temp --exclude ${output_stem}_chr${CHR}.dups --memory 15000 --make-bed --out ${output_stem}_chr${CHR}_temp_nodups > /dev/null ;
	# update name with rs numbers
	plink2 --bfile ${output_stem}_chr${CHR}_temp_nodups --update-name ${snp_names_file}_chr${CHR}.txt --memory 15000 --make-bed --out ${output_stem}_chr${CHR}_temp_nodups_names > /dev/null; 

	#print out numbers of variants
	total_var=$(( $(grep "out of" ${output_stem}_chr${CHR}_temp.log | awk 'BEGIN { ORS="+" } { print $4 }' | sed 's/\(.*\)+/\1 /' ) ))
	afterR2=$( cat ${output_stem}_chr${CHR}_temp.bim | wc -l )
	nomulti=$( cat ${output_stem}_chr${CHR}_temp_nodups.bim | wc -l )
	printf "$total_var variants after imputation, $afterR2 variants with R2>0.8, and $nomulti variants with unique positions for chr ${CHR}\n\n"
	
	grep 'pass filters and QC' ${output}.log
else
    output=$output_stem
    printf "Skipping the conversion to plink format and filtering of the individual chromosome because the -x flag was supplied! Picking up at updating sample IDs and sex in the merged file. Assuming the file stem is $output\n"
fi


#update SNP names, sex, and perform standard SNP filtering
printf "Step 4 : Updating sex in the fam file and applying standard variant filters\n\n"
#update person ids using the race and sex file
output_last=$output
output=${output}_IDs
awk '{ print "0 "$1"_"$2" "$1" "$2 }' $sex_file > ${output_last}_update_ids.txt
plink --bfile $output_last --update-ids ${output_last}_update_ids.txt --memory 15000 --make-bed --out $output > /dev/null

#update SNP names and add sex back into fam file
output_last=$output
output=${output}_sex
plink --bfile ${output_last} --update-sex ${sex_file} --memory 15000 --make-bed --out ${output} > /dev/null

tmp=$( grep "people updated" ${output}.log | awk '{ print $2 }' )
printf "$tmp people whose ids and sex were able to be updated using the sex file.\n"

## it is unclear why this section is here. If sex isn't present in the fam file (or in the supplied file with sex info), there are other issues. 
#n_sex=$( awk '{print $5}' ${output}.fam | sort | uniq -dc | wc -l)
#if [ $n_sex -ne 2 ];
#then
#    printf "\nERROR: sex not proper/updated in fam:
#${output}.fam , 1=Male, 2=Female;
#$(awk '{print $5}' ${output}.fam | sort | uniq -dc )\n"
#    exit 1;
#fi

###### merge back in genotypes #####

printf "\nStep 5: Merging back in the original genotypes.\n"
printf " reference aligned genotype file, ending in *-updated : ${preimputation_geno} \n"

# the genotyped variant ids in TOPMED imputation server format chrX:POS:REF:ALT
# this step was performed in part2, however, keeping it for robustness
printf "Updating genotype var_ID to TOPMED imputation server format chrX:POS:REF:ALT...\n"
awk '{ print "chrX:"$4":"$6":"$5,$2}' ${preimputation_geno}.bim > ${output_folder}/${preimputation_geno##*/}_TOPMED_varID.txt
plink --bfile ${preimputation_geno} --update-name ${output_folder}/${preimputation_geno##*/}_TOPMED_varID.txt 1 2 --make-bed --memory 15000 --out ${output_folder}/${preimputation_geno##*/}_TOPMED_varID  > /dev/null
# update to RSIDs
printf "Now updating genotype var_ID updated to RSID...\n"
plink --bfile ${output_folder}/${preimputation_geno##*/}_TOPMED_varID --update-name ${snp_names_file}_chr${CHR}.txt --make-bed --memory 15000 --out ${output_folder}/${preimputation_geno##*/}_rsid > /dev/null
 
# the genotyped variant ids to exclude
awk '{ print $2 }' ${output_folder}/${preimputation_geno##*/}_rsid.bim > ${output_folder}/${preimputation_geno##*/}_genotyped_variants.txt

# remove variants for which there are genotypes from the bim file
output_last=$output
output=${output}_nogeno
plink --bfile ${output_last} --exclude ${output_folder}/${preimputation_geno##*/}_genotyped_variants.txt --make-bed --memory 15000 \
--out $output > /dev/null
printf "$(grep -e "variants remaining" ${output}.log | awk '{ print $2 }') variants remaining after removing genotyped variants from imputation results.\n"

# merge the genotyped and imputed data
plink --bfile ${output} --bmerge ${output_folder}/${preimputation_geno##*/}_rsid --make-bed --memory 15000 --out ${output}_merged > /dev/null

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
    printf "Getting same position warnings. \n"
    grep "Warning: Variants" ${output}_merged.log | awk  '{ print $3"\n"$5 }' | sed -e "s/'//g" >  ${output_folder}/genotyped_variants_sameposwarnings.txt
    plink --bfile ${output} --exclude ${output_folder}/genotyped_variants_sameposwarnings.txt --make-bed --out ${output}2 > /dev/null
    printf "Confirm variants in ${output_folder}/genotyped_variants_sameposwarnings.txt for same position \n"
    printf "Removing $(wc -l < ${output_folder}/genotyped_variants_sameposwarnings.txt ) variants from the imputed dataset and re-attempting merge.\n"
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
printf "Output file: ${output} \n"

# SNP filters
printf "Step 6 : Applying Hardy-Weinberg equilibrium and MAF filters.\n" 
# Only run the HWE filter if there are females present
if [[ $( awk '{print $5}' ${output}.fam | grep "2" ) ]];
then
    output_last=${output}
    output_females=${output}_1e6femHWE

    printf "Running the HWE filter based on p values in females only...\n"

    # run HWE filter using only females
    plink --bfile ${output_last} --filter-females --hwe 0.000001 --memory 15000 --make-bed --out ${output_females}  > /dev/null

    # extract list of variants to remove from all samples
    awk '{print $2}' ${output_females}.bim > ${output_females}_SNPs.txt

    # output numbers to the log file
    grep -e 'variants loaded from .bim file' ${output_females}.log
    grep -e "hwe: " ${output_females}.log

    # filter those variants from all samples
    output=${output}_femHWE6
    plink --bfile ${output_last} --extract ${output_females}_SNPs.txt --memory 15000 --make-bed --out ${output} > /dev/null
    grep -e 'pass filters and QC' ${output}.log
else
    printf "\nSkipping Hardy Weinberg equilibrium filter because no females are present.\n"
fi

printf "Running MAF filter...\n"
output_last=$output
output=${output}_maf01
plink --bfile ${output_last} --maf 0.01 --memory 15000 --make-bed --out ${output} > /dev/null
grep -e "removed due to minor allele threshold" -e 'pass filters and QC' ${output}.log


# printf "\n\n ##### heterozygosity check & autosomal PC outlier removal ##### \n\n"
# printf " Skipped (PC outlier removal) as these steps are already completed in Part1 \n" 
# #if [ 'false' ]
# #then
# ###### heterozygosity check #####

# if [  $n_sex -eq 1 ] && [ ${sex} -eq 2 ];
# then
#         printf "\nSkipping Step6: #### Heterozygocity check (males only) #### \n"
#         printf "\n Only females or other sex present \n"
# else

# 	printf "\n: #### Heterozygocity check (males only) #### \n"
# 	printf " PS: plink requires autosomal chromosome for --het,  using dog as species to treat X-chromosome as autosomal for this step.
#          and circumvent 'Error: --het requires at least one polymorphic autosomal marker.' \n"
# 	## filter males and prune set (removed --filter-males flag)
# 	plinkset_x_out=${output};
# 	#plink --bfile ${plinkset_x_out} --indep-pairwise 200 100 0.2 --allow-no-sex --out ${plinkset_x_out}_prune --memory 15000 > /dev/null
# 	plink --bfile ${plinkset_x_out} --filter-males --indep-pairwise 200 100 0.2 --allow-no-sex --out ${plinkset_x_out}_prune --memory 15000 > /dev/null
# 	plink --bfile ${plinkset_x_out} --output-missing-phenotype 1 --extract ${plinkset_x_out}_prune.prune.in --make-bed --out ${plinkset_x_out}_pruned > /dev/null
# 	rm ${plinkset_x_out}_prune.*
# 	printf "$( wc -l < ${plinkset_x_out}_pruned.bim ) variants(not males only) out of $( wc -l < ${plinkset_x_out}.bim ) left after pruning.\n"

# 	# using "dog" as species to treat X as autosomal and overcome "Error: --het requires at least one polymorphic autosomal marker."
# 	plink --bfile ${plinkset_x_out}_pruned --het --dog --out ${plinkset_x_out}_pruned_hetcheck > /dev/null
# 	##plot and check for outliers
# 	Rscript plot_het_check_outliers.R ${plinkset_x_out}_pruned_hetcheck
# 	printf "Heterozygocity outliers(X-chr, males only) review figure: ${plinkset_x_out}_pruned_hetcheck.png \n"

# 	#if there are outliers >6 sd from the F stat mean, remove them
# 	if [ -f "${plinkset_x_out}_pruned_hetcheck_outliers.txt" ];
# 	then
# 	    plinkset_x_in=${plinkset_x_out}
# 	    plinkset_x_out=${plinkset_x_out}_nohetout
# 	    plink --bfile ${plinkset_x_in} --remove ${plinkset_x_in}_pruned_hetcheck_outliers.txt --make-bed --out ${plinkset_x_out} --memory 15000 > /dev/null
# 	    printf " Input: ${plinkset_x_in} \n"
# 	    grep -e ' people pass filters and QC' ${plinkset_x_out}.log
# 	    printf "Output: ${plinkset_x_out} \n"
# 	fi

# fi # Skip heterogygocity outlier check
# grep -e ' people pass filters and QC' ${output}.log
# printf "Output file: $output \n"

# ##### PC calculation ####
# #printf "\nStep 7: Calculating post-imputation PCs\n\n"

##### Subset to Overlapping SNPs between sexes #####
# skip overlapping SNP step for single sex dataset
n_sex=$(awk '{print $5}' ${output}.fam | sort | uniq -dc | wc -l)
if [ $n_sex -ne 2 ];
then
        printf "Warning: single sex dataset skipping subset to overlapping SNPs between sexes.\n"
else
	printf "\nStep 8: Subset to overlapping SNPs between sexes \n\n"
	# filter to each sex and use a missingness flag to filter out variants not present for those samples, writing out just a list of variants. 
	plink --bfile ${output} --filter-females --geno 0.01 --write-snplist --out ${output}_females --memory 15000 > /dev/null
	plink --bfile ${output} --filter-males --geno 0.01 --write-snplist --out ${output}_males --memory 15000 > /dev/null
	# sort those files and get overlapping variant IDs
	cat ${output}_males.snplist | sort > ${output}_males_sorted.snplist
	cat ${output}_females.snplist | sort > ${output}_females_sorted.snplist
	comm -12 ${output}_males.snplist ${output}_females.snplist  > ${output}_overlapping_SNPs.txt

	output_last=$output
	output=${output}_overlapping
	# now extract just those variants from the file with all the samples
	plink --bfile ${output_last} --extract ${output_last}_overlapping_SNPs.txt --real-ref-alleles --make-bed --out ${output} --memory 15000 > /dev/null
	printf "Output: ${output} \n"
	printf " female SNPs: $(wc -l < ${output_last}_female_snps.txt) \n"
	printf " male SNPs: $(wc -l < ${output_last}_male_snps.txt) \n"
	printf " Overlapping SNPs in both sexes: $(wc -l < ${output_last}_overlapping_SNPs.txt) \n"
	grep -e ' people pass filters and QC' ${output}.log
fi

# check no hh warnings on females and set hh to missing
printf "\n Confirm no hh warning for females only below: \n"
# will have issues with males only cohort. skip this for male only
plink --bfile ${output} --filter-females --freq --memory 15000 --out $output > /dev/null

# set hh missing
output_last=$output
output=${output}_nohh
printf "\n Set het. haploid genotypes missing in ${output}\n"
plink --bfile ${output_last} --set-hh-missing --memory 15000 --make-bed --out ${output} > /dev/null

### This is the step 

# STEP 9: Remove PC outliers based on final plinkset
##### Step 9: Remove individuals based on autosomal PC outliers) #####
printf "\nStep 9: #### Remove individuals based on autosomal PCs outlier, restricting to final plinkset provided  #### \n"

awk '{print $1,$2}' ${ds_plinkset}.fam > ${output}_ind_to_keep.txt
printf "\n individuals in final autosomal dataset(individual list: ${output}_ind_to_keep.txt): $(wc -l < ${output}_ind_to_keep.txt) \n"

plinkset_x_in=${output}
plinkset_x_out=${plinkset_x_in}_noout
plink --bfile ${plinkset_x_in} --keep ${plinkset_x_in}_ind_to_keep.txt --make-bed --out ${plinkset_x_out} --memory 15000 > /dev/null
    printf " Input: ${plinkset_x_in} \n"
    printf " Output: ${plinkset_x_out} \n"

#make final file
output_last=$plinkset_x_out
output=${output_stem}
plink --bfile ${output_last} --memory 15000 --make-bed --out ${output_stem} > /dev/null
printf "\n Final X-chromosome plinkset: ${output_stem} \n Part Postimputation ... Done! \n"

