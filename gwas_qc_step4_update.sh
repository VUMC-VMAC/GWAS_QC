#!/bin/bash
# step 4 sex check update: minor updates/more verbose on others

##### initial SNP filters: missingness (--geno), minor allele frequency (--maf) ##### VJ 
echo -e "\n\nStep 1: Remove SNPs with >5% missingness or with MAF <0.01\n"
output=$( echo ${output_dir}/${stem}_geno05_maf01 ) # VJ: output_dir should not have '/' in end 
plink --bfile $input_fileset --geno 0.05 --maf 0.01 --make-bed --out $output > /dev/null


grep 'loaded from' $output.log  | sed  's/loaded from.*//' | head -n2
grep -e 'missing genotype data' -e 'minor allele threshold' -e ' people pass filters and QC' $output.log
echo -e "Output file: $output \n"

##### person missingness (--mind) #####
echo -e "\nStep 2: remove subjects w/ >1% missingness\n"
output_last=$output
output=${output}_mind01
plink --bfile $output_last --mind 0.01 --make-bed --out $output > /dev/null


grep -e 'removed due to missing genotype data' -e ' people pass filters and QC' $output.log
echo -e "Output file: $output \n"

##### Relatedness (--genome full unbounded nudge --min 0.20, --remove iid fid text file) #####
# 'full' adds raw pairwise comparison data to the report.
# minimum 0.20 PI_HAT specified for inclusion in --genome report.

echo -e "\nStep 3: Remove related individuals\n"
echo -e "   Identify relateds"
plink --bfile $output --genome full unbounded nudge --min 0.20 --out ${output}_relatedness > /dev/null
echo -e "\n output files: \n$(ls ${output}_relatedness.*)" #VJ

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
echo -e "Rscript: Rscript get_related_ids.R \n input: ${output}_relatedness (stem) \n output: ${output}_relatedness_related_ids.txt" #VJ
Rscript get_related_ids.R ${output}_relatedness


# remove selected ids from the last generated genotype file set
echo -e "   Remove relateds"
echo -e "Removing $( wc -l ${output_last}_relatedness_related_ids.txt ) individuals for relatedness."
output_last=$output
output=${output}_norelated
plink --bfile $output_last --remove ${output_last}_relatedness_related_ids.txt --make-bed --out $output > /dev/null
echo -e "\n output files: \n$(ls ${output}.*)" #VJ
grep ' people pass filters and QC' $output.log
echo -e "Output file: $output \n"



# STEP 4 update

##### sex check #####
# for pheno sex and fam sex mismatch: -- remove ${output}_pheno_fam_sex_mismatch.txt
# for missing sex: --update-sex $race_sex_file 2 (n=2, n+2 th column) 
# 
echo -e "\nStep 4: sex check\n"
echo -e " using sex from file: $race_sex_file" # VJ

# VJ: perhaps its better to wait for updating sex till we have all three sex sources: geno, pheno, and .fam file. 
# VJ check for descripancy in pheno 'race_sex_file' and  fam '*_norelated' sex
# if condition checks for pheno sex provided 
if [ "$( awk '{ print $4 }' $race_sex_file | sort -u )" != 0 ];
then 
  # get mismatched in pheno and fam sex
  join -j1 --check-order <(<${output}.fam awk '{print $1"_"$2, $3, $4, $5, $6}' OFS='\t' | sort -k1,1) <(<$race_sex_file awk '{print $1"_"$2, $3, $4}' OFS='\t' | sort -k1,1) | awk '$4!=$7{print}' | awk '{print $1}' | sed 's/_/\t/g' > ${output}_pheno_fam_sex_mismatch.txt
   
  # VJ: Is there mismatch?
  if [ $( wc -l < ${output}_pheno_fam_sex_mismatch.txt ) -gt 0 ];
  then

    # do sex-check; display pheno, fam, and geno sex; remove the phenotype mismatches
    output_checking_sex=${output}_checking_sex
    plink --bfile $output --check-sex --out ${output}_checking_sex > /dev/null
    echo -e "\n output files: \n$(ls ${output}_checking_sex.*)"  #VJ
    grep -e 'check-sex: ' ${output}_checking_sex.log 

    # join .fam, pheno, sexcheck files  and display
    #echo -e "$(head -n1 ${output}_checking_sex.sexcheck)\t pat\t mat\t fam_sex\t pheno\t eth\t pheno_sex" 
    echo -e "IID\t FID\t fam_sex\t SNP_sex\t \t Ethnicity\t pheno_sex"
    join -j1 --check-order <(<${output}.fam awk '{print $1"_"$2, $3, $4, $5, $6}' OFS='\t' | sort -k1,1) <(<$race_sex_file awk '{print $1"_"$2, $3, $4}' OFS='\t' | sort -k1,1)| join -j1 --check-order <(awk '{print $1"_"$2, $3, $4, $5, $6}' ${output}_checking_sex.sexcheck  OFS='\t' | sort -k1,1) - | sed 's/_/\t/' | awk '$3!=$12{print $1, $2, $3, $4, $11, $12}' OFS='\t'

    # join three fam, pheno, geno (only Problems)
join -j1 --check-order <(<${output}.fam awk '{print $1"_"$2, $3, $4, $5, $6}' OFS='\t' | sort -k1,1) <(<$race_sex_file awk '{print $1"_"$2, $3, $4}' OFS='\t' | sort -k1,1)| join -j1 --check-order <(awk '$5=="PROBLEM"{print $1"_"$2, $3, $4, $5, $6}' ${output}_checking_sex.sexcheck | sort -k1,1) - | sed 's/_/\t/' 

    ## remove phenotype-mismatches: (--remove) ${output}_pheno_fam_sex_mismatch.txt (here or remove with sex_check)
    echo -e " Removing phenotype sex mismatches"
    echo -e "Removing $( wc -l <${output}_pheno_fam_sex_mismatch.txt ) individuals for phenotype sex mismatch."
    output_last=$output
    output=${output}_pheno_mismatch
    plink --bfile $output_last --remove ${output_last}_pheno_fam_sex_mismatch.txt --make-bed --out $output > /dev/null
    grep ' people pass filters and QC' $output.log
    echo -e "Output file: $output \n"
  fi
fi


#only update sex if there is something more than missing values "0" for sex in the provided
 file, provides flexibility to updateonly selected individuals with non-missing sex #VJ
if [ "$( awk '{ print $4 }' $race_sex_file | sort -u )" != 0 ];
then 
    # VJ: output/output_last are only updated if pheno-sex-mismatch is found in previous step
    output_last=$output
    output=${output}_sex
    plink --bfile $output_last --update-sex <(awk '$4>0' $race_sex_file) 2 --make-bed --out $output > /dev/null #VJ
    echo -e "\n output files: \n$(ls ${output}.*)" # VJ
    echo -e "Updated sex from the provided file for non-missing : $race_sex_file"
else
    echo -e "No non-missing sex information in the provided file (${race_sex_file}). Using sex from fam file."
fi
 
# VJ display which iid's are updated (as now part update is an option)
echo -e "$output \n $output_last"
echo -e "IID \t FID \t sex_new \t sex_last"
join -j1 --check-order <(<${output}.fam awk '{print $1"_"$2, $3, $4, $5, $6}' | sort -k1,1) <(<$output_last.fam awk '{print $1"_"$2, $3, $4, $5, $6}' | sort -k1,1) | awk '$4!=$8{print $1, $4, $8}' OFS='\t' | sed 's/_/\t /g' | tee >( echo -e " Number of iid's updated: $(wc -l ) " )
 
 
#do the check (did the check before so can remove this, just need the file name)
plink --bfile $output --check-sex --out ${output}_checking_sex > /dev/null
echo -e "\n output files: \n$(ls ${output}_checking_sex.*)"  #VJ
grep -e 'check-sex: ' ${output}_checking_sex.log 

# print problems
head -n1 ${output}_checking_sex.sexcheck # VJ
awk '$5=="PROBLEM"{print}' ${output}_checking_sex.sexcheck | sort -k1,1 #VJ

# VJ 20200421: phenotype mismatch(output, output_last) check (joining the fam, pheno and sex_checkfile and displaying the options)


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


