#!/bin/bash

#fail if error
set -e

#define usage
display_usage() {
    printf "GWAS QC Part 2 script

Picks up after decisions have been made regarding PC outlier removal. Input dataset can be either the non-Hispanic white subset or include all races present. This script will perform the Hardy-Weinberg filter and prepare the files for upload to the TopMed Imputation Server. Will create all files within the same folder as the current plink fileset.

Usage:
SCRIPTNAME.sh -o [output_stem] -i [input_fileset] -R [ref_file_stem] -n

output_stem = the prefix you want for the files right before imputation
input_fileset = the full path and file stem for the current plink fileset '*[bed,bim,fam]'
ref_file_stem = the full file path and stem for the reference panel. Assumes it is split by chromosome and named STEM_chr*.txt.gz
-n = optional argument set to indicate not to exclude variants for not being in or not matching the reference panel; default is to exclude

-h will show this usage
"
}
#parse arguments
while getopts 'o:i:R:eh' flag; do
  case "${flag}" in
    o) output_stem='${OPTARG}' ;;
    i) input_fileset='${OPTARG}' ;;
    R) ref_file_stem="${OPTARG}" ;;
    n) noexclude='true' ;;
    h) display_usage ; exit ;;
    \?|*) display_usage
       exit 1;;
  esac
done

#validate arguments
if [ -z "$output_stem" ] || [ -z "$input_fileset" ] || [ -z "$ref_file_stem" ] ;
then
    printf "Error: Necessary arguments not present!\n\n"
    display_usage
    exit 1
fi


#check to make sure this is being run in the scripts folder (checking if necessary script is present)
if test ! -f get_related_ids.R ;
then
        printf "Currently in $PWD, but must be in a folder with the necessary scripts to run the GWAS QC! Please move to that folder and run this script again.
Necessary scripts:
get_related_ids.R
check_id_length.R
plot_het_check_outliers.R
plot_PCs_generate_ids_to_keep.R
check_ambig_snps.R
HRC-1000G-check-bim-NoReadKey.pl\n"
        exit 1
fi

#print out inputs
printf "GWAS QC Part 2 Script\n"
printf "Input fileset : "$input_fileset
printf "Reference panel SNP file : "$ref_file
printf "Stem for pre-imputation files :"$output_stem

#get file path for input
input_path=${input_fileset%/*}

#### HWE filter ####

printf "\nRunning the Hardy-Weinberg Equilibrium SNP filter \n"  
output=$( printf ${input_fileset}_hwe6 )
plink --bfile $input_fileset --hwe 0.000001 --make-bed --out $output > /dev/null
grep -e "--hwe:" -e ' people pass filters and QC' $output.log   
printf " outfile: $output \n"

#### liftOver ####

printf "Lifting over the genotype file from build 37 to build 38 to match the TopMed reference panel"
#create bed file to start the liftover
awk '{ print "chr"$1" "$4 -1" "$4" "$2 }' ${output}.bim > ${output}_toliftover.txt
#lift over
./liftOver ${output}_toliftover.txt hg19ToHg38.over.chain.gz ${output}_lifted.txt ${output}_unlifted.txt
#remove "chr" from chromosome column
sed -i 's/^chr//g' ${output}_lifted.txt

#get list of variants to remove (those with non-standard allele codes or which cannot be lifted)
grep -e "random" -e "alt" -e "fix" ${output}_lifted.txt | awk '{ print $4 }' > ${output}_lift_issue_snps_todrop.txt
awk '{ print $4 }' ${output}_unlifted.txt | sort -u >> ${output}_lift_issue_snps_todrop.txt
printf "Removing $( wc -l < ${output}_lift_issue_snps_todrop.txt ) variants which fail liftOver ($( grep "random" ${output}_lifted.txt | wc -l ) 'random', $( grep "alt" ${output}_lifted.txt | wc -l ) 'alt', $( grep "fix" ${output}_lifted.txt | wc -l ) 'fix', and $( awk '{ print $4 }' ${output}_unlifted.txt | sort -u | wc -l ) not present in b38)."

#drop 
plink --bfile $output --exclude ${output}_lift_issue_snps_todrop.txt --make-bed --out ${output}_noprobSNPs > /dev/null
#update the chr and BP positions for remaining
plink --bfile ${output}_noprobSNPs --update-chr ${output}_lifted.txt 1 4 --make-bed --out ${output}_noprobSNPs_chr > /dev/null
plink --bfile ${output}_noprobSNPs_chr --update-map ${output}_lifted.txt 3 4 --make-bed --out ${output}_noprobSNPs_chr_bplifted > /dev/null
output=${output}_noprobSNPs_chr_bplifted

#### Imputation prep ####

printf "Comparing with the reference panel and preparing files for imputation for each chromosome"

#run imputation prep for each chromosome
for i in {1..22} ;
do
    printf "Beginning chr ${i}..."
    #split out this chr
    plink --bfile ${output} --allow-no-sex --chr $i --freq --make-bed --out ${output_stem}_chr$i > /dev/null

    #run the imputation checking script
    perl HRC-1000G-check-bim-NoReadKey.pl -b ${output_stem}_chr$i.bim  -f ${output_stem}_chr$i.frq -r ${ref_file_stem}_chr${i}.txt.gz -h -n -c > /dev/null

    #if noexclude, then skip the exclusion step in the plink file
    if [ "$noexclude" == true ];
    then
	sed -i -e '1s/.*/#&/' -e "s|${input_path}/TEMP1|${stem}_chr${i}|g" ${input_path}/Run-plink.sh
    fi

    #run created script with all the plink commands; breaks into files for each chr
    sed -i "s|rm TEMP|rm ${input_path}/TEMP|" ${input_path}/Run-plink.sh
    sh ${input_path}/Run-plink.sh > /dev/null

    #update chr code to have chr# rather than just #
    bcftools annotate --rename-chrs update_chr_names_b38.txt ${output_stem}_chr${i}-updated-chr${i}.vcf -o ${output_stem}_chr${i}-updated-chr${i}_chrupdate.vcf -O v 
    #create a sorted *.vcf.gz file using VCFtools and tabix
    vcf-sort ${output_stem}_chr${i}-updated-chr${i}_chrupdate.vcf | bgzip -c > ${output_stem}_chr${i}-updated-chr${i}_chrupdate.vcf.gz;

    printf "Chr ${i} complete..."
done

#print out number of variants actually excluded or which would have been excluded
if [ "$noexclude" == true ];
then 
    printf "Would have removed $( cat ${input_path}/Exclude-* | wc -l ) variants for mismatch with the reference panel, being palindromic with MAF > 0.4, or being absent from the reference panel leaving $( cat ${output_stem}_chr*-updated-chr*.bim | wc -l ) for imputation, but the no-exclude option was specified.\n"
else
    printf "Removed $( cat ${input_path}/Exclude-* | wc -l ) variants for mismatch with the reference panel, being palindromic with MAF > 0.4, or being absent from the reference panel leaving $( cat ${output_stem}_chr*-updated-chr*.bim | wc -l ) for imputation.\n"
fi

#remove the intermediate .vcf files
rm ${input_path}/*.vcf

printf "\nConversion complete! Upload the files (${output_stem}-updated-chr*_chrupdate.vcf.gz) to the imputation server.\n"

