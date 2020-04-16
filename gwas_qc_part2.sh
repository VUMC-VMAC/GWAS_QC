#!/bin/bash

stem=$1
input_fileset=$2
ref_file_stem=$3

#fail if error
set -e

#define usage
display_usage() {
    printf "GWAS QC Part 2 script
Picks up after decisions have been made regarding PC outlier removal. Input dataset can be either the non-Hispanic white subset or include all races present. This script will perform the Hardy-Weinberg filter and prepare the files for upload to the TopMed Imputation Server. Will create all files within the same folder as the current plink fileset.

Usage:
SCRIPTNAME.sh [stem] [input_fileset] [ref_file_stem]

stem = the prefix you want for the files right before imputation
input_fileset = the full path and file stem for the current plink fileset '*[bed,bim,fam]'
ref_file_stem = the full file path and stem for the reference panel. Assumes it is split by chromosome and named STEM_chr*.txt.gz
"
        }
if [  $# -lt 2 ]
        then
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
Run_PCcalc_infolder.sh
plot_PCs_generate_ids_to_keep.R
check_ambig_snps.R
HRC-1000G-check-bim-NoReadKey.pl\n"
        exit 1
fi

#print out inputs
echo -e "GWAS QC Part 2 Script\n"
echo "Input fileset : "$input_fileset
echo "Reference panel SNP file : "$ref_file
echo "Stem for pre-imputation files :"$stem

#get file path for input
input_path=${input_fileset%/*}
#add input path to the stem for the final files
stem=${input_path}/${stem}

#### HWE filter ####

echo -e "\nRunning the Hardy-Weinberg Equilibrium SNP filter \n"  
output=$( echo ${input_fileset}_hwe6 )
plink --bfile $input_fileset --hwe 0.000001 --make-bed --out $output > /dev/null
grep -e "--hwe:" -e ' people pass filters and QC' $output.log   
echo -e " outfile: $output \n"

#### liftOver ####

echo -e "Lifting over the genotype file from build 37 to build 38 to match the TopMed reference panel"
#create bed file to start the liftover
awk '{ print "chr"$1" "$4 -1" "$4" "$2 }' ${output}.bim > ${output}_toliftover.txt
#lift over
./liftOver ${output}_toliftover.txt hg19ToHg38.over.chain.gz ${output}_lifted.txt ${output}_unlifted.txt
#remove "chr" from chromosome column
sed -i 's/^chr//g' ${output}_lifted.txt

#get list of variants to remove (those with non-standard allele codes or which cannot be lifted)
grep -e "random" -e "alt" -e "fix" ${output}_lifted.txt | awk '{ print $4 }' > ${output}_lift_issue_snps_todrop.txt
awk '{ print $4 }' ${output}_unlifted.txt | sort -u >> ${output}_lift_issue_snps_todrop.txt
echo -e "Removing $( wc -l < ${output}_lift_issue_snps_todrop.txt ) variants which fail liftOver ($( grep "random" ${output}_lifted.txt | wc -l ) 'random', $( grep "alt" ${output}_lifted.txt | wc -l ) 'alt', $( grep "fix" ${output}_lifted.txt | wc -l ) 'fix', and $( awk '{ print $4 }' ${output}_unlifted.txt | sort -u | wc -l ) not present in b38)."

#drop 
plink --bfile $output --exclude ${output}_lift_issue_snps_todrop.txt --make-bed --out ${output}_noprobSNPs > /dev/null
#update the chr and BP positions for remaining
plink --bfile ${output}_noprobSNPs --update-chr ${output}_lifted.txt 1 4 --make-bed --out ${output}_noprobSNPs_chr > /dev/null
plink --bfile ${output}_noprobSNPs_chr --update-map ${output}_lifted.txt 3 4 --make-bed --out ${output}_noprobSNPs_chr_bplifted > /dev/null
output=${output}_noprobSNPs_chr_bplifted

#### Imputation prep ####

echo -e "Comparing with the reference panel and preparing files for imputation for each chromosome"

#run imputation prep for each chromosome
for i in {1..22} ;
do
    echo -e "Beginning chr ${i}..."
    #split out this chr
    plink --bfile ${output} --allow-no-sex --chr $i --freq --make-bed --out ${stem}_chr$i > /dev/null

    #run the imputation checking script
    perl HRC-1000G-check-bim-NoReadKey.pl -b ${stem}_chr$i.bim  -f ${stem}_chr$i.frq -r ${ref_file_stem}_chr${i}.txt.gz -h -n -c > /dev/null

    #run created script with all the plink commands; breaks into files for each chr
    sed -i "s|rm TEMP|rm ${input_path}/TEMP|" ${input_path}/Run-plink.sh
    sh ${input_path}/Run-plink.sh > /dev/null

    #update chr code to have chr# rather than just #
    bcftools annotate --rename-chrs update_chr_names_b38.txt ${stem}_chr${i}-updated-chr${i}.vcf -o ${stem}_chr${i}-updated-chr${i}_chrupdate.vcf -O v 
    #create a sorted *.vcf.gz file using VCFtools and tabix
    vcf-sort ${stem}_chr${i}-updated-chr${i}_chrupdate.vcf | bgzip -c > ${stem}_chr${i}-updated-chr${i}_chrupdate.vcf.gz;

    echo -e "Chr ${i} complete..."
done

echo -e "Removed $( cat Exclude-* | wc -l ) variants for mismatch with the reference panel leaving $( cat ${stem}_chr${i}-updated-chr${i}.bim | wc -l ) for imputation.\n"

#remove the intermediate .vcf files
rm ${input_path}/*.vcf

echo -e "\nConversion complete! Upload the files (${output}-updated-chr*_chrupdate.vcf.gz) to the imputation server.\n"

