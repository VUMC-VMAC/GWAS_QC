#!/bin/bash

#fail if error
set -e

# VJ: singularity image and directories bound for reproducibility
[ ! -z "$SINGULARITY_CONTAINER" ] && printf "\nSingularity image: $SINGULARITY_CONTAINER"
[ ! -z "$SINGULARITY_COMMAND" ] && printf "\nSingularity shell: $SINGULARITY_COMMAND"
[ ! -z "$SINGULARITY_BIND" ] && printf "\nMapped directories:\n $SINGULARITY_BIND\n" | sed 's/,/\n /g'

#define usage
display_usage() {
    printf "GWAS QC Part 2 script

Picks up after decisions have been made regarding PC outlier removal. Input dataset can be either the non-Hispanic white subset or include all races present. This script will perform the Hardy-Weinberg filter and prepare the files for upload to the TopMed Imputation Server. Will create all files within the same folder as the current plink fileset.

Usage:
SCRIPTNAME.sh -o [output_stem] -i [input_fileset] -R [ref_file] -b [input genome build] -n -c


output_stem = the prefix you want for the files right before imputation

input_fileset = the full path and file stem for the current plink fileset '*[bed,bim,fam]'

ref_file = the full file path and name for the reference panel

-b = optional argument indicating the build of the input dataset; default assumption is b37 but can supply b36 or b38. This will impact the liftOver step.

-n = optional argument set to indicate not to exclude variants for not being in or not matching the reference panel; default is to exclude

-c = optional argument to skip clean-up

-h will show this usage
"
}
#set default to do the exclusion and build to b37
noexclude='false'
build='b37'
skip_cleanup='false'
#parse arguments
while getopts 'o:i:R:b:nch' flag; do
  case "${flag}" in
    o) output_stem="${OPTARG}" ;;
    i) input_fileset="${OPTARG}" ;;
    R) ref_file="${OPTARG}" ;;
    b) build="${OPTARG}" ;;
    n) noexclude='true' ;;
    c) skip_cleanup='true' ;;
    h) display_usage ; exit ;;
    \?|*) display_usage
       exit 1;;
  esac
done

#validate arguments
if [ -z "$output_stem" ] || [ -z "$input_fileset" ] || [ -z "$ref_file" ] ;
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
printf "GWAS QC Part 2 Script

Input fileset : $input_fileset 
Reference panel SNP file : $ref_file
Stem for pre-imputation files :$output_stem
"

#get file path for input
input_path=${input_fileset%/*}
output_path=${output_stem%/*}

#get input file name to be the stem for the output
input_stem=${input_fileset##*/}

#### remove palindromic ####
output=$( printf ${output_path}/${input_stem} )
printf "\nStep 2: Removing palindromic variants \n"

#get palindromic variants
awk '{ if( ($5 == "A" && $6 == "T") || ($5 == "T" && $6 == "A") || ($5 == "C"  && $6 == "G") || ($5 == "G" && $6 == "C")) print }'  ${output}.bim > ${output_path}/palindromic_snps.txt

#remove them
output_last=$output
output=${output}_nopal
plink --bfile $output_last --exclude ${output_path}/palindromic_snps.txt --output-chr M --make-bed --out $output > /dev/null
printf "Removed $( wc -l < ${output_path}/palindromic_snps.txt ) palindromic variants.\n"
grep ' people pass filters and QC' ${output}.log
printf "Output file: $output \n"

#### liftOver ####

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
	printf "Warning: single sex dataset skipping Step: Subset to Overlapping SNPs between sexes \n"
	printf "N sexes: $n_sex, sex code: $sex \n"
else
	plink --bfile ${output} --filter-females --make-bed --out ${output}_females --memory 15000 > /dev/null
	plink --bfile ${output} --filter-males --make-bed --out ${output}_males --memory 15000 > /dev/null
	awk '{print $2}' ${output}_females.bim | sort > ${output}_female_snps.txt
	awk '{print $2}' ${output}_males.bim | sort > ${output}_male_snps.txt
	comm -12 ${output}_male_snps.txt ${output}_female_snps.txt  > ${output}_overlapping_SNPs.txt
fi

output_last=$output
output=${output}_overlap
if [ $n_sex -ne 2 ];
then
	printf "Creating ${output} plinkset for consistency \n"
	plink --bfile ${output_last} --real-ref-alleles --make-bed --chr 23 --out ${output} --memory 15000 > /dev/null
	printf "Output: ${output} \n"
	grep -e ' people pass filters and QC' ${output}.log
else
	# limit ot chr 23, Keep the allele order defined in the .bim file,instead of forcing A2 to be the major allele.--real-ref-alleles also removes 'PR' from the INFO values emitted by --recode vcf{-fid/-iid}.
	plink --bfile ${output_last} --extract ${output_last}_overlapping_SNPs.txt --real-ref-alleles --make-bed --chr 23 --out ${output} --memory 15000 > /dev/null
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

plink --bfile ${output_last} --update-name ${output_last}_TOPMED_varID.txt 1 2 --make-bed --memory 15000 \
--out ${output}  > /dev/null

# create VCF
plink --bfile ${output} --real-ref-alleles --recode-vcf --chr 23 --out ${output} --memory 15000 > /dev/null

printf "Output: ${output} \n"
grep -e ' people pass filters and QC' ${output}.log


#run imputation prep for each chromosome
printf "Beginning Chr X...\n"

#update chr code to have chr# rather than just #, sort and output vcf
mkdir ${output_path}/tmpX/
bcftools annotate --rename-chrs /scripts_local/update_chr_name.txt ${output}.vcf | bcftools sort --temp-dir  ${output_path}/tmpX/ -O z -o ${output_stem}_chrX.vcf.gz > /dev/null
printf "\nChr X complete...\n"

#print out number of variants actually excluded or which would have been excluded
if [ "$noexclude" = true ];
then 
    printf "Would have removed $( cat ${output_path}/Exclude-* | wc -l ) variants for mismatch with the reference panel, being palindromic with MAF > 0.4, or being absent from the reference panel leaving $( cat ${output_stem}-updated-chr23_overlap.bim | wc -l ) for imputation, but the no-exclude option was specified.\n"
else
    printf "Removed $( cat ${output_path}/Exclude-* | wc -l ) variants for mismatch with the reference panel, being palindromic with MAF > 0.4, or being absent from the reference panel leaving $( cat ${output_stem}-updated-chr23_overlap.bim | wc -l ) for imputation.\n"
fi

#remove the intermediate .vcf and .bed files
##make sure that it only removes the files from running this script, not the initial files or files for other datsets
if [ "$skip_cleanup" = false ];
then
    rm ${output_path}/*.vcf
    files_to_remove=$( find ${output_path}/${input_stem}*.bed | grep -v "${output}.bed" | grep -v "${output_stem}-updated.bed" )
    rm $files_to_remove
else
    printf "Skipping cleanup. Please manually remove unnecessary files.\n"
fi

printf "\nConversion complete! Upload the file (${output_stem}_chrX.vcf.gz) to the imputation server.\n"