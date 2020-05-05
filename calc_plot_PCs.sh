#!/bin/sh

#fail on error
set -e

#define usage                                                                                                                                          
display_usage() {
    printf "PC calculation script

Use this script to calculate and plot PCs with or without 1000G samples.

Usage: SCRIPTNAME.sh -i [input_stem] -r [race_sex_file] -G [stem_1000G] -n


input_stem = the full path and file stem for the plink set for PC calculations '*[bed,bim,fam]'

race_sex_file (optional) = a file with FID and IID (corresponding to the fam file), 1 column indicating both race and ethnicity for PC plots, and another indicating sex for the sex check (1 for males, 2 for females, 0 if unknown), with NO header. Non-hispanic whites need to be indicated with 'White.' No other values in the race column must be fixed; however, the race column must not include spaces.

stem_1000G (optional) = the full path and stem to the 1000G genotype files in plink format. There must also be a file with FID, IID, race with the same stem and _race.txt as the suffix (ie for a plink file set like this: all_1000G.bed, all_1000G.bim, all_1000G.fam the race file would be like this all_1000G_race.txt)

-n (optional) = indicates to not create the exclusion file. Default is to create it

-h will show this usage
"
}


#parse arguments
while getopts 'i:r:G:nh' flag; do
  case "${flag}" in
    i) input_stem="${OPTARG}" ;;
    r) race_sex_file="${OPTARG}" ;;
    G) stem_1000G="${OPTARG}" ;;
    n) exclusion_file="no" ;;
    h) display_usage ; exit ;;
    \?|*) display_usage
       exit 1;;
  esac
done

#check to make sure necessary argument is present
if [ -z "$input_stem" ];
then
    printf "Error: Must specify genotype dataset to calculate PCs for!\n\n"
    display_usage ; exit 1
else
    printf "Input dataset for PC calculation : $input_stem\n"
fi
#notify user about other arguments
##race/eth
if [ -z "$race_sex_file" ];
then
    printf "No file with race/ethnicity was specified so plots will not color on those groups.\n"
    race_sex_file=NA
else
    printf "Race/ethnicity file : $race_sex_file\n"
fi
##1000G data
if [ ! -z "$stem_1000G" ];
then
    printf "1000G file stem : $stem_1000G
    Since 1000G files were specified, those samples will be included in PC calculation and plots.\n\n"
fi
##exclusion file to be created?
if [ "$exclusion_file" = "no" ];
then
    printf "No file for default exclusions will be created.\n"
else
    exclusion_file="yes"
fi


########################################### start of prep for PC calculation #######################################

#if 1000G data was supplied, then merge it with the input dataset
if [ ! -z $stem_1000G ];
then
    #### prep for PCs with 1000G data ###
    printf "Merging with 1000G for PC calculation to make ancestry filtering decisions\n"
    #get variants in current dataset
    awk '{ print $2 }' ${input_stem}.bim > ${input_stem}_snps.txt
    
    #subset 1000G dataset to only variants present here
    plink --bfile ${stem_1000G}  --output-missing-phenotype 1 --extract ${input_stem}_snps.txt --make-bed --out ${input_stem}_1000G_overlapping > /dev/null
    stem_1000G_formerge=${input_stem}_1000G_overlapping
    num_1000G_snps=$( wc -l < ${stem_1000G_formerge}.bim )
    printf "$num_1000G_snps overlap between 1000G and the current dataset\n\n"
    
    #attempt merge (keep the script from failing when this command fails)
    pcainput=${input_stem}_1000G_merged
    plink --bfile ${input_stem} --bmerge ${stem_1000G_formerge} --allow-no-sex --make-bed --out ${pcainput} || true &> /dev/null


    if [ -f "${pcainput}-merge.missnp" ];
    then

	#remove the mismatching variants
	plink --bfile ${input_stem} --exclude ${pcainput}-merge.missnp --allow-no-sex --make-bed --out ${input_stem}_formerge > /dev/null
	plink --bfile ${stem_1000G_formerge} --exclude ${pcainput}-merge.missnp --allow-no-sex --make-bed --out ${stem_1000G_formerge}_formerge  > /dev/null
	plink --bfile ${input_stem}_formerge --bmerge ${stem_1000G_formerge}_formerge --allow-no-sex --make-bed --out ${pcainput}  > /dev/null

	#filter to overlapping
	plink --bfile ${pcainput} --geno 0.05 --allow-no-sex --make-bed --out ${pcainput}_geno05  > /dev/null
	pcainput=${pcainput}_geno05
	printf "$( wc -l < ${pcainput}.bim ) variants in dataset merged with 1000G\n\n"
    fi
else
    pcainput=${input_stem}
fi

printf "Prune genotypes\n"

#prune and set phenotypes to 1
plink --bfile $pcainput --indep-pairwise 200 100 0.2 --allow-no-sex --out  ${pcainput}_forprune > /dev/null
plink --bfile $pcainput --output-missing-phenotype 1 --extract  ${pcainput}_forprune.prune.in --make-bed --out ${pcainput}_pruned > /dev/null
printf "$( wc -l < ${pcainput}_pruned.bim ) variants out of $( wc -l < ${pcainput}.bim ) left after pruning.\n\n"

#cleanup
rm ${pcainput}*forprune*

printf "Running smartpca...\n"
#create input file for smartpca
printf "genotypename: ${pcainput}_pruned.bed
snpname: ${pcainput}_pruned.bim
indivname: ${pcainput}_pruned.fam
evecoutname: ${pcainput}_pruned.pca.evec
evaloutname: ${pcainput}_pruned.eigenvalues
altnormstyle: NO
numoutevec: 10
numoutlieriter: 0
numoutlierevec: 10
outliersigmathresh: 6
qtmode: 0" > ${pcainput}_pccalc.par
#run smartpca
smartpca -p ${pcainput}_pccalc.par > ${pcainput}_pccalc.log

printf "smartpca complete! Plotting...\n\n"

#plot
if [ -z $stem_1000G ];
then
    #arguments: PCA file, race/sex file, optional file with 1000G race, yes/no to indicate whether to output a file for outlier exclusion
    Rscript plot_PCs_generate_ids_to_keep.R ${pcainput}_pruned.pca.evec $race_sex_file NA $exclusion_file
    printf "Plots of PCs are saved here: ${pcainput}_pruned_PCplots.pdf."
else  
    #arguments: PCA file, race/sex file, optional file with 1000G race, yes/no to indicate whether to output a file for outlier exclusion
    Rscript plot_PCs_generate_ids_to_keep.R ${pcainput}_pruned.pca.evec $race_sex_file ${stem_1000G}_race.txt $exclusion_file
    printf "Plots of PCs calculated with 1000G samples are saved here: ${pcainput}_pruned_PCplots.pdf."
fi

if [ "$exclusion_file" = "yes" ];
then
    printf "A file with ids for NHW who were not PC outliers (>5sd from the mean) is written out for your convenience if all outliers should be removed: ${pcainput}_pruned_nooutliers.txt\n"
fi

