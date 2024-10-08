output_stem=$1

output_folder=${output_stem%/*}

# get the initial post-imputation numbers
total_var=$(( $(grep "out of" ${output_stem}_chr*_temp.log | awk 'BEGIN { ORS="+" } { print $4 }' | sed 's/\(.*\)+/\1 /' ) ))
R2=$(( $total_var - $( cat ${output_stem}_chr*_temp.bim | wc -l ) ))
multi=$(( $total_var - $afterR2 - $( cat ${output_stem}_chr*_temp_nodups.bim | wc -l ) ))

# get the number of variants after these initial filters
samples=$( grep "pass filters and QC" ${output_stem}.log | awk '{ print $4 }' )
variants=$( grep "pass filters and QC" ${output_stem}.log | awk '{ print $1 }' )

# output these 
printf "$samples samples\n$total_var variants\n\nRemoved $R2 variants for R2<0.8, $multi duplicated or multi-allelic SNPs.\n\n$samples samples\n$variants variants\n\n"

#### this is where SNPWeights is run

# check for which sets were run in SNPWeights
output_stem=${output_stem}_IDs_sex

printf "Subset by ancestry using SNPWeights\n\n"

for i in EUR AFR AMR AMR2admx AMR3admx ; 
do 
  if [ -e "${output_stem}_${i}_geno01.bed" ]; 
  then 
    printf "For $i subset: \n\n"
    
    output=${output_stem}_${i}_geno01
    
    # first just get the number of samples dropped to get to this subset
    samples=$( grep "pass filters and QC" ${output}.log | awk '{ print $4 }' )
    printf "$samples samples\n$variants variants\n\n"
    
    # now get the number of variants dropped for missingess
    vargeno=$(  grep "removed due to missing genotype data" ${output}.log | awk '{print $1;}' )
    variants=$( grep "pass filters and QC" ${output}.log | awk '{ print $1 }' )
    printf "Removed $vargeno variants for missingness>0.01.\n\n$samples samples\n$variants variants\n\n"
    
    output=${output}_nogeno
    # get number of genotyped variants dropped from imputed set
    droppedgeno=$(( $variants - $( grep "variants remaining" ${output}.log | awk '{ print $2 }' ) ))
    # get number of unimputed genotyped variants
    justgeno=$(( $( wc -l < ${output_folder}/${i}_genotyped_variants.txt ) - $droppedgeno ))
    variants=$( grep "pass filters and QC" ${output}.log | awk '{ print $1 }' )
    
    # check whether there were any other variants that needed to be removed
    if [ -e "${output}2.bim" ]; 
    then 
      samepos=$(( $variants - $( grep "variants remaining" ${output}2.log | awk '{ print $2 }' ) ))
      printf "Replaced ${droppedgeno} imputed genotyped SNPs, removed $samepos imputed variants which were at the same position as genotyped variants, and included ${justgeno} unimputed genotyped SNPs.\n\n"
    else
      printf "Replaced ${droppedgeno} imputed genotyped SNPs and included ${justgeno} unimputed genotyped SNPs.\n\n"
    fi
    
    # get final imp/geno merged values
    output=${output}_merged
    samples=$( grep "pass filters and QC" ${output}.log | awk '{ print $4 }' )
    variants=$( grep "pass filters and QC" ${output}.log | awk '{ print $1 }' )
    printf "$samples samples\n$variants variants\n\n"
  
    # get hwe and maf filter numbers
    output=${output}_hwe6_maf01
    varhwe=$(  grep "removed due to Hardy-Weinberg exact test" ${output}.log | awk '{print $2;}' )
    varmaf=$(  grep "removed due to minor allele threshold(s)" ${output}.log | awk '{print $1;}' )
    samples=$( grep "pass filters and QC" ${output}.log | awk '{ print $4 }' )
    variants=$( grep "pass filters and QC" ${output}.log | awk '{ print $1 }' )
    printf "Removed ${varhwe} SNPs for HWE p<1e-6, ${varmaf} for MAF<0.01.\n\n$samples samples\n$variants variants\n\n"

    # heterozygosity check
    if [ -f "${output}_nohetout.bim" ];
    then
        output=${output}_nohetout
        sampleshet=$( grep "people remaining" ${output}.log | awk '{print $2;}' )
        printf '%s\n' "Removed $sampelshet individuals as heterozygosity outliers."
        ## get the resulting number of samples and variants
        variants=$( grep "pass filters and QC" ${output}.log | awk '{print $1;}' )
        samples=$( grep "pass filters and QC" ${output}.log | awk '{print $4;}' )
        printf "\n\n$samples samples\n$variants variants\n\n"
    else
        printf "No heterozygosity outliers identified.\n\n$samples samples\n$variants variants.\n\n"
    fi

  fi
done
