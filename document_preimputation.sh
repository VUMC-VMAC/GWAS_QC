output=$1

set -e

# get number of samples and variants to begin
variants=$( grep "loaded from .bim" ${output}.log | awk '{print $1;}' )
samples=$( grep "loaded from .fam" ${output}.log | awk '{print $1;}' )

printf "$samples samples
$variants variants\n\n"

# geno and maf
output=${output}_geno05_maf01
## get the number of variants dropped for initial variant filters
vargeno=$(  grep "removed due to missing genotype data" ${output}.log | awk '{print $1;}' )
varmaf=$(  grep "removed due to minor allele threshold(s)" ${output}.log | awk '{print $1;}' )

printf '%s\n' "Removed $vargeno SNPs for >5% missingness and $varmaf SNPs for MAF <0.01."

## get the resulting number of samples and variants
variants=$( grep "pass filters and QC" ${output}.log | awk '{print $1;}' )
samples=$( grep "pass filters and QC" ${output}.log | awk '{print $4;}' )

printf "\n$samples samples
$variants variants\n\n"

# missingness
output=${output}_mind01
## get the number of samples removed for missingness
samplesgeno=$( grep "removed due to missing genotype data" ${output}.log | awk '{print $1;}' )
printf '%s\n' "Removed $samplesgeno people for >1% missingness."
## get the resulting number of samples and variants
variants=$( grep "pass filters and QC" ${output}.log | awk '{print $1;}' )
samples=$( grep "pass filters and QC" ${output}.log | awk '{print $4;}' )
printf "\n$samples samples
$variants variants\n\n"

# relatedness
output=${output}_norelated
## get the number of samples dropped for relatedness
samplesrelated=$(( $samples - $( grep "people remaining" ${output}.log | awk '{print $2;}' ) ))
printf '%s\n' "Removed $samplesrelated related individuals.
"
## get the resulting number of samples and variants
variants=$( grep "pass filters and QC" ${output}.log | awk '{print $1;}' )
samples=$( grep "pass filters and QC" ${output}.log | awk '{print $4;}' )

printf "$samples samples
$variants variants\n\n"

# sex check
## test whether sex in the fam file was updated
if test -f ${output}_sex.fam ; 
then 
    output=${output}_sex
    # test whether the people were dropped
    if test -f ${output}_nomismatchedsex.fam ; 
    then 
	output=${output}_nomismatchedsex
	## get the number of samples which fail sex check
	samplessex=$(($samples - $( grep "people remaining" ${output}.log | awk '{print $2;}' )))
	printf '%s\n' "Removed $sampelssex individuals with mismatched sex."
	## get the resulting number of samples and variants
	variants=$( grep "pass filters and QC" ${output}.log | awk '{print $1;}' )
	samples=$( grep "pass filters and QC" ${output}.log | awk '{print $4;}' )
	printf "\n\n$samples samples\n$variants variants\n"
    else
	printf "Removed 0 individuals with mismatched sex.\n\n$samples samples\n$variants variants.\n\n"
    fi
else
    printf "Removed 0 individuals with mismatched sex.\n\n$samples samples\n$variants variants.\n\n"
fi

# restrict to autosomes
output_last=$output
output=${output}_keep_autosomes
## get the number of variants removed 
# total
non_autosomal_var=$(( $(wc -l <${output_last}.bim) - $(wc -l <${output}.bim) ))
# unmapped (chr=0)
unmapped_var=$( awk '{ if($1 == 0) print }' $output_last.bim | wc -l)
# sex variants (chr=23,24,25)
sex_var=$( awk '{ if($1 == 23 || $1 == 24 || $1 == 25) print }' $output_last.bim | wc -l)
# mitochondrial (chr=26)
mito_var=$(awk '{ if($1 == 26) print }' $output_last.bim | wc -l)

printf '%s\n' "Removed $non_autosomal_var variants total ($unmapped_var unmapped, $sex_var sex chr, $mito_var mitochondrial)."
printf "\n"

## get the resulting number of samples and variants
variants=$( grep "pass filters and QC" ${output}.log | awk '{print $1;}' )
samples=$( grep "pass filters and QC" ${output}.log | awk '{print $4;}' )
printf "\n$samples samples\n$variants variants\n\n"

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

# HWE filter
output=${output}_hwe6
varhwe=$(  grep "removed due to Hardy-Weinberg exact test" ${output}.log | awk '{print $2;}' )
printf "Removed $varhwe SNPs for HWE p<1e-6.\n"
## get the resulting number of samples and variants
variants=$( grep "pass filters and QC" ${output}.log | awk '{print $1;}' )
samples=$( grep "pass filters and QC" ${output}.log | awk '{print $4;}' )
printf "\n$samples samples
$variants variants\n\n"

# palindromic variants
output=${output}_nopal
variants_old=$variants
## identify the number of variants remaining
variants=$( grep "pass filters and QC" ${output}.log | awk '{print $1;}' )
samples=$( grep "pass filters and QC" ${output}.log | awk '{print $4;}' )
## log
varspalindromic=$(( $variants_old - $variants ))

#get list of variants to remove (those with non-standard allele codes or which cannot be lifted)
if [ -f ${output}_toliftover.txt ]; 
then
    varsfaillift=$( wc -l < ${output}_lift_issue_snps_todrop.txt )
    output=${output}_noprobSNPs_chr_bplifted
else 
    varsfaillift=0
fi

if [ "$( wc -l < ${output}_samepos_vars.txt )" -gt 0 ];
then
    varssamepos=$( wc -l < ${output}_samepos_vars.txt )
    output=${output}_nosamepos
else
    varssamepos=0
fi

# get the number of variants that didn't match reference
output_path=${output%/*}
varsmismatchref=$( cat ${output_path}/Exclude-* | wc -l )

printf "Removing $varspalindromic palindromic, $varsfaillift variants which fail liftOver, $varssamepos same position SNPs, and $varsmismatchref SNPs not matching the reference panel.\n"

## get the resulting number of samples and variants
variants=$( grep "pass filters and QC" ${output}*-updated.log | awk '{print $1;}' )
samples=$( grep "pass filters and QC" ${output}*-updated.log | awk '{print $4;}' )
printf "\n$samples samples
$variants variants\n"
