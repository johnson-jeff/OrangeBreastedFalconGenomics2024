#!/bin/bash

###### Overview ######

## this script will quantify annotated variants generated for each individual sample using snpEff
## assumes snpEff was run with single sample VCF (not a combined sample VCF) and ANN= (or EFF=) annotations are present for each variant

## assumes that each *ann.VCF file is zipped with bgzip (if not, replace all the zcat commands with cat)
## if executing on maxos, replace zcat with gzcat

## PRIOR to running script, prepare name.list with VCF file names in a single column

## example command line
## bash ./snpEff_summary.sh genic.annotated.VCFs.list > OUT.txt

##### Command Line arguments #####
#$ $1 name of sample list file, eg, genic.annotated.VCFs.list


### note that the command used below to count the number of each variant type only counts each variant
### once (eg, zcat ${I} | grep HIGH | wc -l) even if multiple transcripts with the same variant type are present at that particular site
### use the following example to count total number of transcripts with a particular variant type:
###
### ## example for quantifying the number of HIGH Impact variant transcripts
### zcat ${i} | cut -f 8 | tr ";" "\n" | grep ^ANN= | cut -f 2 -d = | tr "," "\n" | grep HIGH | wc -l
### 
### ## example for quantifying the number of heterozygote HIGH Impact variant transcripts
### zcat ${i} | grep "0/1:" | cut -f 8 | tr ";" "\n" | grep ^ANN= | cut -f 2 -d = | tr "," "\n" | grep HIGH | wc -l)


###### MAIN SCRIPT ######

NAMES=( `cat "$1" `)

for i in "${NAMES[@]}"
do

echo "##############################################"
echo "################  ${i}  ################"
echo "##############################################"

### for estimating genetic load for each individual
# To get the allele counts from the genotypes #
# 0/0 = 0 allele ; 0/1 = 1 allele; 1/1 = 2 alleles

######### LOF #########

## number of LOF genic variants in individual i
LOF=$(zcat ${i} | grep -v WARNING | grep LOF | wc -l)
step1=$(python -c "print($LOF)")
echo "$step1 = number of LOF genic variants in ${i}"

## number of LOF homozygous genic variants & mutations in individual i
LOF_homo1a=$(zcat ${i} | grep -v WARNING | grep LOF | grep "1/1:" | wc -l)
LOF_homo1b=$(zcat ${i} | grep -v WARNING | grep LOF | grep "1|1:" | wc -l)
step2a=$(python -c "print($LOF_homo1a + $LOF_homo1b)")
echo "$step2a = number of LOF homozygous 1/1 + 1|1 genic variants in ${i}"
step2b=$(python -c "print($step2a * 2)")
echo "$step2b = number of LOF homozygous 1/1 + 1|1 genic mutations in ${i}"

LOF_homo2a=$(zcat ${i} | grep -v WARNING | grep LOF | grep "0/0:" | wc -l)
LOF_homo2b=$(zcat ${i} | grep -v WARNING | grep LOF | grep "0|0:" | wc -l)
step2c=$(python -c "print($LOF_homo2a + $LOF_homo2b)")
echo "$step2c = number of LOF homozygous 0/0 + 0|0 genic variants in ${i}"
step2d=$(python -c "print($LOF_homo1a + $LOF_homo1b + $LOF_homo2a + $LOF_homo2b)")
echo "$step2d = number of LOF homozygous 0/0 + 0|0 + 1/1 + 1|1 genic variants in ${i}"

## number of LOF heterozygous genic variants & mutations in individual i
LOF_het1a=$(zcat ${i} | grep -v WARNING | grep LOF | grep "0/1:" | wc -l)
LOF_het1b=$(zcat ${i} | grep -v WARNING | grep LOF | grep "0|1:" | wc -l)
step2e=$(python -c "print($LOF_het1a + $LOF_het1b)")
echo "$step2e = number of LOF heterozygous 0/1 + 0|1 genic variants & mutations in ${i}"

######### HIGH Impact Variant #########

## number of HIGH Impact genic variants in individual i
HIGH=$(zcat ${i} | grep -v WARNING | grep HIGH | wc -l)
step3=$(python -c "print($HIGH)")
echo "$step3 = number of HIGH Impact genic variants in ${i}"

## number of HIGH homozygous genic variants & mutations in individual i
HIGH_homo1a=$(zcat ${i} | grep -v WARNING | grep HIGH | grep "1/1:" | wc -l)
HIGH_homo1b=$(zcat ${i} | grep -v WARNING | grep HIGH | grep "1|1:" | wc -l)
step4a=$(python -c "print($HIGH_homo1a + $HIGH_homo1b)")
echo "$step4a = number of HIGH Impact homozygous 1/1 + 1|1 genic variants in ${i}"
step4b=$(python -c "print($step4a * 2)")
echo "$step4b = number of HIGH Impact homozygous 1/1 + 1|1 genic mutations in ${i}"

HIGH_homo2a=$(zcat ${i} | grep -v WARNING | grep HIGH | grep "0/0:" | wc -l)
HIGH_homo2b=$(zcat ${i} | grep -v WARNING | grep HIGH | grep "0|0:" | wc -l)
step4c=$(python -c "print($HIGH_homo2a + $HIGH_homo2b)")
echo "$step4c = number of HIGH Impact homozygous 0/0 + 0|0 genic variants in ${i}"
step4d=$(python -c "print($HIGH_homo1a + $HIGH_homo1b + $HIGH_homo2a + $HIGH_homo2b)")
echo "$step4d = number of HIGH Impact homozygous 0/0 + 0|0 + 1/1 + 1|1 genic variants in ${i}"

## number of HIGH heterozygous genic variants & mutations in individual i
HIGH_het1a=$(zcat ${i} | grep -v WARNING | grep HIGH | grep "0/1:" | wc -l)
HIGH_het1b=$(zcat ${i} | grep -v WARNING | grep HIGH | grep "0|1:" | wc -l)
step4e=$(python -c "print($HIGH_het1a + $HIGH_het1b)")
echo "$step4e = number of HIGH Impact heterozygous 0/1 + 0|1 genic variants & mutations in ${i}"

######### MODERATE Impact Variant #########

## number of MOD Impact genic variants in individual i
MOD=$(zcat ${i} | grep -v WARNING | grep MODERATE | wc -l)
step5=$(python -c "print($MOD)")
echo "$step5 = number of MOD Impact genic variants in ${i}"

## number of MOD homozygous genic variants & mutations in individual i
MOD_homo1a=$(zcat ${i} | grep -v WARNING | grep MODERATE | grep "1/1:" | wc -l)
MOD_homo1b=$(zcat ${i} | grep -v WARNING | grep MODERATE | grep "1|1:" | wc -l)
step6a=$(python -c "print($MOD_homo1a + $MOD_homo1b)")
echo "$step6a = number of MOD Impact homozygous 1/1 + 1|1 genic variants in ${i}"
step6b=$(python -c "print($step6a * 2)")
echo "$step6b = number of MOD Impact homozygous 1/1 + 1|1 genic mutations in ${i}"

MOD_homo2a=$(zcat ${i} | grep -v WARNING | grep MODERATE | grep "0/0:" | wc -l)
MOD_homo2b=$(zcat ${i} | grep -v WARNING | grep MODERATE | grep "0|0:" | wc -l)
step6c=$(python -c "print($MOD_homo2a + $MOD_homo2b)")
echo "$step6c = number of MOD Impact homozygous 0/0 + 0|0 genic variants in ${i}"
step6d=$(python -c "print($MOD_homo1a + $MOD_homo1b + $MOD_homo2a + $MOD_homo2b)")
echo "$step6d = number of MOD Impact homozygous 0/0 + 0|0 + 1/1 + 1|1 genic variants in ${i}"

## number of MOD heterozygous genic variants & mutations in individual i
MOD_het1a=$(zcat ${i} | grep -v WARNING | grep MODERATE | grep "0/1:" | wc -l)
MOD_het1b=$(zcat ${i} | grep -v WARNING | grep MODERATE | grep "0|1:" | wc -l)
step6e=$(python -c "print($MOD_het1a + $MOD_het1b)")
echo "$step6e = number of MOD Impact heterozygous 0/1 + 0|1 genic variants & mutations in ${i}"

######### LOW Impact Variant #########

## number of LOW Impact genic variants in individual i
LOW=$(zcat ${i} | grep -v WARNING | grep LOW | wc -l)
step7=$(python -c "print($LOW)")
echo "$step7 = number of LOW Impact genic variants in ${i}"

## number of LOW homozygous genic variants & mutations in individual i
LOW_homo1a=$(zcat ${i} | grep -v WARNING | grep LOW | grep "1/1:" | wc -l)
LOW_homo1b=$(zcat ${i} | grep -v WARNING | grep LOW | grep "1|1:" | wc -l)
step8a=$(python -c "print($LOW_homo1a + $LOW_homo1b)")
echo "$step8a = number of LOW Impact homozygous 1/1 + 1|1 genic variants in ${i}"
step8b=$(python -c "print($step8a * 2)")
echo "$step8b = number of LOW Impact homozygous 1/1 + 1|1 genic mutations in ${i}"

LOW_homo2a=$(zcat ${i} | grep -v WARNING | grep LOW | grep "0/0:" | wc -l)
LOW_homo2b=$(zcat ${i} | grep -v WARNING | grep LOW | grep "0|0:" | wc -l)
step8c=$(python -c "print($LOW_homo2a + $LOW_homo2b)")
echo "$step8c = number of LOW Impact homozygous 0/0 + 0|0 genic variants in ${i}"
step8d=$(python -c "print($LOW_homo1a + $LOW_homo1b + $LOW_homo2a + $LOW_homo2b)")
echo "$step8d = number of LOW Impact homozygous 0/0 + 0|0 + 1/1 + 1|1 genic variants in ${i}"

## number of LOW heterozygous genic variants & mutations in individual i
LOW_het1a=$(zcat ${i} | grep -v WARNING | grep LOW | grep "0/1:" | wc -l)
LOW_het1b=$(zcat ${i} | grep -v WARNING | grep LOW | grep "0|1:" | wc -l)
step8e=$(python -c "print($LOW_het1a + $LOW_het1b)")
echo "$step8e = number of LOW Impact heterozygous 0/1 + 0|1 genic variants & mutations in ${i}"

######### NoImpact Variant #########

## number of NoImpact genic variants in individual i
NoImpact=$(zcat ${i} | grep -v WARNING | grep MODIFIER | wc -l)
step9=$(python -c "print($NoImpact)")
echo "$step9 = number of NoImpact genic variants in ${i}"

## number of NoImpact homozygous genic variants & mutations in individual i
NoImpact_homo1a=$(zcat ${i} | grep -v WARNING | grep MODIFIER | grep "1/1:" | wc -l)
NoImpact_homo1b=$(zcat ${i} | grep -v WARNING | grep MODIFIER | grep "1|1:" | wc -l)
step10a=$(python -c "print($NoImpact_homo1a + $NoImpact_homo1b)")
echo "$step10a = number of NoImpact homozygous 1/1 + 1|1 genic variants in ${i}"
step10b=$(python -c "print($step10a * 2)")
echo "$step10b = number of NoImpact homozygous 1/1 + 1|1 genic mutations in ${i}"

NoImpact_homo2a=$(zcat ${i} | grep -v WARNING | grep MODIFIER | grep "0/0:" | wc -l)
NoImpact_homo2b=$(zcat ${i} | grep -v WARNING | grep MODIFIER | grep "0|0:" | wc -l)
step10c=$(python -c "print($NoImpact_homo2a + $NoImpact_homo2b)")
echo "$step10c = number of NoImpact homozygous 0/0 + 0|0 genic variants in ${i}"
step10d=$(python -c "print($NoImpact_homo1a + $NoImpact_homo1b + $NoImpact_homo2a + $NoImpact_homo2b)")
echo "$step10d = number of NoImpact homozygous 0/0 + 0|0 + 1/1 + 1|1 genic variants in ${i}"

## number of NoImpact heterozygous genic variants & mutations in individual i
NoImpact_het1a=$(zcat ${i} | grep -v WARNING | grep MODIFIER | grep "0/1:" | wc -l)
NoImpact_het1b=$(zcat ${i} | grep -v WARNING | grep MODIFIER | grep "0|1:" | wc -l)
step10e=$(python -c "print($NoImpact_het1a + $NoImpact_het1b)")
echo "$step10e = number of NoImpact heterozygous 0/1 + 0|1 genic variants & mutations in ${i}"

######### Total number of genic Variants #########

## Total number of genic variants in individual i
TOTAL=$(zcat ${i} | grep -v WARNING | wc -l)
step11=$(python -c "print($TOTAL)")
echo "$step11 = TOTAL number of genic variants in ${i}"

TOTAL_homo1a=$(zcat ${i} | grep -v WARNING | grep "1/1:" | wc -l)
TOTAL_homo1b=$(zcat ${i} | grep -v WARNING | grep "1|1:" | wc -l)
step12a=$(python -c "print($TOTAL_homo1a + $TOTAL_homo1b)")
echo "$step12a = TOTAL number of homozygous 1/1 + 1|1 genic variants in ${i}"
step12b=$(python -c "print($step12a * 2)")
echo "$step12b = TOTAL number of homozygous 1/1 + 1|1 genic mutations in ${i}"

TOTAL_homo2a=$(zcat ${i} | grep -v WARNING | grep "0/0:" | wc -l)
TOTAL_homo2b=$(zcat ${i} | grep -v WARNING | grep "0|0:" | wc -l)
step12c=$(python -c "print($TOTAL_homo2a + $TOTAL_homo2b)")
echo "$step12c = Total number of homozygous 0/0 + 0|0 genic variants in ${i}"
step12d=$(python -c "print($TOTAL_homo1a + $TOTAL_homo1b + $TOTAL_homo2a + $TOTAL_homo2b)")
echo "$step12d = Total number of homozygous 0/0 + 0|0 + 1/1 + 1|1 genic variants in ${i}"

## number of LOW heterozygous genic variants & mutations in individual i
TOTAL_het1a=$(zcat ${i} | grep -v WARNING | grep "0/1:" | wc -l)
TOTAL_het1b=$(zcat ${i} | grep -v WARNING | grep "0|1:" | wc -l)
step12e=$(python -c "print($TOTAL_het1a + $TOTAL_het1b)")
echo "$step12e = TOTAL number of heterozygous 0/1 + 0|1 genic variants & mutations in ${i}"

## Total number of genic mutations in individual i
step12f=$(python -c "print($step12b + $step12e)")
echo "$step12f = TOTAL number of genic mutations in ${i}"

######### Synonymous Variants #########

## number of synonymous genic variants in individual i
SYN1=$(zcat ${i} | grep -v WARNING | grep synonymous | wc -l)
step13=$(python -c "print($SYN1)")
echo "$step13 = number of synonymous genic variants in ${i}"

## number of synonymous homozygous genic variants in individual i
SYN_homo1a=$(zcat ${i} | grep -v WARNING | grep synonymous | grep "1/1:" | wc -l)
SYN_homo1b=$(zcat ${i} | grep -v WARNING | grep synonymous | grep "1|1:" | wc -l)
step14a=$(python -c "print($SYN_homo1a + $SYN_homo1b)")
echo "$step14a = number of synonymous homozygous 1/1 + 1|1 genic variants in ${i}"

SYN_homo2a=$(zcat ${i} | grep -v WARNING | grep synonymous | grep "0/0:" | wc -l)
SYN_homo2b=$(zcat ${i} | grep -v WARNING | grep synonymous | grep "0|0:" | wc -l)
step14b=$(python -c "print($SYN_homo2a + $SYN_homo2b)")
echo "$step14b = number of synonymous homozygous 0/0 + 0|0 genic variants in ${i}"
step14c=$(python -c "print($SYN_homo1a + $SYN_homo1b + $SYN_homo2a + $SYN_homo2b)")
echo "$step14c = number of synonymous homozygous 0/0 + 0|0 + 1/1 + 1|1 genic variants in ${i}"

## number of synonymous heterozygous genic variants in individual i
SYN_het1a=$(zcat ${i} | grep -v WARNING | grep synonymous | grep "0/1:" | wc -l)
SYN_het1b=$(zcat ${i} | grep -v WARNING | grep synonymous | grep "0|1:" | wc -l)
step14d=$(python -c "print($SYN_het1a + $SYN_het1b)")
echo "$step14d = number of synonymous heterozygous 0/1 + 0|1 genic mutations in ${i}"


done