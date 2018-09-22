#!/bin/sh

args=("$@")
pid=${args[0]}
HOME=$(pwd)


snpeff=${HOME}/bin/snpEff/SnpEff.jar
sconfig=${HOME}/bin/snpEff/snpEff.config
snpsift=${HOME}/bin/snpEff/SnpSift.jar
dbSnp=${HOME}/bin/dbSnp/All_20180418.vcf.gz

## Annotate using SnpEff and fields from another VCF file (e.g. dbSnp, 1000 Genomes projects, ClinVar, ExAC, etc.).
#cd ${HOME}/data/${pid}/step6/
#java -Xmx4g -jar $snpeff -c $sconfig -v GRCh37.75 output.vcf > output_ann.vcf
#java -jar $snpsift annotate $dbSnp output_ann.vcf > output_ann_dbsnp.vcf



# Annotate using SnpEff
outputDir=${HOME}/data/${pid}/snpsift && mkdir $outputDir
input=${HOME}/data/${pid}/step7/output.vcf
output=${outputDir}/${pid}_annotated.vcf

cd $outputDir
java -Xmx4g -jar $snpeff -c $sconfig -v GRCh37.75 $input > $output
cd $HOME



# Annotate using fields from another VCF file (e.g. dbSnp, 1000 Genomes projects, ClinVar, ExAC, etc.).

input=$output
output=${outputDir}/${pid}_annotated_dbsnp.vcf 

java -jar $snpsift annotate $dbSnp $input > $output



# Applying filters
input=$output
output=${outputDir}/${pid}_annotated_dbsnp_filtered.vcf

cat $input | java -jar $snpsift filter "(FILTER = 'PASS')" > $output

