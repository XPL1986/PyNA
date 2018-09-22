#!/bin/sh

# Pipeline for Calling variants in DNAseq
# Adapted from Broad Institute at https://software.broadinstitute.org/gatk/documentation/article.php?id=3891

args=("$@")
pid=${args[0]}
HOME=$(pwd)
export samtools=${HOME}/bin/samtools-1.8/bin/samtools

mkdir ${HOME}/data/${pid}

#cd ${HOME}/data/${pid}
#input=/Volumes/qmu/NG2016/rGBM_DNA/${pid}.bam
#output=${HOME}/data/${pid}/${pid}.sort.bam
#samtools sort $input -@ 8 -o $output
#samtools index ${pid}.sort.bam


# 2. Add read groups, sort, mark duplicates, and create index
cd ${HOME}
picardDir=${HOME}/bin/picard/picard.jar
outputDir=${HOME}/data/${pid}/step2 && mkdir $outputDir

#input=${HOME}/data/${pid}/${pid}.sort.bam
input=/Volumes/qmu/NG2016/rGBM_DNA/${pid}.bam
output=${outputDir}/rg_added_sorted.bam
java -jar $picardDir AddOrReplaceReadGroups I=$input O=$output SO=coordinate RGID=id RGLB=library RGPL=platform RGPU=machine RGSM=sample VALIDATION_STRINGENCY=LENIENT

input=$output
output=${outputDir}/dedupped.bam
java -jar $picardDir MarkDuplicates I=$input O=$output  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics

# 3. Split'N'Trim and reassign mapping qualities (skip)
# 4. Indel Realignment (skip)
# 5. Base Recalibration (skip)


# 6. Variant calling ||GATK requires Java 1.8 and reference files||
export JAVA_HOME=`/usr/libexec/java_home -v 1.8`
#ref=${HOME}/bin/hg38/hg38.fa
ref=${HOME}/bin/GRCh37/reference.fa
gatkDir=${HOME}/bin/gatk/GenomeAnalysisTK.jar

outputDir=${HOME}/data/${pid}/step6 && mkdir $outputDir

input=$output
output=${outputDir}/output.vcf
java -Xmx8g -jar $gatkDir -T HaplotypeCaller -R $ref -I $input -o $output -stand_call_conf 20.0


# 7. Variant filtering ||GATK requires Java 1.8 and reference files||
outputDir=${HOME}/data/${pid}/step7 && mkdir $outputDir

input=$output
output=${outputDir}/output.vcf
java -jar $gatkDir -T VariantFiltration -R $ref -V $input -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o $output



