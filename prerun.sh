#!/bin/sh

# Pipeline for Calling variants in RNAseq
# Adapted from Broad Institute at https://software.broadinstitute.org/gatk/documentation/article.php?id=3891

args=("$@")
pid=${args[0]}
HOME=$(pwd)



#--outSAMunmapped Within

# Ensure quality control via FastQC and MultiQC
export fastqc=${HOME}/bin/fastqc/FastQC.app/Contents/MacOS/fastqc
outputDir=${HOME}/data/${pid}/fastqc && mkdir $outputDir

input_1=${HOME}/data/${pid}/${pid}_1.fastq.gz
input_2=${HOME}/data/${pid}/${pid}_2.fastq.gz
#$fastqc $input_1 $input_2 --outdir=$outputDir

# pip install multiqc
inputDir=$outputDir
outputDir=${HOME}/data/${pid}/multiqc && mkdir $outputDir
cd $outputDir && multiqc $inputDir && cd $HOME



# 1. Mapping to the reference
export STAR=${HOME}/bin/STAR-master/bin/Linux_x86_64/STAR

genomeDir=${HOME}/bin/hg38/genomeDir && mkdir $genomeDir
cd $genomeDir
ref=${HOME}/bin/hg38/hg38.fa
#STAR --runMode genomeGenerate --genomeDir $genomeDir --genomeFastaFiles $ref  --runThreadN 12


outputDir=${HOME}/data/${pid}/step1 && mkdir $outputDir
runDir=${outputDir}/1pass && mkdir $runDir
cd $runDir
#STAR --genomeDir $genomeDir --readFilesIn $input_1 $input_2 --readFilesCommand gunzip -c --runThreadN 12 


genomeDir=${HOME}/data/${pid}/step1/hg38_2pass && mkdir $genomeDir
cd $genomeDir
#STAR --runMode genomeGenerate --genomeDir $genomeDir --genomeFastaFiles $ref --sjdbFileChrStartEnd ${runDir}/SJ.out.tab --sjdbOverhang 75 --runThreadN 12


runDir=${HOME}/data/${pid}/step1/2pass && mkdir $runDir && cd $runDir
#STAR --genomeDir $genomeDir--readFilesIn $input_1 $input_2 --readFilesCommand gunzip --runThreadN 12

cd $HOME



# 1.5 Conversion to bam, sort, index
input=${HOME}/data/${pid}/step1/2pass/Aligned.out.sam
