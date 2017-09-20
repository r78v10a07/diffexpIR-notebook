#!/bin/bash

sample=$1
sDir=$2
index=$3

star="/panfs/pan1.be-md.ncbi.nlm.nih.gov/alt_splicing/bioNotebook-DiffExpIR/bin/STAR"
samtools="/usr/local/samtools/1.3.1/bin/samtools"

if [ ! -e "${sample}_sorted.bam" ]
then
  mkdir -p $sample
  cd $sample
  samples=""
  if [ -e "$sDir/${sample}_1.fastq.gz" ]
  then
    samples="$sDir/${sample}_1.fastq.gz"
  fi
  if [ -e "$sDir/${sample}_2.fastq.gz" ]
  then
    samples="$sDir/${sample}_1.fastq.gz $sDir/${sample}_2.fastq.gz"
  fi
  if [ -z "$samples" ]
  then
    echo "ERROR: No sample"
    exit -1
  fi
  $star --readFilesCommand zcat --runThreadN 16 --genomeDir $index --outSAMtype BAM Unsorted --outStd BAM_Unsorted --readFilesIn ${samples} > ${sample}.bam
  $samtools flagstat ${sample}.bam > ${sample}.stats
  $samtools sort --threads 16 -o ${sample}_sorted.bam ${sample}.bam
  $samtools index ${sample}_sorted.bam ${sample}_sorted.bam.bai
  mv ${sample}_sorted.bam ${sample}.stats ${sample}_sorted.bam.bai ../
  cd ..
  rm -rf $sample
fi
exit 0
