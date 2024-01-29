#!/bin/bash
core=4

rsemDir="Software/RSEM/bin"
refDir="SupplementaryData/hg38/fasta"
gecodeDir="SupplementaryData/hg38/GENCODE/"
rsemIndex="SupplementaryData/hg38/RSEMIndex_v21"
starPath="STAR/2.7.7a-GCC-10.2.0/bin"
inPath="Unaligned/Project_lhanafi/${sampleName}"
suppPath="SupplementaryData/hg38/exon"

ml STAR/2.7.7a-GCC-10.2.0 #STAR/2.7.3a-foss-2016b


## generate reference genome index
## $rsemDir/rsem-prepare-reference --gtf  $gecodeDir/gencode.v33.annotation.gtf  --star --star-path $starPath -p 20 $refDir/hg38.fa $rsemIndex/hg38
$rsemDir/rsem-prepare-reference --gtf  $gecodeDir/gencode.v21.annotation.gtf  --star --star-path $starPath -p 6 $refDir/hg38.fa $rsemIndex/hg38

# optional, merge multiple lanes per read end
mkdir -p $outPath/$sampleName/FASTQ/
mkdir -p $outPath/$sampleName/RSEM/
mkdir -p $outPath/$sampleName/log/

# cat $inPath/*_R1_001.fastq.gz >$outPath/$sampleName/FASTQ/R1.fastq.gz
# cat $inPath/*_R2_001.fastq.gz >$outPath/$sampleName/FASTQ/R2.fastq.gz

## RSEM alignment
$rsemDir/rsem-calculate-expression -p $core --star --star-path $starPath --estimate-rspd --append-names --output-genome-bam --star-gzipped-read-file --paired-end $fastqPath/${sampleName}/FASTQ/R1.fastq.gz $fastqPath/${sampleName}/FASTQ/R2.fastq.gz $rsemIndex/hg38 $outPath/${sampleName}/RSEM > $outPath/${sampleName}/log/${sampleName}.log 2>&1

## Overlap with exon regions
ml BEDTools/2.30.0-GCC-10.2.0  #BEDTools/2.29.1-foss-2016b
bedtools intersect -abam $outPath/$sampleName/RSEM.genome.bam -b $suppPath/hg38_exon_v32_codeExon.bed >$outPath/$sampleName/RSEM.genome.exon.bam

## Summary
ml SAMtools/1.14-GCC-11.2.0 #SAMtools/1.10-foss-2016b
echo "Aligned reads to exon:" >>$outPath/$sampleName/RSEM.summary
samtools view -F 4 $outPath/$sampleName/RSEM.genome.exon.bam | awk '{print $1}' | uniq | wc -l >>$outPath/$sampleName/RSEM.summary
