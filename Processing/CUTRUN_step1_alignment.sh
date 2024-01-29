#!/bin/bash

ml Bowtie2/2.4.1-GCCcore-8.3.0
ml picard/2.21.6-Java-11


if [ ! -e "$outPath" ]; then
    mkdir -p $outPath
fi


picardCMD="java -jar $EBROOTPICARD/picard.jar"

## alignment
bowtie2 --local --very-sensitive-local --no-unal --no-mixed --dovetail --no-discordant --phred33 -I 10 -X 700 -p $cores -x $ref -1 $inPath/R1.fastq.gz -2 $inPath/R2.fastq.gz -S $outPath/bowtie2_align.sam 

## sort by coordinate
$picardCMD SortSam I=$outPath/bowtie2_align.sam O=$outPath/bowtie2_align.sorted.sam SORT_ORDER=coordinate


## ========== Stringent Setting =======
## mark duplicates
$picardCMD MarkDuplicates I=$outPath/bowtie2_align.sorted.sam O=$outPath/bowtie2_align.sorted.dupMarked.sam METRICS_FILE=$outPath/picard.dupMark.txt

## remove duplicates
$picardCMD MarkDuplicates I=$outPath/bowtie2_align.sorted.sam O=$outPath/bowtie2_align.sorted.rmDup.sam METRICS_FILE=$outPath/picard.rmDup.txt REMOVE_DUPLICATES=true

## filter by the fragment length
# awk '
# BEGIN { FS="\t"; SIZE=120; S2=SIZE*SIZE } 
# /^@/ { print $0; next }
# { if ($9*$9 < S2) print $0} 
# ' $outPath/bowtie2_align.sorted.rmDup.sam >$outPath/bowtie2_align.sorted.rmDup.fragLenFilter.sam




