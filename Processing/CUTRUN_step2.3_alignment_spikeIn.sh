#!/bin/bash

ml Bowtie2/2.4.1-GCCcore-8.3.0
ml SAMtools/1.10-GCCcore-8.3.0
ml BEDTools/2.29.2-GCC-9.3.0


if [ ! -e "$outPath" ]; then
    mkdir -p $outPath
fi

## alignment
bowtie2 --local --very-sensitive-local --no-unal --no-mixed --dovetail --no-discordant --phred33 -I 10 -X 700 -p $cores -x $refSpikeIn -1 $inPath/R1.fastq.gz -2 $inPath/R2.fastq.gz -S $outPath/bowtie2_align_spikeIn.sam 

seqDepth=`samtools view -F 0x04 $outPath/bowtie2_align_spikeIn.sam | wc -l`
echo $seqDepth >$outPath/bowtie2_align_spikeIn.seqDepth
echo $seqDepth

if [[ "$seqDepth" -gt "1" ]]; then
    scale_factor=`echo "20000 / $seqDepth" | bc -l`
    echo $scale_factor
    bedtools genomecov -bg -scale $scale_factor -i $outPath/bowtie2_align.fragments.bed -g $chromSize >$outPath/bowtie2_align_norm.fragments.bedgraph
fi


		
## sort by coordinate
$picardCMD SortSam I=$outPath/bowtie2_align.sam O=$outPath/bowtie2_align.sorted.sam SORT_ORDER=coordinate


## ========== Stringent Setting =======
## mark duplicates
$picardCMD MarkDuplicates I=$outPath/bowtie2_align.sorted.sam O=$outPath/bowtie2_align.sorted.dupMarked.sam METRICS_FILE=$outPath/picard.dupMark.txt

## remove duplicates
$picardCMD MarkDuplicates I=$outPath/bowtie2_align.sorted.sam O=$outPath/bowtie2_align.sorted.rmDup.sam METRICS_FILE=$outPath/picard.rmDup.txt REMOVE_DUPLICATES=true

## filter by the fragment length
awk '
BEGIN { FS="\t"; SIZE=120; S2=SIZE*SIZE } 
/^@/ { print $0; next }
{ if ($9*$9 < S2) print $0} 
' $outPath/bowtie2_align.sorted.rmDup.sam >$outPath/bowtie2_align.sorted.rmDup.fragLenFilter.sam




