#!/bin/bash

ml MACS2/2.2.6-foss-2019b-Python-3.7.4
ml SAMtools/1.10-GCCcore-8.3.0
ml BEDTools/2.29.2-GCC-9.3.0
ml picard/2.21.6-Java-11


picardCMD="java -jar $EBROOTPICARD/picard.jar"
## ========== Stringent Setting =======
## call peaks by MACS2
$picardCMD SortSam I=$outPath/$treatFile/alignment/bowtie2_align.sorted.rmDup.sam O=$outPath/$treatFile/alignment/bowtie2_align.sorted.rmDup.sortName.sam SORT_ORDER=queryname
samtools view -bS $outPath/$treatFile/alignment/bowtie2_align.sorted.rmDup.sortName.sam >$outPath/$treatFile/alignment/bowtie2_align.sorted.rmDup.sortName.bam

$picardCMD SortSam I=$outPath/$controlFile/alignment/bowtie2_align.sorted.rmDup.sam O=$outPath/$controlFile/alignment/bowtie2_align.sorted.rmDup.sortName.sam SORT_ORDER=queryname
samtools view -bS $outPath/$controlFile/alignment/bowtie2_align.sorted.rmDup.sortName.sam >$outPath/$controlFile/alignment/bowtie2_align.sorted.rmDup.sortName.bam
## rm -rf $outPath/$controlFile/alignment/bowtie2_align.sorted.rmDup.sortName.sam

macs2 callpeak \
      -t $outPath/$treatFile/alignment/bowtie2_align.sorted.rmDup.sortName.bam \
      -c $outPath/$controlFile/alignment/bowtie2_align.sorted.rmDup.sortName.bam \
      -g hs -f BAMPE -n macs2_peak --outdir $resultsPath/stringent/MACS2 -q 0.1 -B --SPMR --keep-dup all 2>$resultsPath/stringent/MACS2/macs2Peak.txt


## =========== Relaxed Setting =========
## call peaks by SEACR

## prepare input bedgraph file for seacr %% according to seacr instructions
for file in $controlFile #$treatFile
do
    samtools view -bS $outPath/$file/alignment/bowtie2_align.sam >$outPath/$file/alignment/bowtie2_align.bam
    samtools view -bS $outPath/$file/alignment/bowtie2_align.sorted.dupMarked.sam >$outPath/$file/alignment/bowtie2_align.sorted.dupMarked.bam
    samtools view -bS $outPath/$file/alignment/bowtie2_align.sorted.rmDup.sam >$outPath/$file/alignment/bowtie2_align.sorted.rmDup.bam
    # rm -rf $outPath/$file/alignment/bowtie2_align.sam
    # rm -rf $outPath/$file/alignment/bowtie2_align.sorted.sam
    # rm -rf $outPath/$file/alignment/bowtie2_align.sorted.dupMarked.sam
    # rm -rf $outPath/$file/alignment/bowtie2_align.sorted.rmDup.sam

    bedtools bamtobed -i $outPath/$file/alignment/bowtie2_align.bam -bedpe >$outPath/$file/alignment/bowtie2_align.bed

    awk '$1==$4 && $6-$2 < 1000 {print $0}' $outPath/$file/alignment/bowtie2_align.bed >$outPath/$file/alignment/bowtie2_align.clean.bed
    
    cut -f 1,2,6 $outPath/$file/alignment/bowtie2_align.clean.bed | sort -k1,1 -k2,2n -k3,3n  >$outPath/$file/alignment/bowtie2_align.fragments.bed

    bedtools genomecov -bg -i $outPath/$file/alignment/bowtie2_align.fragments.bed -g $chromSize >$outPath/$file/alignment/bowtie2_align.fragments.bedgraph


    
    ## calculate fragment length
    awk '{print $3-$2}' $outPath/$file/alignment/bowtie2_align.fragments.bed | sort | uniq -c | awk -v OFS="\t" '{print $2, $1}' >$outPath/$file/alignment/bowtie2_align.fragLen.freq
    

done

bash $seacr $outPath/$treatFile/alignment/bowtie2_align.fragments.bedgraph \
     $outPath/$controlFile/alignment/bowtie2_align.fragments.bedgraph \
     norm stringent $resultsPath/relax/SEACR/seacr_control.peaks

bash $seacr $outPath/$treatFile/alignment/bowtie2_align.fragments.bedgraph \
     0.01 \
     non stringent $resultsPath/relax/SEACR/seacr_top0.01.peaks


## call peaks by MACS2 on sorted sam file
macs2 callpeak -t $outPath/$treatFile/alignment/bowtie2_align.bam \
      -c $outPath/$controlFile/alignment/bowtie2_align.bam \
      -g hs -f BAMPE -n macs2_peak_q0.1 --outdir $resultsPath/relax/MACS2 -q 0.1 --keep-dup all 2>$resultsPath/relax/MACS2/macs2Peak.txt

macs2 callpeak -t $outPath/$treatFile/alignment/bowtie2_align.bam \
      -c $outPath/$controlFile/alignment/bowtie2_align.bam \
      -g hs -f BAMPE -n macs2_peak_p0.001 --outdir $resultsPath/relax/MACS2 -p 0.001 --keep-dup all 2>$resultsPath/relax/MACS2/macs2Peak.txt
