#!/bin/bash

mkdir -p $outPath/$name
scp $inPath/${name}_${end}.fastq.gz $outPath/$name/
gunzip $outPath/${name}/${name}_${end}.fastq.gz

/fh/fast/gottardo_r/yezheng_working/Software/FastQC/fastqc -o $outPath/$name -f fastq $outPath/$name/${name}_${end}.fastq
rm -rf $outPath/$name/${name}_${end}.fastq
mv $outPath/$name/${name}_${end}_fastqc.html $outPath
