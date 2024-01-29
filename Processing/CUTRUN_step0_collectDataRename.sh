#!/bin/bash

inPath="CART_CUTRUN_Project"

fileNameL=($(awk -F, '{ if (NR>1) {print $1}}' ${inPath}/meta/metadata_hd7.csv))
folderNameL=($(awk -F,  '{ if (NR > 1) {print $2}}' ${inPath}/meta/metadata_hd7.csv))
rawPathL=($(awk -F, '{ if (NR >1) {print $14}}' ${inPath}/meta/metadata_hd7.csv))

for i in $(seq 0 ${#fileNameL[@]})
do
    fileName=${fileNameL[$i]}
    folderName=${folderNameL[$i]}
    rawPath=${rawPathL[$i]}

    ## check folderName if it start with "Sample"
    head="$(cut -d'_' -f1 <<<"$folderName")"
    if [ "$head" != "Sample" ]; then
	folderName="Sample_${folderName}"
    fi

    count1=`ls $rawPath/$folderName/*R1* | wc -l`
    count2=`ls $rawPath/$folderName/*R2* | wc -l`
    echo "${folderName}: ${count1}-${count2}"
    if [ "$count1" != "$count2" ] || [ "$count1" -lt "1" ] || [ "$count2" -lt "1" ]; then
	echo "!!!!!!!!Warning ${folderName}: ${count1}-${count2}"
    fi
    
    mkdir -p $inPath/data/CUTANDRUN/$fileName
    cat $rawPath/$folderName/*R1*.fastq.gz >$inPath/data/CUTANDRUN/$fileName/R1.fastq.gz
    cat $rawPath/$folderName/*R2*.fastq.gz >$inPath/data/CUTANDRUN/$fileName/R2.fastq.gz
    ## mv $inPath/data/CUTANDRUN/$folderName $inPath/data/CUTANDRUN/$fileName
    
done
