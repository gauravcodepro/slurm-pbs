#!/usr/bin/env bash
# -*- coding:  utf-8 -*-
# Author: Gaurav Sablok
# date: 2023-08-09
# MIT License
# from bacterial genome sequencing to assembly
echo "hi there i am a denovo assembler script which will check your fastq files, 
clean then and assemble them coming from the illumina platform. I take MiSeq, NextSeq, 
NovaSeq reads, you can just clean, or clean and map or clean and assemble.As of the 
current it supports the cleaning through the fastp, mapping through the bowtie2 and 
assembly using the spades I will add more assembler to this interface and also multiple 
assembly configurations. Currently, it supports 45, 55, 65 kmers you can change the 
kmer according to the assembly preferences. Please kindly remove the echo from line 32,
44, 50, 51, 63, 69, 70, 76"
read -p "enter the name of the directory: " dirname
read -p "do you want to clean the reads:" cleanreads
read -p "enter the number of the threads:" thread
read -p "do you want to map the reads:" mapper
read -p "enter the name of the bowtie index:" index
read -p "do you want to assemble the reads:" assemble
if
    [[ -d "${dirname}" ]] &&
        [[ $cleanreads == "yes" ]]
    [[ $mapper == "no" ]] &&
        [[ $assemble == "no" ]] &&
        [[ $index == "no" ]]
then
    for i in ${dirname}/*.R1.fastq; do
        for j in ${dirname}/*.R2.fastq; do
            echo "Processing the cleaning of the reads file from the sequencing runs"
        done
        echo "fastp --in1 $i --out1 $i.cleaned.R1 --in2 $j --out2 $j.cleaned.R2 --threads $thread"
    done
fi
if [[ -d "${dirname}" ]] &&
    [[ $cleanreads == "yes" ]] &&
    [[ $mapper == "yes" ]] &&
    [[ $index == "$index" ]] &&
    [[ $assemble == "no" ]]; then
    for i in ${dirname}/*.R1.fastq; do
        for j in ${dirname}/*.R2.fastq; do
            echo "Processing the cleaning of the reads file from the sequencing runs"
        done
        echo "fastp --in1 $i --out1 $i.cleaned.R1 --in2 $j --out2 $j.cleaned.R2 --threads $thread"
    done
    for i in "${dirname}"/*.cleaned.R1; do
        for j in "${dirname}"/*cleaned.R2; do
            echo "Processing the mapping of the reads file from the sequencing runs"
        done
        echo "bowtie2-build $index $index"
        echo "bowtie2 -t -x $index -p $thread --very-sensitive-local -1 $i -2 $j -S $index.sam --no-unal --al-conc $index.aligned.fastq"
    done
fi
if [[ -d "${dirname}" ]] &&
    [[ $cleanreads == "yes" ]] &&
    [[ $mapper == "yes" ]] &&
    [[ $index == "$index" ]] &&
    [[ $assemble == "yes" ]]; then
    for i in ${dirname}/*.R1.fastq; do
        for j in ${dirname}/*.R2.fastq; do
            echo "Processing the cleaning of the reads file from the sequencing runs"
        done
        echo "fastp --in1 $i --out1 $i.cleaned.R1 --in2 $j --out2 $j.cleaned.R2 --threads $thread"
    done
    for i in "${dirname}"/*.cleaned.R1; do
        for j in "${dirname}"/*cleaned.R2; do
            echo "Processing the mapping of the reads file from the sequencing runs"
        done
        echo "bowtie2-build $index $index"
        echo "bowtie2 -t -x $index -p $thread --very-sensitive-local -1 $i -2 $j -S $index.sam --no-unal --al-conc $index.aligned.fastq"
    done
    for i in "${dirname}"/*.aligned.R1.fastq; do
        for j in "${dirname}"/*.aligned.R2.fastq; do
            echo "Processing genome assembly"
        done
        echo "spades.py -1 $i -2 $j --careful --threads $thread --tmp-dir $index_tmp -k 45,55,65,75 -o ${dirname}/$index_assembly"
    done
fi

```
shell_command
hi there i am a denovo assembler script which will check your fastq files, \
clean then and assemble them coming from the illumina platform. I take MiSeq, NextSeq,\
NovaSeq reads, you can just clean, or clean and map or clean and assemble.As of the \
current it supports the cleaning through the fastp, mapping through the bowtie2 and \
assembly using the spades I will add more assembler to this interface and also multiple \
assembly configurations. Currently, it supports 45, 55, 65 kmers you can change the \
kmer according to the assembly preferences.\
# if from cleaning only 
main ✗ $ sh genome_cleaner.sh
enter the name of the directory: /Users/gauravsablok/Desktop/GitHub/shell_builder \
do you want to clean the reads:yes \
enter the number of the threads:1 \
do you want to map the reads:no \
enter the name of the bowtie index:no \
do you want to assemble the reads:no \
Processing the cleaning of the reads file from the sequencing runs
fastp --in1 /Users/gauravsablok/Desktop/GitHub/shell_builder/text.txt \
        --out1 /Users/gauravsablok/Desktop/GitHub/shell_builder/text.txt.cleaned.R1 \
        --in2 /Users/gauravsablok/Desktop/GitHub/shell_builder/text1 \
        --out2 /Users/gauravsablok/Desktop/GitHub/shell_builder/text.txt.cleaned.R2 --threads 1 \
# if cleaning to remapping files 
enter the name of the directory: /Users/gauravsablok/Desktop/GitHub/shell_builder \
do you want to clean the reads:yes \
enter the number of the threads:1 \
do you want to map the reads:yes \
enter the name of the bowtie index:gaurav \
do you want to assemble the reads:no \
Processing the cleaning of the reads file from the sequencing runs \
fastp --in1 /Users/gauravsablok/Desktop/GitHub/shell_builder/text.txt \
        --out1 /Users/gauravsablok/Desktop/GitHub/shell_builder/text.txt.cleaned.R1 \
        --in2 /Users/gauravsablok/Desktop/GitHub/shell_builder/text1 \
        --out2 /Users/gauravsablok/Desktop/GitHub/shell_builder/text.txt.cleaned.R2 --threads 1 \
bowtie2-build gaurav gaurav \
bowtie2 -t -x gaurav -p 1 --very-sensitive-local -1 \
       /Users/gauravsablok/Desktop/GitHub/shell_builder/text.txt.cleaned.R1 \
         -2 /Users/gauravsablok/Desktop/GitHub/shell_builder/text3.txt.cleaned.R2 \
                 -S gaurav.sam --no-unal --al-conc gaurav.aligned.fastq \                          
# if from cleaning to assembly configuration \
enter the name of the directory: /Users/gauravsablok/Desktop/GitHub/shell_builder \
do you want to clean the reads:yes \
enter the number of the threads:1 \
do you want to map the reads:yes \
enter the name of the bowtie index:gaurav \
do you want to assemble the reads:yes \
Processing the cleaning of the reads file from the sequencing runs \
fastp --in1 /Users/gauravsablok/Desktop/GitHub/shell_builder/text.txt \
        --out1 /Users/gauravsablok/Desktop/GitHub/shell_builder/text.txt.cleaned.R1 \
        --in2 /Users/gauravsablok/Desktop/GitHub/shell_builder/text1 \
        --out2 /Users/gauravsablok/Desktop/GitHub/shell_builder/text.txt.cleaned.R2 --threads 1 \
bowtie2-build gaurav gaurav \
bowtie2 -t -x gaurav -p 1 --very-sensitive-local -1 \
       /Users/gauravsablok/Desktop/GitHub/shell_builder/text.txt.cleaned.R1 \
         -2 /Users/gauravsablok/Desktop/GitHub/shell_builder/text3.txt.cleaned.R2 \
                 -S gaurav.sam --no-unal --al-conc gaurav.aligned.fastq \
spades.py -1 /Users/gauravsablok/Desktop/GitHub/shell_builder/text.aligned.R1.fastq \
       -2 /Users/gauravsablok/Desktop/GitHub/shell_builder/*.aligned.R2.fastq \
                 --careful --threads 1 --tmp-dir  -k 45,55,65,75 \
                                 -o /Users/gauravsablok/Desktop/GitHub/shell_builder/ \