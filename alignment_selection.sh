PS3 = "Genomic alignment of the PacBio and the Illumina reads (#)? "; \ 
select choice in blasr lra_sam lra_paf lra_bed S-conLSH_paf S-conLSH_sam bowtie star; do
if [[ $choice == bowtie ]]
   then echo "bowtie2-build $genome $genome \
        paste Illumina.R1.txt Illumina.R2.txt \
                        | while read col1 col2; \
                     do echo bowtie2 -t -x $genome \
                  -p $thread --very-sensitive-local \
                -1 ${col1} -2 ${col2} -S $genome.sam \
                    --no-unal --al-conc $genome.aligned.fastq"; break; fi
          
if [[ $choice == star ]]
    then echo "STAR --runThreadN 6 \
                     --runMode genomeGenerate \
                     --genomeDir $genome_dir \
                     --genomeFastaFiles $reference \
                     --sjdbGTFfile $gtf \
                     --sjdbOverhang 99 \
               STAR --genomeDir $dir \
                      --runThreadN $threads \
                      --readFilesIn $IlluminaR1.txt $IlluminaR2.txt \
                        --outFileNamePrefix $genome \
                          --outSAMtype BAM SortedByCoordinate \
                            --outSAMunmapped Within \
                              --outSAMattributes Standard"; break; fi                 
if [[ $choice == blasr ]]
   then echo "blasr $reads $genome --bam \
         $genome.bam --unaligned $genome.unaligned.fasta \
              --nproc $threads"; break; fi
if [[ $choice == lra_sam ]]
   then echo " lra index -CCS/CLR/ONT/CONTIG $reference \
         lra align --CCS/CLR/ONT/CONTIG $reference $reads \
             -t $threads -p -s > $reference.sam"; break; fi
if [[ $choice == lra_paf ]]
   then echo " lra index -CCS/CLR/ONT/CONTIG $reference \
         lra align --CCS/CLR/ONT/CONTIG $reference $reads \
             -t $threads -p -p > $reference.paf"; break; fi
if [[ $choice == lra_bed ]]
   then echo "lra index -CCS/CLR/ONT/CONTIG $reference \
         lra align --CCS/CLR/ONT/CONTIG $reference $reads \
             -t $threads -p -b > $reference.bed"; break; fi
if [[ $choice == S-conLSH_paf ]]
    then echo " S-conLSH $reference $reads > $reference.paf"; break; fi
if [[ $choice == S-conLSH_sam ]]
     then echo " S-conLSH $reference $reads --align 1 > $reference.sam"; break; fi
     done
