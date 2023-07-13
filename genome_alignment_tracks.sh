# a streamline analysis workflow for generating the chromosome 
# or the contigs tracts for the metgenomics with the variants
#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#SBATCH -J self.name 
#SBATCH -p constraint="snb|hsw 
#SBATCH -p self.queue 
#SBATCH -n self.threads 
#SBATCH -c self.core 
#SBATCH --mem=self.memory 
#SBATCH --workdir = self.change 
#SBATCH --mail = self.mail 
#SBATCH --mail-type=END'
files = PACBIO.fasta
directory = 'PATH'
genome_size = 'genome_size'
species_prefix = 'prefix'
readlength = 'read_length'
threads = 'num_of_threads'
canu gridOptions="--time=24:00:00" -corMhapSensitivity=high \
    -p $species_prefix \
    -d $directory genomeSize=$genome_size -minReadLength=$readlength \
    -merylMemory=64 -gnuplotImageFormat=png \
    -ovsThreads=$threads \
    -ovbThreads=$threads \
    -pacbio-raw $files
#reference based assembly
reference = reference.fasta
mecat2ref -d $files -r $reference -o \
    $files.reference.sam \
        -w $files_intermediate \
        -t 24 -n 20 -n 50 -m 2 -x 0 
#error correction and pacbio and oxford nanopore assembly
mkdir fastq_files
mv *.R1.fastq ./fastq_files 
mv *.R2.fastq ./fastq_files
for f in $(pwd)/fastq_files/*.R1.fastq; \
       do echo $f; done > fastq.R1.txt
for f in $(pwd)/fastq_files/*.R1.fastq; \
       do echo $f; done > fastq.R1.txt
bowtie2-build $species_prefix $species_prefix
paste fastq.R1.txt fastq.R2.txt \
    | while read col2 col3; \
        do echo bowtie2 -t -x $species_prefix \
            -p $thread --very-sensitive-local \
                -1 ${col1} -2 ${col2} -S $species_prefix.sam \
                    --no-unal --al-conc $species_prefix.aligned.fastq; done
# calling variants and building the variants tracks
samtools view -bS $species_prefix.sam -o $species_prefix.bam 
samtools sort $species_prefx.bam $species_prefix.sorted.bam
bamtools -index $species_prefix.bam
samtools index faidx $species_prefix.fasta
samtools mpileup -g -f $species_prefix.fasta > $species_prefix.bcf
bcftools view -bvcg $species_prefix.bcf > $species_prefix_variant.bcf
bcftools view $species_prefix_variant.bcf \
       | vcfutils.pl varFilter -d 10 -d 30 > $species_prefix_variant.vcf
tar zcvf $species_prefix_variant.vcf
# annotating the variants with the read depth 
gatk AnnotateVcfWithBamDepth --input $species_prefix.sorted.bam \
       --reference $species_prefix.fasta \
       --output $species_prefix_filter_variants.vcf \
       --variant $species_prefix_variant.vcf
tar zcvf $species_prefix_filter_variants.vcf
# writing the alignment track for the jbrowse
{
  "trackId": "long_read_alignment_track",
  "name": "organelle_genome_alignment",
  "assemblyNames": ["$species_prefix"],
  "type": "AlignmentsTrack",
  "adapter": {
    "type": "BamAdapter",
    "bamLocation": {
      "uri": "http://127.0.0.1/$species_prefix.bam"
    },
    "index": {
      "location": {
        "uri": "http://127.0.0.1/$species_prefix.bam.bai"
      }
    }
  }
}
# importing jbrowse
from jbrowse_jupyter import launch, create
browser = create ("CGV")
browser.set_assembly("$species_prefix.fa.gz", aliases=["$species_prefix.fa.gz"])
browser.add_track("$species_prefix_variant.vcf.tar.gz", trackID)
browser.add_track("$species_prefix_filter_variants.vcf.tar.gz", trackID)
browser.set_default_session(["trackID"])
cgv_conf = browser.get_config()
launch(cgv_conf,dash_comp="CGV",height=$height, port=8081)