#!/usr/bin/env bash
# -*- coding: utf-8 -*-
# an end to end scaling genome assembly for pacbio and oxford nanopore reads on cluster
# pacbio and oxford nanopore denovo assembly
def __init__(self, name, queue, threads, core, memory, user, change, mail,filename):
        self.name = name
        self.queue = queue
        self.threads = threads
        self.core = core
        self.memory = memory
        self.dir = os.getcwd()
        self.change = os.chroot(change)
        self.mail = mail
        self.user = user
        self.filename = filename
def read_file(self):
        self.filecontents = open(self.filename, 'rb')
        for i in self.filecontents.readlines():
            if i.startswith('#SBATCH'):
                raise FileExistsError('Already a configuration file')
            else:
                print(f'the_name_of_the_filename')
def getCluster(self):
        return print(f'#SBATCH -J self.name \
                     \n#SBATCH -p constraint="snb|hsw \
                     \n#SBATCH -p self.queue \
                     \n#SBATCH -n self.threads \
                     \n#SBATCH -c self.core \
                     \n#SBATCH --mem=self.memory \
                     \n#SBATCH --workdir = self.change \
                     \n#SBATCH --mail = self.mail \
                    \n#SBATCH --mail-type=END')
def writeCluster(self):
            self.filecontent = open(self.filename, 'rb')
            self.filecontents.write(self.getCluster())
            self.filecontents.write(self.getcwd())
            self.filecontents.close()
            print(f'the_cluster_file_for_the_configuration_run:{self.filecontents}')
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
lordfast index $reference
lordfast --search $reference --seq $files \
      --out $reference.$files.sam --threads $threads --numMap 10
samtools view -bS $reference.$files.sam $reference.$files.bam
bamtools convert -in $reference.$files.bam -format fasta
lordfast-correct -T $threads -S $files.stats.out \
     -i $files.mapping.bam.fasta -2 $files.entire.fasta \
     -k 21 -o $files.corrected.fasta -s 3
canu gridOptions="--time=24:00:00" -corMhapSensitivity=high \
    -p $species_prefix \
    -d $directory genomeSize=$genome_size -minReadLength=$readlength \
    -merylMemory=64 -gnuplotImageFormat=png \
    -ovsThreads=$threads \
    -ovbThreads=$threads \
    -pacbio-raw $files.corrected.fasta   

