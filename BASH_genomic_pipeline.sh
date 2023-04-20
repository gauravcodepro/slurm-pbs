# End_to_end_genomics_pipeline
#assembly version
#/usr/bin/python 
#!/bin/bash -l
# write the cluster configuration
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
#folder_set_up
mkdir files_sequencing_1 files_sequencing_2 tag
cd files_sequencing_1 && wget 'link'
cd files_sequencing_2 && wget 'link'
for file in files_sequencing/*; do echo $file >> file1_tags.txt
for file in file_sequencing/*; do echo $file >> file2_tags.txt
#generating the length tags
for file in files_sequencing/*.fastq; do echo $file && wc -l $file; done >> file_size_tags_1.txt
for file in files_sequencing/*.fastq; do echo $file && wc -l $file; done >> file_size_tags_2.txt
#generating the length count tags
#python way 
line_tags_R1=[]
with open(file_size_tags_1.txt, 'rb') as fname:
    for line in fname.readline():
        if line.startswith(r'[0-9]'):
            line_tags.append(line)
            fname.write(line_tags)
            fname.close()
line_tags_R2=[]
with open(file_size_tags_2.txt, 'rb') as fname:
    for line in fname.readline():
        if line.startswith(r'[0-9]'):
            line_tags.append(line)
            fname.write(line_tags)
            fname.close()
print(map(lambda n: ''.join(n),line_tags_R1))
print(map(lambda n: ''.join(n),line_tags_R2))
#grep way 
for file in file_sequencing/*fastq; do echo wc -l $file; done >> file_size.txt
#running the alignments for the genome_filtering
#checking if bowtie2 is installed or not
echo $PATH | grep bowtie2
#adding the aligner to the PATH environment variable
PATH='path_to_aligner'
export PATH=$PATH:$PATH
#listing all the files 
files1 = file_size_tags1.txt
files2 = file_size_tags2.txt 
filecontent = file3_tags.txt # species names only
thread = indicate the thread
paste file1_tags.txt file2_tags.txt file3_tags.txt \
    | while read col1 col2 col3; \
        do echo bowtie2 -t -x "${col3}" \
            -p $thread --very-sensitive-local \
                -1 ${col1} -2 ${col2} -S "${col3}".sam \
                    --no-unal --al-conc ${col3}.aligned.fastq; done 
                    echo 'mapping_of_reads_finished'
# separating the samfiles
# and converting them to bamfiles
mkdir samfiles
for f in *.sam; do samtools view -bam $f -sam $f; done
mkdir right_aligned
mkdir left_aligned
mv *.aligned.R1.fastq ./right_aligned 
mv *.aligned.R2.fastq ./left_aligned
for file in *.aligned.R1.fastq; do echo $f; done > R1_files.txt
for file in *.aligned.R2.fastq; do echo $f; done > R2_files.txt
# counting the aligned reads insert size
sudo apt-get install bamtools 
echo $PATH
paste R1_files.txt R2_files.txt \
  | while read R1 R2; \
    do echo bamtools -1 ${R1} -2 ${R2} --insert; done
           echo 'insert_size_calculation_finished'
# install the assembler of your choice
line_tags_R1 = []
with open(file1_assembly.txt, 'rb+') as fname:
    for line in fname.readline():
            line_tags.append(line)
            fname.write(line_tags_R1)
            fname.close()
line_tags_R2 = []
with open(file2_assembly.txt, 'rb+') as fname:
    for line in fname.readline():
            line_tags.append(line)
            fname.write(line_tags_R2)
            fname.close()
# python way
import subprocess
for i in range(len(line_tags_R1):
               subprocess.run('assembler', \
                        '--assembly.fasta', \
                        '--R1', \
                        'line_tags_R1[i]', \
                        '--R2', 'line_tags_R2[j]', \
                        text=True, \
                        capture_output=True)
# bash way 
paste line_tags_R1 line_tags_R2 \
  | while read R1 R2;\
    do echo assembler\
        --assembly.fasta\ 
        --R1 {$R1}\
        --R2 {$R2}\
        done
