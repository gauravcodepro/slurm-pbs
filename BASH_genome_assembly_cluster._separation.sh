# a server cluster file for read-alignment 
# your BAM files and computing stats for
# generating all the results from the sequencing runs
#after the cluster separation
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
fastq_files1 = file_1.txt
fastq_files2 = file_2.txt 
file_information = information.txt
file_names=file_names.txt 
file = tags.txt # species names only
thread = indicate_the_thread
mapping_file = file.fasta
# processing the single file for the alignment 
bowtie2-build $mapping_file $mapping_file
paste file1.txt file2.txt \
    | while read col1 col2 col3; \
        do echo bowtie2 -t -x $mapping_file \
            -p $thread --very-sensitive-local \
                -1 ${col1} -2 ${col2} -S "${mapping_file}".sam \
                    --no-unal --al-conc ${mapping_file}.aligned.fastq; done 
                    echo 'mapping_of_reads_finished'
#processing an array of the files for the alignment
for f in */*.fasta; do echo $f; done >> file_information.text
# spliting the names at the . as the .fasta is not needed
cat file_information.txt | cut -f -1 -d "." > filenames.txt
cat filename.txt | while read line; do echo bowtie2-build \
    $line $line; done
    echo "bowtie2_indexes_are_build"
for f in *.bowtie; do echo $f; done >> bowtie.indices.txt 
paste bowtie.indices.txt  file1.txt file2.txt \
    | while read col1 col2 col3; \
        do echo bowtie2 -t -x ${col1} \
            -p $thread --very-sensitive-local \
                -1 ${col2} -2 ${col3} -S "${col1}".sam \
                    --no-unal --al-conc ${col1}.aligned.fastq; done 
                    echo 'mapping_of_reads_finished'
# generating the alignment and the insert stats 
for f in *.sam; do echo samtools -in $f -out $f.bam; done
for f in *.bam; do echo $f; bamtools stats -in \
             $f --insert >> $f.insert; done 
for f in *.bam; do echo $f; bamtools coverage \
            -in $f -out $f.coverage.txt; done
for f in *.bam; do echo $f; bamtools count -in $f \
     > $f.count.read.alignment.txt; done
# you can modify this according to your choice 
touch insert_coverage_alignment.txt | paste $f.insert.txt \
    $f.coverage.txt $f.count.read.alignment.txt 
