# a genomic class for the genome size estimation from the multi sequencing reads
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
threads = NUM_THREADS
overlap = SET_THE_OVERLAP
memory = SET_MEMORY

mkdir counts
mv counts

pip install biopython=version
export PATH=/PATH/TO/JELLYFISH

paste tags.txt file1.txt file2.txt \
    | while read col1 col2 col3; \
        do echo jellyfish count -m $overlap \
        -s $memory -t $threads -o ${col1}.jf
        -C ${col2} ${col3}; done 
        echo "jellyfish_count_finished"

for f in $(pwd)/*.jf; do echo $f; done >> count_file.txt

paste count_file.txt tags.txt \
    | while read col1 col2; \
        do echo jellyfish histo -o ${col2}.histo \
        -t $threads -l 1 \
        -h 100000000 -i 1 \
        ${col1}; done 
        echo "histogram_files_prepared"

paste tags.txt file1.txt file2.txt \
    | while read col1 col2 col3; \
        do echo draw_kmer_distribution_from_fastq.py \
        -m 23 -t $threads -b -s $memory -e png \
        -o ${col1} -i ${col2},${col3} \
        -w 10 -g 100 -d \
        echo "cluster_genome_finished"
