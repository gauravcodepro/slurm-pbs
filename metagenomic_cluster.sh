# a long read metagenomic assembly analysis for the cluster
#plus a streamlit app for the visualization of the frequencies
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
long_read = 'PATH_DIRECTORY'
short_read = 'PATH_DIRECTORY'
long_read_final = 'PATH_DIRECTORY'
threads = 'THREAD_COUNT' or $SLURM_TASK or $PBS
long_clean = 'PATH_DIRECTORY'
threshold = 'VALUE'
barcodes = 'barcodes'
databasedir = 'PATH_DIRECTORY'
for f in long_read/*.gz; do echo md5sum $f; done > ./long_read_file.txt
for f in short_read/*.R1.fastq.gz; do echo \
      gunzip $f; done >> unzipping.short.read.R1.sh 
for f in short_read/*.R2.fastq.gz; do echo \
      gunzip $f; done >> unzipping.short.read.R2.sh 
sh unzipping.short.read.R1.sh && echo "unzipping_done" done;
sh unzipping.short.read.R2.sh && echo "unzipping_done" done;
for f in long_read/*.gz; do echo gunzip $f; \
       done >> long_read_all.sh
sh long_read_all.sh
cat *.fasta >> long_read_all.fasta
porechop -i long_read_all.fasta -o long_read_all.fasta 
            --threads $threads 
porechop -i long_read_all.fasta -b $long_clean 
              --barcode_threshold $value --$barcodes
krakenuniq-download --db $databasedir --threads $dir \
                  --dust refseq/bacteria refseq/archaea
krakenuniq-build --db $databasedir  
krakenuniq --db $databasedir --threads $threads \
      --report-file report.tsv > classification.tsv
import streamlit as st
import matplotlib.pyplot as plt
from wordcloud import WordCloud
st.subheader("Text of classification")
classify_text=st.text_area("classify text")
if st.button("Generate frequency chart"):
	wc=WordCloud().generate(classify_text)
	plt.imshow(wc)
	plt.axis('off')
	st.pyplot(plt)



