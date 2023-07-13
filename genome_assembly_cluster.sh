# running an genomic assembly on the cluster, 
#           remapping the assembly,filtering the assembly 
#                and genome assessment
# i am using the glimmer for the predictions
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
file = tags.txt # species names only
thread = indicate the thread
filter_files = file.fasta
e_value = file.e_value
e_type = e.type
# preprocessing the file for the alignment 
bowtie2-build $filter_fasta $filter_fasta
paste tags.txt file1.txt file2.txt \
    | while read col1 col2 col3; \
        do echo bowtie2-build ${tags}
# aligning the reads for the alignment file
paste tags.txt file1.txt file2.txt \
    | while read col1 col2 col3; \
        do echo bowtie2 -t -x "${col1}" \
            -p $thread --very-sensitive-local \
                -1 ${col2} -2 ${col3} -S "${col1}".sam \
                    --no-unal --al-conc ${col1}.aligned.fastq; done 
                    echo 'mapping_of_reads_finished'
# making the assembly 
paste tags.txt file1.txt file2.txt \
    | while read col1 col2 col3; \
        do echo w2rap-contigger -t $thread \
          -m 1000 -r ${col2},${col3} \
          -o ${col1}_output_dir_w2rap \
          -p ${col1}_kmer_prefix \
          -tmp_dir ${col1}_tmp_dir \
          -s 2; done 
          echo "w2rap_genome_assembly_finished"
for dir in *; do echo \
       $(pwd)/$dir; done >> listing_directories_assembly_path.txt
# making a common directory for all the assemblies
mkdir common_dir
for f in ${col3}_*_*/*.fasta; \
   do echo $f; mv $f ./common_dir; done
for f in common_dir; \
   do echo $f; grep ">" -c $f > $f.assembly.count.txt; done
   echo "Finished_processing_the_assembly_count"
for f in *.count.txt;
    do echo $f; cut -f 1 -d "\t" $f; > ${col1}.count.txt; done
    echo "assemnled_gene_count_are_present_in_the_count_file"
# filtering the genes
assembly_file.txt = mv common_dir; for f in $(pwd)/*.fasta; \
            do echo $f; done >> assmebly.file.txt
paste assembly_file.txt \
    | while read col1; \
        do echo blastn -db database -query ${col1} \
        -outfmt "6 qseqid staxids bitscore std" \
        -max_target_seqs 10 -max_hsps 1 -evalue 1e-25 -num_threads 75 -out \
        ${col1}.blast.out; done
        echo "BLAST_assembly_finished"
# gene completeness
for f in common_dir; \
    do echo $(pwd)$f; done >> file_assembled.txt
cat file_assembled_txt | \
  while readline do echo quast $line \
      -r $reference_genome_name \
      -g $reference_gene.gff3 \
      --fragmented \
      -o $line.gene.txt \
      -m 500 \
      -t $thread \
      -- circos \
      --glimmer \
      done; >> genome_completeness.sh  > /dev/null
# run_BUSCO 
for f in common_dir; \
    do echo $(pwd)$f; done >> file_assembled.txt
cat file_assembled_txt | \
  while readline do echo run_BUSCO.py -i $line \
      -c $threads \
      -o $line.busco.out \
      -e $e_value \
      -m $e.type \
      -l -sp \
      done; >> genome_completeness.sh  > /dev/null
