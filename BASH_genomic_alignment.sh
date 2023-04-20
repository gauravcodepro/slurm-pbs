# genomic pipeline for genomic parsing
# on clusters 

# Author: Firstname Lastname
#
alignment_tools = "muscle,\
fasttools,\
samtools,\
clipkit,\
makealignment\"
# tasks
# download all the alignment tools
# set a common folder 
# generate alignments
# common matrix
# run phylogeny
# make genomics R based plots 



mkdir fasta_files = "$pwd"
mkdir download_alignment = "$pwd/download_alignment"
for i in $alignment_tools; do
    echo "Downloading $i"
    wget -O $i.tar.gz 



#!/usr/bin/python
# -*- coding: utf-8 -*-
# CONDA ENV SOFTWARE: python (3.7), samtools (v1.9), gatk4 (v4.1.4.0)

import os, sys, subprocess, shutil, argparse
from core import FindSupplementaryFile, CAPTURE

pipe_script = os.path.realpath(__file__) # extract script path
script_name, *arguments = sys.argv # extract command line arguments

split_script = f'{os.path.dirname(pipe_script)}/fasta_split.sh' # specify module script
assert(os.path.exists(split_script)), f'Problem finding the module "{split_script}".'# check module exists

parser = argparse.ArgumentParser(description='Merge contigs', prog=script_name, usage='%(prog)s [options]', epilog='see the readme file for further details.')

inputs = parser.add_argument_group(description='user inputs:') # user inputs
inputs.add_argument('-fasta', metavar='</path/to/directory>', type=str, required=True, help='specify fasta path')
inputs.add_argument('-contigs', metavar='<number>', type=int, default=500, help='specify final contig number')
inputs.add_argument('-delimiter', metavar='<character>', type=str, default='N', help='specify delimiter symbol')
inputs.add_argument('-length', metavar='<number>', type=int, default=500, help='specify delimiter number')

fasta_dir, total_merged_contigs, delimiter_character, delimiter_length = vars(parser.parse_args()).values() # define user inputs

fasta_dir = f'/{fasta_dir.strip("/")}' # ensure correct path format
assert(os.path.exists(fasta_dir)), f'Problem finding the reference directory "{fasta_dir}".'# check path exists

# file info

original_fasta_file = FindSupplementaryFile(fasta_dir, '.fasta') # find reference fasta

original_fasta_basename = os.path.basename(original_fasta_file) # extract reference basename

split_name = 'split'  # specify split fasta subdir name

split_dir = f'{fasta_dir}/{split_name}' # specify split fasta subdirectory path

file_names = [ # specify output file names
    'smallest_contigs_combined'
    ,f'sorted_{original_fasta_basename}'
    ,f'merged_{original_fasta_basename}' ]

output_files = [ f'{fasta_dir}/{name}' for name in file_names ] # specify output file paths

merged_contig_name, sorted_fasta_name, rebuilt_fasta_name = file_names # extract individual file names

merged_small_contig_file, sorted_fasta_file, rebuilt_fasta_file = output_files # extract individual file paths

large_contig_subdir, small_contig_subdir = size_subdirs = [ f'{split_dir}/{label}' for label in ['large','small'] ] # specify contig size subdirectory paths

merge_cmd = ( # specify command to merge smallest contigs
f'''awk -v delimiter=`printf \'{delimiter_character}%.0s\' {{1..{delimiter_length}}}` \\
\'BEGIN {{ print \">{merged_contig_name}\" }} \\
FNR>1 {{ printf \"%s%s\", seq, $0 ; seq=delimiter }} \\
END {{ print \"\" }} \' \\
{small_contig_subdir}/*.fasta > {merged_small_contig_file}'''
            )
recombine_sorted_cmd = (f'cat {large_contig_subdir}/*.fasta {small_contig_subdir}/*.fasta > {sorted_fasta_file}') # specify command to recombine fasta with contigs sorted by size

recombine_merged_cmd = (f'cat {large_contig_subdir}/*.fasta {merged_small_contig_file} > {rebuilt_fasta_file}') # specify command to recombine fasta with largest contigs seperate & smallest contigs merged

total_contigs = int(CAPTURE(f'grep -c ">" {original_fasta_file}')) # calculate total number of contigs

print(f'\n   ORIGINAL CONTIGS: {total_contigs:,}')

print(f'\n   MERGED CONTIGS: {total_merged_contigs:,}')

if total_contigs <= total_merged_contigs: print(f'\n   REFERENCE ALREADY HAS THE REQUIRED NUMBER OF CONTIGS\n')
else:

    
    os.makedirs(split_dir, exist_ok=True) # create subdirectory as required

    print('\n   INDEXING REFERENCE...')
    subprocess.call(f'samtools faidx {original_fasta_file}', shell=True) # index reference fasta
    print('\tCOMPLETE')

    print('\n   SPLITTING REFERENCE...')
    subprocess.run(f'bash {split_script} {original_fasta_file} {split_dir}', shell=True, stdout=subprocess.DEVNULL) # split reference fasta via module
    print('\tCOMPLETE')


    contigs_to_retain = total_merged_contigs-1 # specify number of original contigs to retain i.e. not to be merged
    
    [ os.makedirs(subdir, exist_ok=True) for subdir in size_subdirs ] # create subdirectories as required

    fai_file = FindSupplementaryFile(fasta_dir, '.fai') # find reference index

    faidx_info = [ line.strip(' \n').split('\t') for line in open(fai_file, 'r').readlines() ] # extract index info

    contig_info = [ ( name, int(length) ) for name, length, *_ in faidx_info ] # extract contig names & lengths

    sorted_contig_info = sorted( contig_info, key=lambda info: info[1], reverse=True ) # sort contig info by length; largest -> smallest

    contig_names, *_ = zip(*sorted_contig_info) # discard contig lengths

    contig_split_by_size = [ # split contig info by size
        ( contig_names[ :contigs_to_retain ], large_contig_subdir ) # largest contig info
        ,( contig_names[ contigs_to_retain: ], small_contig_subdir ) # smallest contig info
                ]

    padding = len(str(total_contigs)) # calculate leading zeroes for index label

    contig_files = { contents.name.replace('.fasta',''):contents.path for contents in os.scandir(split_dir) if contents.name.endswith('.fasta') } # list contig files

    print('\n   ORGANISING CONTIGS...')
    for contigs, subdir in contig_split_by_size:

        for name in contigs:

            contig_file = contig_files[name] # extract relevant contig file

            new_name = str( contig_names.index(name)+1 ).zfill( padding ) # specify contig index label

            renamed_contig_file = f'{subdir}/{new_name}.fasta' # specify renamed contig file
            
            subprocess.call(f'mv {contig_file} {renamed_contig_file}', shell=True) # move contig to relevant subdirectory & rename with index label

    subprocess.call( merge_cmd, shell=True ) # merge smallest contigs together
    print('\tCOMPLETE')

    print('\n   REBUILDING REFERENCE...')
    subprocess.call( recombine_sorted_cmd, shell=True ) # recombine sorted reference fasta

    subprocess.call( recombine_merged_cmd, shell=True ) # recombine merged reference fasta
    print('\tCOMPLETE')

    print(f'\n   ALL STAGES COMPLETE\n')

paste R1_files.txt R2_files.txt names.txt | while read col1 col2 col3 ; do echo megahit -1 ${col1} -2 ${col2} --k-list 79,99,119 --memory 0.7 --num-cpu-threads 70 --out-dir "${col3}"_LSU_megahit_assembly --tmp-dir /mnt/long_term_storage/fungal_mitochondrial_assemblies/gene_wise_mapping/LSU_SSU/LSU_B_plate/tmp; done
canu -corMhapSensitivity=normal -corMinCoverage=3 -corErrorRate=0.013 -corOutCoverage=100 -p femsjoenia -d /mnt/gasablok2/femsjoenia_additional_assembly/femsjoenia_canu_assembly genomeSize=60m -minReadLength=500 -merylMemory=256 -gnuplotImageFormat=png -batMemory=300 -cnsMemory=300 -pacbio-raw /mnt/gasablok2/femsjoenia_additional_assembly/All.Femzonia.Pacbio.fasta
flye --pacbio-raw /mnt/gasablok2/blasia_assembly_2020/ALL.BLASIA.PACBIO_CLEANED_2020.fasta --genome-size 540000000 --threads 40 --out-dir /mnt/gasablok2/blasia_assembly_2020/blasia_assembly_2020_cleaned_FLYE --min-overlap 1000
bowtie2-build 110.Bridgeoporus_assembly_spades.scaffolds.fasta 110.Bridgeoporus_assembly_spades
bowtie2 -t -x 110.Bridgeoporus_assembly_spades -p 60 --very-sensitive-local -1 /mnt/gasablok2/fungal_ITS_assembly/fungal_raw_reads/110_Bridgeoporus_nobilissimus_S109_L001_R1_001.fastq -2 4/mnt/gasablok2/fungal_ITS_assembly/fungal_raw_reads/110_Bridgeoporus_nobilissimus_S109_L001_R2_001.fastq -S 110.Bridgeoporus_assembly_spades.sam --no-unal --al-conc 110.Bridgeoporus_assembly_spades.aligned.fastq
for f in *.fasta; do echo python funannotate.py predict -i $f --force -s "$f"_annotations --name "$f"_AN -o "$f"_ANF --protein_evidence All.proteins_agaromycotina.fasta --augustus_species yarrowia_lipolytica --cpus 40; done
orthofinder -t 74 -a 60 -M dendroblast -S diamond  -M msa -A mafft -T fasttree -p ./ -n auricularales_index_with_annotation -f /mnt/gasablok2/protein_directory
canu -corMhapSensitivity=normal -corMinCoverage=3 -corErrorRate=0.013 -corOutCoverage=100 -p simplicillium -d /mnt/gasablok2/femsjoenia_additional_assembly/simplicillium_lamellicola_updated_genome genomeSize=40m -minReadLength=500 -merylMemory=256 -gnuplotImageFormat=png -batMemory=300 -cnsMemory=300 -pacbio-raw /mnt/gasablok2/femsjoenia_additional_assembly/simplicillium_lamellicola_updated_genome/All.simplicillium.fasta
paste R1_files.txt R2_files.txt names.txt | while read col1 col2 col3 ; do echo bowtie2 -t -x Final_TEF1_Phylogeny_headers -p 60 --very-sensitive-local -1 ${col1} -2 ${col2} -S ${col3}.sam --no-unal --al-conc ${col3}.aligned.TEF1.fastq; done
paste R1_files.txt R2_files.txt names.txt | while read col1 col2 col3 ; do echo bowtie2 -t -x Final_RPB2_Phylogeny_headers -p 60 --very-sensitive-local -1 ${col1} -2 ${col2} -S ${col3}.sam --no-unal --al-conc ${col3}.aligned.RPB2.fastq; done
paste R1_files.txt R2_files.txt names.txt | while read col1 col2 col3 ; do echo bowtie2 -t -x Final_RPB1_Phylogeny_headers -p 60 --very-sensitive-local -1 ${col1} -2 ${col2} -S ${col3}.sam --no-unal --al-conc ${col3}.aligned.RPB1.fastq; done
megahit -1 /mnt/long_term_storage/chaetospiridium_updated_assembly_all_pacbio_all_illumina/Chaeto_ALL.R1.fastq -2 /mnt/long_term_storage/chaetospiridium_updated_assembly_all_pacbio_all_illumina/Chaeto_ALL.R2.fastq --k-list 79,99,119 --memory 0.7 --num-cpu-threads 70 --out-dir /mnt/long_term_storage/chaetospiridium_updated_assembly_all_pacbio_all_illumina/chaetospidium_megahit --min-contig-len 500 --tmp-dir /mnt/long_term_storage/chaetospiridium_updated_assembly_all_pacbio_all_illumina/tmp_megahit
paste list.R1.txt list.R2.txt names.txt | while read col1 col2 col3 ; do echo megahit -1 ${col1} -2 ${col2} --k-list 79,99,119 --memory 0.7 --num-cpu-threads 70 --out-dir "${col3}"_assembly --tmp-dir /mnt/long_term_storage/fungal_previous_sequencing_runs_MEGAHIT/tmp_megahit
./orthofinder -t 74 -a 60 -M dendroblast -S diamond  -M msa -A mafft -T fasttree -p ./ -n cryptomonas_phylogeny -f /mnt/gasablok2/cryptomonas_mallomonas_phylogenetic_reconstruction/genmark_based_phylogeny
orthofinder -t 74 -a 60 -M dendroblast -S diamond  -M msa -A mafft -T fasttree -p ./ -n host_parasite_reslationship -f /mnt/storage/fungal_sequencing/fungal_genome_annotations/funcannotation/funannotate/host_parasite_relationship_orthology
paste R1_files.txt R2_files.txt names.txt | while read col1 col2 col3 ; do echo megahit -1 /mnt/gasablok2/C_plate_genomes/C_plate_reads/C_plate_cleaned_reads/${col1}.cleaned.fastq -2 /mnt/gasablok2/C_plate_genomes/C_plate_reads/C_plate_cleaned_reads/${col2}.cleaned.fastq --k-list 79,99,119 --memory 0.7 --num-cpu-threads 70 --out-dir "${col3}"_megahit_assembly --tmp-dir /mnt/gasablok2/C_plate_genomes/C_plate_reads/C_plate_cleaned_reads/temp_megahit_assembly_C_plate; done
paste 1_files.txt 2_files.txt names.txt | while read col1 col2 col3 ; do echo bowtie2 -t -x Fungal_MT_genomes -p 60 --very-sensitive-local -1 ${col1} -2 ${col2} -S "${col3}".genome.mapping.sam --no-unal --al-conc ${col3}.genome.aligned.fastq; done
paste R1_files.txt R2_files.txt names.txt | while read col1 col2 col3 ; do echo bowtie2 -t -x Fungal_MT_CDS -p 60 --very-sensitive-local -1 ${col1} -2 ${col2} -S "${col3}".CDS.mapping.sam --no-unal --al-conc ${col3}.CDS.aligned.fastq; done
paste R1_files.txt R2_files.txt names.txt | while read col1 col2 col3 ; do echo bowtie2 -t -x Fungal_MT_CDS -p 60 --very-sensitive-local -1 ${col1} -2 ${col2} -S "${col3}".CDS.mapping.sam --no-unal --al-conc ${col3}.CDS.aligned.fastq; done
paste list.R1.txt list.R2.txt names.txt | while read col1 col2 col3 ; do echo megahit -1 ${col1} -2 ${col2} --k-list 79,99,119 --memory 0.7 --num-cpu-threads 70 --out-dir "${col3}"_assembly --tmp-dir /mnt/long_term_storage/fungal_previous_sequencing_runs_MEGAHIT/tmp_megahit
for f in *.fasta; do echo $f; done > fasta_files.txt
for f in *.fasta.column.txt; do echo $f; done > fasta_column.txt
cat fasta_files.txt | cut -f 1 -d "." > names.txt
paste fasta_files.txt fasta_column.txt names.txt | while read col1 col2 col3 ; do echo pullseq -i ${col1} -n ${col2} $$ "${col3}".extraction.fasta; done
gmes_petap.pl --ES --max_intron 3000 --soft_mask 2000 --cores 60 --sequence A006_COMBINED_ABYSS_SPADES_polished_3.fasta.header.masked.fasta.genome.masked.fasta.softmasked.fa --fungus
sed -i -e 's/^>jgi|/>/' 45_Hastodontia_aff_hastata_assembly.extraction.fasta
for f in *.fasta; do echo java -jar -Xmx100g macse_v2.03.jar -prog alignSequences -gc_def 12 -seq $f -out_AA $f.AA -out_NT $f.NT @ $f.MASCSE.run.txt; done
for f in *.fasta; do echo muscle -in $f -out $f.alignment.fasta; done
for f in *.fasta; do echo mafft --thread 40 --threadtb 30 --threadit 0 --reorder --auto $f @ $f.MAFFT.fasta; done
for f in *.fasta; do echo prank -d=$f -o=$f.PRANK.alignment.fasta; done
for f in *.fasta; do echo quast $f -o $f.output.quast -m 500 -t 40 -e --circos -f --rna-finding; done
for f in *.fa; do echo lastz RPB1_Final.fasta[multiple] $f --output=$f.RPB1.txt --ambiguous=iupac --coverage=70 --format=general; done
for f in *.fasta; do echo iqtree2 --seqtype DNA -s $f --alrt 1000 -b 1000 -T 40; done 
iqtree2 --seqtype DNA -s Final_ITS_Phylogeny_headers.fasta.PRANK.alignment.fasta.best.fas.phy --alrt 1000 -b 1000 -T 40
iqtree2 --seqtype DNA -s Final_RPB2_Phylogeny_headers.fasta.PRANK.alignment.fasta.best.fas.phy --alrt 1000 -b 1000 -T 40
iqtree2 --seqtype DNA -s Final_TEF1_Phylogeny_headers.fasta.PRANK.alignment.fasta.best.fas.phy --alrt 1000 -b 1000 -T 40
iqtree2 --seqtype DNA -s Revised_Final_RPB1_Phylogeny_headers.fasta.PRANK.alignment.fasta.best.fas.phy --alrt 1000 -b 1000 -T 40
for f in *.fasta; do echo modeltest-ng -i $f -d aa -o $f.model.fasta -p 40 -t ml; done 
for f in *.1.fastq; do echo $f; done >> R1_files.txt
for f in *.2.fastq; do echo $f; done >> R2_files.txt
cat R1_files.txt | cut -f 1,2 -d "." > names.txt
paste R1_files.txt R2_files.txt names.txt | while read col1 col2 col3 ; do echo megahit -1 ${col1} -2 ${col2} --k-list 59,69,79  --memory 0.7 --num-cpu-threads 70 --out-dir "${col3}"_SSU_megahit_assembly --tmp-dir /mnt/gasablok2/C_plate_genomes/C_plate_reads/C_plate_mitochondrial_genomes/C_plate_mitochondrial_gene_specific_assembly/LSU_SSU_assembly/SSU_aligned/tmp; done
paste names1.txt names2.txt | while read col1 col2; do echo lastz ${col2} {col1} --ambiguous=iupac --coverage=70 --output="$line".output.alignment --format=general; done
paste R1_files.txt R2_files.txt list.txt | while read col1 col2 col3 ; do echo abyss-pe j=10 k=75 name=${col3} 'in={col1} {col2}' done > {col3}_abyss_assembly_K75_log.txt; done 
paste R1_files.txt R2_files.txt list.txt | while read col1 col2 col3 ; do echo abyss-pe j=10 k=85 name=${col3} 'in={col1} {col2}' done > {col3}_abyss_assembly_K85_log.txt; done 
paste R1_files.txt R2_files.txt list.txt | while read col1 col2 col3 ; do echo abyss-pe j=10 k=95 name=${col3} 'in={col1} {col2}' done > {col3}_abyss_assembly_K95_log.txt; done 
paste R1_files.txt R2_files.txt names.txt | while read col1 col2 col3 ; do echo bowtie2 -t -x hyme -p 60 --very-sensitive-local -1 ${col1} -2 ${col2} -S ${col3}.sam --no-unal --al-conc ${col3}.aligned.fastq; done
for f in *.fasta; do echo python funannotate.py predict -i $f --force -s "$f"_annotations --name "$f"_AN -o "$f"_ANF --protein_evidence All.proteins_agaromycotina.fasta --augustus_species yarrowia_lipolytica --cpus 40; done
for f in *.fa; do echo python funannotate.py predict -i $f --force -s "$f"_annotations --name "$f"_AN -o "$f"_ANF --protein_evidence All.proteins_agaromycotina.fasta --augustus_species yarrowia_lipolytica --cpus 40; done 
raxmlHPC-PTHREADS-SSE3 -f a -# 100 -T 40 -m PROTGAMMAAUTO -p 12345 -x 12345 -s mycoglea_complete_ortholog_inpoaralog_concatenation.clipkit.fasta.phylip -n mycoglea_complete_ortholog_inpoaralog_concatenation.clipkit.fasta.tree

cat names.txt | while read line; do echo bowtie2-build "$line"_round2_polished_genome.fasta $line; done && cat names.txt | while read line; do echo bowtie2 -t -p 75 -x $line --very-sensitive-local -1 /mnt/long_term_storage/fungal_gene_based_assembly/original_fastq_files/"$line"_L001_R1_001.fastq -2 /mnt/long_term_storage/fungal_gene_based_assembly/original_fastq_files/"$line"_L001_R2_001.fastq --no-unal -S "$line".sam; done  && cat names.txt | while read line; do echo samtools view -bS --threads 60 "$line".sam -o "$line".bam; done && cat names.txt | while read line; do echo samtools sort --threads 60 "$line".bam -o "$line"_sort.bam; done && cat names.txt | while read line; do echo samtools index "$line"_sort.bam; done && cat names.txt | while read line; do echo java -Xmx100g -jar pilon-1.24.jar --genome "$line"_round1_polished_CDS.fasta --frags "$line"_sort.bam --output "$line"_round2_polished_CDS --threads 60 --changes --fix all --mindepth 10; done
cat names.txt | while read line; do echo bowtie2-build "$line"_round1_polished_genome.fasta $line; done && cat names.txt | while read line; do echo bowtie2 -t -p 75 -x $line --very-sensitive-local -1 /mnt/long_term_storage/fungal_mitochondrial_assemblies/genome_aligned_fastq_files/"$line".genome.aligned.1.fastq -2 /mnt/long_term_storage/fungal_mitochondrial_assemblies/genome_aligned_fastq_files/"$line".genome.aligned.2.fastq --no-unal -S "$line".sam; done  && cat names.txt | while read line; do echo samtools view -bS --threads 60 "$line".sam -o "$line".bam; done && cat names.txt | while read line; do echo samtools sort --threads 60 "$line".bam -o "$line"_sort.bam; done && cat names.txt | while read line; do echo samtools index "$line"_sort.bam; done && cat names.txt | while read line; do echo java -Xmx100g -jar pilon-1.24.jar --genome "$line"_round1_polished_genome.fasta --frags "$line"_sort.bam --output "$line"_round2_polished_genome --threads 60 --changes --fix all --mindepth 10; done
cat names.txt | while read line; do echo bowtie2-build "$line".fasta $line; done && cat names.txt | while read line; do echo bowtie2 -t -p 75 -x $line --very-sensitive-local -1 /mnt/long_term_storage/fungal_mycoglea_phylogeny/"$line"_L001_R1_001.fastq -2 /mnt/long_term_storage/fungal_mycoglea_phylogeny/"$line"_L001_R2_001.fastq --no-unal -S "$line".sam; done  && cat names.txt | while read line; do echo samtools view -bS --threads 60 "$line".sam -o "$line".bam; done && cat names.txt | while read line; do echo samtools sort --threads 60 "$line".bam -o "$line"_sort.bam; done && cat names.txt | while read line; do echo samtools index "$line"_sort.bam; done && cat names.txt | while read line; do echo java -Xmx100g -jar pilon-1.24.jar --genome "$line".scaffolds.fasta --frags "$line"_sort.bam --output "$line"_round1_polished_genome --threads 60 --changes --fix all --mindepth 10; done
cat names.txt | while read line; do echo bowtie2-build "$line".final.ABYSS.fasta.1000bp.fasta $line; done && cat names.txt | while read line; do echo bowtie2 -t -p 75 -x $line --very-sensitive-local -1 /mnt/long_term_storage/fungal_mycoglea_phylogeny/"$line"_L001_R1_001.fastq -2 /mnt/long_term_storage/fungal_mycoglea_phylogeny/"$line"_L001_R2_001.fastq --no-unal -S "$line".sam; done  && cat names.txt | while read line; do echo samtools view -bS --threads 60 "$line".sam -o "$line".bam; done && cat names.txt | while read line; do echo samtools sort --threads 60 "$line".bam -o "$line"_sort.bam; done && cat names.txt | while read line; do echo samtools index "$line"_sort.bam; done && cat names.txt | while read line; do echo java -Xmx100g -jar pilon-1.24.jar --genome "$line".final.ABYSS.fasta.1000bp.fasta --frags "$line"_sort.bam --output "$line"_round1_polished_genome --threads 60 --changes --fix all --mindepth 10; done
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}'

grep ^orf output.glimmer | awk '{OFS="\t"}{strang = "+"}{if($4 < 0) strang="-"}{gsub(/[+-]/," ")}{print "FASTA_HEADER", "GLIMMER", "gene" , $2 , $3, $5, strang , $4, "ID="$1"; NOTE:GLIMMER ORF prediction;"}'



paste R1_files.txt R2_files.txt names.txt | while read col1 col2 col3 ; do echo spades.py -1 ${col1} -2 ${col2} --careful --threads 60 --tmp-dir ${col3}_tmp -k 47,57,67,77,87 -o ${col3}_assembly_C_plate_spades; done 


fastp --in1 --out1 --in2 --out2 --html --thread


#!/usr/bin/env python
02
import sys
03
from Bio.SeqRecord import SeqRecord
04
from Bio import SeqIO
05
 
06
glimmer_file = sys.argv[1]
07
fasta_file = sys.argv[2]
08
 
09
# Read the sequence file
10
seq_record = SeqIO.parse(open(fasta_file),"fasta").next()
11
 
12
outseqrecords = []
13
# Read the glimmer file, record by record
14
for inline in file(glimmer_file):
15
    if '>' in inline:
16
        seqname = inline.split()[0][1:]
17
        outfilename = "%s_g3.tfa" % (seqname)
18
        continue
19
    if "orf" not in inline:
20
        continue
21
    orfname, sbegin, send, rf, score = inline.strip().split()
22
    sbegin = int(sbegin)
23
    send = int(send)
24
    rf = int(rf)
25
    # reverse complement
26
    if rf < 0:
27
        sbegin, send = send, sbegin
28
    sbegin -= 1     # Python indexes start a 0
29
    score = float(score)
30
    # split the sequence record
31
    newseq = seq_record.seq[sbegin:send]
32
    if rf < 0:
33
        newseq = newseq.reverse_complement()
34
    # Add a sequence record to the output
35
    seqrecord_description = "begin=%d end=%d rf=%d score=%.2f" % (sbegin+1, send, rf, score)
36
    outseqrecords.append(SeqRecord(newseq,id=seqname+"_"+orfname, description=seqrecord_description))
37
 
38
SeqIO.write(outseqrecords,open(outfilename,"w"),"fasta")


awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}'awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}'

ls *.fa | parallel prank -d={} -o={.}.aligned.fasta; done 
ls *.fa | parallel Gblocks {} -t=p -e=-gb1 -b4=6 -s=y -p=y -d=y; done
selecton -i {} -e 0.1  -n 8 -l {.}.log -r {.}.results -o {.}.likelihood -t {.}.rasmol -c {.}.color
paste fasta_files.txt list.txt | while read col1 col2 ; do echo pullseq -i ${col1} -n ${col2} $$ "${col3}".extraction.fasta; done
paste fasta_files.txt list.txt | while read col1 col2 ; do echo sed 's/^>/>${col2}_/' ${col1}; done 
for f in *.fasta; do echo $f; done > fasta_files.txt
cat fasta_files.txt| cut -f 1,2,3 -d "_" > list.txt
paste fasta_files.txt list.txt | while read col1 col2 ; do echo sed "s/^>/>${col2}_/d" ${col1}; done
paste R1_files.txt R2_files.txt names.txt | while read col1 col2 col3 ; do echo bowtie2 -t -x LSU -p 60 --very-sensitive-local -1 ${col1} -2 ${col2} -S ${col3}.LSU.sam --no-unal --al-conc ${col3}.aligned.LSU.fastq; done



/software/SOAPdenovo2-r241/SOAPdenovo-fusion -D -s config -p 40 -K 41 -g k41 -c ../final.contigs.fa
/software/SOAPdenovo2-r241/SOAPdenovo-127mer map -s config -p 40 -g k41
/software/SOAPdenovo2-r241/SOAPdenovo-127mer scaff -p 40 -g k41

/software/SOAPdenovo2-r241/SOAPdenovo-fusion -D -s config -p 40 -K 41 -g k41_1 -c ../final.contigs.fa
/software/SOAPdenovo2-r241/SOAPdenovo-127mer map -s config -p 40 -g k41_1
/software/SOAPdenovo2-r241/SOAPdenovo-127mer scaff -p 40 -g k41_1 -F


#maximal read length
max_rd_len=151
[LIB]
avg_ins=500
reverse_seq=0
asm_flags=2
#in which order the reads are used while scaffolding
rank=1
# cutoff of pair number for a reliable connection (at least 3 for short insert size)
pair_num_cutoff=3
#minimum aligned length to contigs for a reliable read location (at least 32 for short insert size)
map_len=32
#a pair of fastq file, read 1 file should always be followed by read 2 file
q1=/DenovoSeq/trimmomatic/500B_R_1P.fastq
q2=/DenovoSeq/trimmomatic/500B_R_2P.fastq
[LIB]
#average insert size
avg_ins=800
#if sequence needs to be reversed
reverse_seq=0
#in which part(s) the reads are used
asm_flags=2
#in which order the reads are used while scaffolding
rank=3
# cutoff of pair number for a reliable connection (at least 3 for short insert size)
pair_num_cutoff=3
#minimum aligned length to contigs for a reliable read location (at least 32 for short insert size)
map_len=32
#a pair of fastq file, read 1 file should always be followed by read 2 file
q1=/DenovoSeq/trimmomatic/800B_R_1P.fastq
q2=/DenovoSeq/trimmomatic/800B_R_2P.fastq
[LIB]
avg_ins=3000
reverse_seq=1
asm_flags=2
rank=3
# cutoff of pair number for a reliable connection (at least 5 for large insert size)
pair_num_cutoff=4
#minimum aligned length to contigs for a reliable read location (at least 35 for large insert size)
map_len=35
q1=/DenovoSeq/trimmomatic/3k_1_R_1P.fastq
q2=/DenovoSeq/trimmomatic/3k_1_R_2P.fastq
[LIB]
avg_ins=5000
reverse_seq=1
asm_flags=2
rank=4
# cutoff of pair number for a reliable connection (at least 5 for large insert size)
pair_num_cutoff=5
#minimum aligned length to contigs for a reliable read location (at least 35 for large insert size)
map_len=35
q1=/DenovoSeq/trimmomatic/5k-1_R_1P.fastq
q2=/DenovoSeq/trimmomatic/5k-1_R_2P.fastq
[LIB]
avg_ins=5000
reverse_seq=1
asm_flags=2
rank=4
# cutoff of pair number for a reliable connection (at least 5 for large insert size)
pair_num_cutoff=5
#minimum aligned length to contigs for a reliable read location (at least 35 for large insert size)
map_len=35
q1=/DenovoSeq/trimmomatic/5k-2_R_1P.fastq
q2=/DenovoSeq/trimmomatic/5k-2_R_2P.fastq
[LIB]
avg_ins=10000
reverse_seq=1
asm_flags=2
rank=5
# cutoff of pair number for a reliable connection (at least 5 for large insert size)
pair_num_cutoff=5
#minimum aligned length to contigs for a reliable read location (at least 35 for large insert size)
map_len=35
q1=/DenovoSeq/trimmomatic/10k_R_1P.fastq
q2=/DenovoSeq/trimmomatic/10k_R_2P.fastq


Example usage:

MEM=128 # 128gb
BASE=STRAINX
BBTools - https://jgi.doe.gov/data-and-tools/bbtools/
Trimmomatic - http://www.usadellab.org/cms/?page=trimmomatic (Optional)
bowtie2 - http://bowtie-bio.sourceforge.net/bowtie2/index.shtml (Optional)
bwa - https://github.com/lh3/bwa
Pilon - https://github.com/broadinstitute/pilon/wiki
sourmash - https://sourmash.readthedocs.io/ (install via conda/pip)
NCBI BLAST+ - ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST
minimap2 - https://github.com/lh3/minimap2 Assemblers
SPAdes - http://cab.spbu.ru/software/spades/
megahit - https://github.com/voutcn/megahit
dipspades - (SPAdes 3.11.1 - note it is not part of later SPAdes packages) http://cab.spbu.ru/files/release3.11.1/dipspades_manual.html



READSDIR=reads
TRIMREAD=reads_trimmed
CPU=8
AAFTF trim --method bbduk --memory $MEM -c $CPU \
 --left $READSDIR/${BASE}_R1.fq.gz --right $READSDIR/${BASE}_R2.fq.gz \
  -o $TRIMREAD/${BASE}
# this step make take a lot of memory depending on how many filtering libraries you use
AAFTF filter -c $CPU --memory $MEM --aligner bbduk \
	  -o $TRIMREAD/${BASE} --left $TRIMREAD/${BASE}_1P.fastq.gz --right $TRIMREAD/${BASE}_2P.fastq.gz
Assembly

The specified assembler can be made through the --method option. The full set of options are below.

usage: AAFTF assemble [-h] [-q] [--method METHOD] -o OUT [-w WORKDIR]
                      [-c cpus] [-m MEMORY] [-l LEFT] [-r RIGHT] [-v]
                      [--tmpdir TMPDIR] [--assembler_args ASSEMBLER_ARGS]
                      [--haplocontigs] [--pipe]

Run assembler on cleaned reads

optional arguments:
  -h, --help            show this help message and exit
  -q, --quiet           Do not output warnings to stderr
  --method METHOD       Assembly method: spades, dipspades, megahit
  -o OUT, --out OUT     Output assembly FASTA
  -w WORKDIR, --workdir WORKDIR
                        assembly output directory
  -c cpus, --cpus cpus  Number of CPUs/threads to use.
  -m MEMORY, --memory MEMORY
                        Memory (in GB) setting for SPAdes. Default is 32
  -l LEFT, --left LEFT  Left (Forward) reads
  -r RIGHT, --right RIGHT
                        Right (Reverse) reads
  -v, --debug           Print Spades stdout to terminal
  --tmpdir TMPDIR       Assembler temporary dir
  --assembler_args ASSEMBLER_ARGS
                        Additional SPAdes/Megahit arguments
  --haplocontigs        For dipSPAdes take the haplocontigs file
  --pipe                AAFTF is running in pipeline mode
CPU=24
MEM=96
LEFT=$TRIMREAD/${BASE}_filtered_1.fastq.gz
RIGHT=$TRIMREAD/${BASE}_filtered_2.fastq.gz
WORKDIR=working_AAFTF
OUTDIR=genomes
ASMFILE=$OUTDIR/${BASE}.spades.fasta
mkdir -p $WORKDIR $OUTDIR
AAFTF assemble -c $CPU --mem $MEM \
	  --left $LEFT --right $RIGHT  \
	   -o $ASMFILE -w $WORKDIR/spades_$BASE
vectrim

CPU=16
MEM=16
LEFT=$TRIMREAD/${BASE}_filtered_1.fastq.gz
RIGHT=$TRIMREAD/${BASE}_filtered_2.fastq.gz
WORKDIR=working_AAFTF
OUTDIR=genomes
ASMFILE=$OUTDIR/${BASE}.spades.fasta
VECTRIM=$OUTDIR/${BASE}.vecscreen.fasta
mkdir -p $WORKDIR $OUTDIR
AAFTF vecscreen -c $CPU -i $ASMFILE -o $VECTRIM


raxmlHPC-PTHREADS-SSE3 -f a -# 100 -T 70 -m PROTGAMMAAUTO -p 12345 -x 12345 -s Auricularales_SpeciesTreeAlignment.trimmed.fasta.phylip -n Auricularales_SpeciesTreeAlignment.trimmed.fasta.phylip.PROTGAMMAAUTO.tree
raxmlHPC-PTHREADS-SSE3 -f a -# 100 -T 70 -m PROTGAMMAIGTR -p 12345 -x 12345 -s Auricularales_SpeciesTreeAlignment.trimmed.fasta.phylip -n Auricularales_SpeciesTreeAlignment.trimmed.fasta.phylip.PROTGAMMAIGTRtree
raxmlHPC-PTHREADS-SSE3 -f a -# 100 -T 70 -m PROTGAMMAAUTO -p 12345 -x 12345 -s Hymneochaeteles_SpeciesTreeAlignment.trimmed.fasta.phylip -n Hymneochaeteles_SpeciesTreeAlignment.trimmed.fasta.phylip.PROTGAMMAAUTO.tree
raxmlHPC-PTHREADS-SSE3 -f a -# 100 -T 70 -m PROTGAMMAIGTR -p 12345 -x 12345 -s Hymneochaeteles_SpeciesTreeAlignment.trimmed.fasta.phylip -n Hymneochaeteles_SpeciesTreeAlignment.trimmed.fasta.phylip.phylip.PROTGAMMAIGTRtree


abyss			cairo			fastani			hdf5			libharu			mpdecimal		pcre			raxml-ng		stringtie
adam			cannoli			fastp			hh-suite		libidn			mpfr			pcre2			rcorrector		subread
any2fasta		canu			fastq-pair		hisat2			libidn2			mreps			peat			readline		superlu
apache-spark		cap3			fastqc			hmmer			libmagic		mummer			perl			recon			swi-prolog
aragorn			capnp			fasttree		htop			libmpc			muscle			phlawd			recon-ng		swig
arcs			cd-hit			fftw			htsbox			libmuscle		mysql@5.7		phylobayes		rename			szip
arks			cegma			filtlong		htslib			libomp			nanofilt		phyml			repeatmasker		tagdust
armadillo		cgns			flash			hwloc			libpng			nanoflann		phyx			rmblast			taxonkit
arpack			circlator		flex			icu4c			libssh2			nanopolish		picard-tools		rtg-tools		tbb
augustus		circos			flye			idba			libtiff			ncurses			pilercr			salmid			tbl2asn
autoconf		clapack			fontconfig		igraph			libtool			newick-utils		pilon			salt			tcl-tk
autoconf-archive	cluster-picker		fraggenescan		inetutils		libunistring		newicktools		pirate			samclip			transdecoder
automake		cmake			freebayes		infernal		libxml2			ngmaster		pixman			samtools		transpose
bali-phy		consel			freetype		isl			libyaml			ngmerge			pkg-config		samtools@0.1		treepl
bamtools		cpanminus		gappa			ispcr			links-scaffolder	nlopt			plink2			scipy			trf
barrnap			cufflinks		gatk			ivar			lmdb			node			porechop		sepp			trimal
bazam			cutadapt		gblocks			jellyfish		lrsim			notung			pplacer			seq-gen			trimmomatic
bbtools			dehomopolymerate	gcc			jpeg			lsd2			ntcard			prank			seqan			trnascan
bcftools		delly			gcc@8			k8			lz4			ntedit			prodigal		seqan@3			unicycler
beast2			dextractor		gd			kaiju			lzo			nthits			prokka			seqkit			unikmer
bedtools		diamond			gdbm			kalign			m4			ntjoin			proteinortho		seqtk			unixodbc
berkeley-db		diffr			gemma			kallisto		mafft			numpy			protobuf		sextractor		vcf2phylip
berokka			diffutils		geneid			kat			maftools		nxrepair		pullseq			sga			vcflib
bfc			disty			genewise		kent-tools		maker			nxtrim			pybind11		sibelia			veclibfort
bifrost			docbook			gepard			kma			mapcaller		oma			python@2		skesa			velvet
biomake			dos2unix		gettext			kmc			mash			open-mpi		python@3.8		skewer			virulign
bioperl			doxygen			git			kmergenie		matplotlib		openblas		python@3.9		smalt			vmatch
blast			dsh-bio			glib			kounta			maxbin2			openjdk			quast			snap			vsearch
blast@2.2		easel			glimmer3		kraken			mcl			openssl@1.1		quickmerge		snippy			vt
blat			edena			glimmerhmm		kraken2			megahit			orthofinder		quicktree		snoscan			webp
bless			eigen			glpk			last			meme			ossp-uuid		r			snp-dists		wget
boost			elph			gmap-gsnap		lastz			minced			paml			r8s			snp-sites		wiggletools
bowtie			emboss			gmcloser		libarchive		minia			pandaseq		racon			snpeff			xmatchview
bowtie2			epa-ng			gmp			libb2			miniasm			panito			rampart			spades			xssp
bracken			exabayes		google-sparsehash	libbigwig		minigraph		parallel		rasusa			sqlite			xz
busco			exonerate		gperftools		libevent		minimap			parsnp			rate4site		squeakr			zeromq
bwa			expat			gsl			libffi			minimap2		pathd8			raven-assembler		sratoolkit		zlib
bwa-mem2		express			harvest-tools		libgit2			mir-prefer		pcap			raxml			star-aligner		zstd


export PATH="/mnt/gasablok2/fungal_genome_completeness_w2rap_assemblies/augustus-3.2.1/bin:$PATH"
export PATH="/mnt/gasablok2/fungal_genome_completeness_w2rap_assemblies/augustus-3.2.1/scripts:$PATH"
export AUGUSTUS_CONFIG_PATH="/mnt/gasablok2/fungal_genome_completeness_w2rap_assemblies/augustus-3.2.1/config/"


/mnt/gasablok2/fungal_genome_completeness_w2rap_assemblies/augustus-3.2.1/bin:$PATH
/mnt/gasablok2/fungal_genome_completeness_w2rap_assemblies/augustus-3.2.1/config:$PATH
/mnt/gasablok2/fungal_genome_completeness_w2rap_assemblies/augustus-3.2.1/scripts



quast blasia_2020_hybrid_polishing.fasta -r Mpolymorpha_320_v3.0.hardmasked.fa -g Mpolymorpha_320_v3.1.gene.gff3 --fragmented -o blasia_2020_hybrid_polishing.fasta.output.quast.gene.finding -m 500 -t 40 -e --circos --glimmer
quast blasia_canu_version_2020_polishing.fasta -r Mpolymorpha_320_v3.0.hardmasked.fa -g Mpolymorpha_320_v3.1.gene.gff3 --fragmented -o blasia_canu_version_2020_polishing.fasta.output.quast.gene.finding -m 500 -t 40 -e --circos --glimmer
quast blasia_wtdbg2_assembly_2020_polishing.fasta -r Mpolymorpha_320_v3.0.hardmasked.fa -g Mpolymorpha_320_v3.1.gene.gff3 --fragmented -o blasia_wtdbg2_assembly_2020_polishing.fasta.output.quast.gene.finding -m 500 -t 40 -e --circos --glimmer
--glimmer






run_BUSCO.py -i blasia_2020_hybrid_polishing.fasta -c 30 -o blasia_2020_busco -e 1e-10 -m DNA -l -sp ; done  


Options:
-o  --output-dir  <dirname>       Directory to store all result files [default: quast_results/results_<datetime>]
-r                <filename>      Reference genome file
-g  --features [type:]<filename>  File with genomic feature coordinates in the reference (GFF, BED, NCBI or TXT)
                                  Optional 'type' can be specified for extracting only a specific feature type from GFF
-m  --min-contig  <int>           Lower threshold for contig length [default: 500]
-t  --threads     <int>           Maximum number of threads [default: 25% of CPUs]

Advanced options:
-s  --split-scaffolds                 Split assemblies by continuous fragments of N's and add such "contigs" to the comparison
-l  --labels "label, label, ..."      Names of assemblies to use in reports, comma-separated. If contain spaces, use quotes
-L                                    Take assembly names from their parent directory names
-e  --eukaryote                       Genome is eukaryotic (primarily affects gene prediction)
    --fungus                          Genome is fungal (primarily affects gene prediction)
    --large                           Use optimal parameters for evaluation of large genomes
                                      In particular, imposes '-e -m 3000 -i 500 -x 7000' (can be overridden manually)
-k  --k-mer-stats                     Compute k-mer-based quality metrics (recommended for large genomes)
                                      This may significantly increase memory and time consumption on large genomes
    --k-mer-size                      Size of k used in --k-mer-stats [default: 101]
    --circos                          Draw Circos plot
-f  --gene-finding                    Predict genes using GeneMarkS (prokaryotes, default) or GeneMark-ES (eukaryotes, use --eukaryote)
    --mgm                             Use MetaGeneMark for gene prediction (instead of the default finder, see above)
    --glimmer                         Use GlimmerHMM for gene prediction (instead of the default finder, see above)
    --gene-thresholds <int,int,...>   Comma-separated list of threshold lengths of genes to search with Gene Finding module
                                      [default: 0,300,1500,3000]
    --rna-finding                     Predict ribosomal RNA genes using Barrnap
-b  --conserved-genes-finding         Count conserved orthologs using BUSCO (only on Linux)
    --operons  <filename>             File with operon coordinates in the reference (GFF, BED, NCBI or TXT)
    --est-ref-size <int>              Estimated reference size (for computing NGx metrics without a reference)
    --contig-thresholds <int,int,...> Comma-separated list of contig length thresholds [default: 0,1000,5000,10000,25000,50000]
-u  --use-all-alignments              Compute genome fraction, # genes, # operons in QUAST v1.* style.
                                      By default, QUAST filters Minimap's alignments to keep only best ones
-i  --min-alignment <int>             The minimum alignment length [default: 65]
    --min-identity <float>            The minimum alignment identity (80.0, 100.0) [default: 95.0]
-a  --ambiguity-usage <none|one|all>  Use none, one, or all alignments of a contig when all of them
                                      are almost equally good (see --ambiguity-score) [default: one]
    --ambiguity-score <float>         Score S for defining equally good alignments of a single contig. All alignments are sorted
                                      by decreasing LEN * IDY% value. All alignments with LEN * IDY% < S * best(LEN * IDY%) are
                                      discarded. S should be between 0.8 and 1.0 [default: 0.99]
    --strict-NA                       Break contigs in any misassembly event when compute NAx and NGAx.
                                      By default, QUAST breaks contigs only by extensive misassemblies (not local ones)
-x  --extensive-mis-size  <int>       Lower threshold for extensive misassembly size. All relocations with inconsistency
                                      less than extensive-mis-size are counted as local misassemblies [default: 1000]
    --scaffold-gap-max-size  <int>    Max allowed scaffold gap length difference. All relocations with inconsistency
                                      less than scaffold-gap-size are counted as scaffold gap misassemblies [default: 10000]
    --unaligned-part-size  <int>      Lower threshold for detecting partially unaligned contigs. Such contig should have
                                      at least one unaligned fragment >= the threshold [default: 500]
    --skip-unaligned-mis-contigs      Do not distinguish contigs with >= 50% unaligned bases as a separate group
                                      By default, QUAST does not count misassemblies in them
    --fragmented                      Reference genome may be fragmented into small pieces (e.g. scaffolded reference)
    --fragmented-max-indent  <int>    Mark translocation as fake if both alignments are located no further than N bases
                                      from the ends of the reference fragments [default: 85]
                                      Requires --fragmented option
    --upper-bound-assembly            Simulate upper bound assembly based on the reference genome and reads
    --upper-bound-min-con  <int>      Minimal number of 'connecting reads' needed for joining upper bound contigs into a scaffold
                                      [default: 2 for mate-pairs and 1 for long reads]
    --est-insert-size  <int>          Use provided insert size in upper bound assembly simulation [default: auto detect from reads or 255]
    --plots-format  <str>             Save plots in specified format [default: pdf].
                                      Supported formats: emf, eps, pdf, png, ps, raw, rgba, svg, svgz
    --memory-efficient                Run everything using one thread, separately per each assembly.
                                      This may significantly reduce memory consumption on large genomes
    --space-efficient                 Create only reports and plots files. Aux files including .stdout, .stderr, .coords will not be created.
                                      This may significantly reduce space consumption on large genomes. Icarus viewers also will not be built
-1  --pe1     <filename>              File with forward paired-end reads (in FASTQ format, may be gzipped)
-2  --pe2     <filename>              File with reverse paired-end reads (in FASTQ format, may be gzipped)
    --pe12    <filename>              File with interlaced forward and reverse paired-end reads. (in FASTQ format, may be gzipped)
    --mp1     <filename>              File with forward mate-pair reads (in FASTQ format, may be gzipped)
    --mp2     <filename>              File with reverse mate-pair reads (in FASTQ format, may be gzipped)
    --mp12    <filename>              File with interlaced forward and reverse mate-pair reads (in FASTQ format, may be gzipped)
    --single  <filename>              File with unpaired short reads (in FASTQ format, may be gzipped)
    --pacbio     <filename>           File with PacBio reads (in FASTQ format, may be gzipped)
    --nanopore   <filename>           File with Oxford Nanopore reads (in FASTQ format, may be gzipped)
    --ref-sam <filename>              SAM alignment file obtained by aligning reads to reference genome file
    --ref-bam <filename>              BAM alignment file obtained by aligning reads to reference genome file
    --sam     <filename,filename,...> Comma-separated list of SAM alignment files obtained by aligning reads to assemblies
                                      (use the same order as for files with contigs)
    --bam     <filename,filename,...> Comma-separated list of BAM alignment files obtained by aligning reads to assemblies
                                      (use the same order as for files with contigs)
                                      Reads (or SAM/BAM file) are used for structural variation detection and
                                      coverage histogram building in Icarus
    --sv-bedpe  <filename>            File with structural variations (in BEDPE format)

Speedup options:
    --no-check                        Do not check and correct input fasta files. Use at your own risk (see manual)
    --no-plots                        Do not draw plots
    --no-html                         Do not build html reports and Icarus viewers
    --no-icarus                       Do not build Icarus viewers
    --no-snps                         Do not report SNPs (may significantly reduce memory consumption on large genomes)
    --no-gc                           Do not compute GC% and GC-distribution
    --no-sv                           Do not run structural variation detection (make sense only if reads are specified)
    --no-gzip                         Do not compress large output files
    --no-read-stats                   Do not align reads to assemblies
                                      Reads will be aligned to reference and used for coverage analysis,
                                      upper bound assembly simulation, and structural variation detection.
                                      Use this option if you do not need read statistics for assemblies.
    --fast                            A combination of all speedup options except --no-check

Other:
    --silent                          Do not print detailed information about each step to stdout (log file is not affected)
    --test                            Run QUAST on the data from the test_data folder, output to quast_test_output
    --test-sv                         Run QUAST with structural variants detection on the data from the test_data folder,
                                      output to quast_test_output
-h  --help                            Print full usage message
-v  --version                         Print version


bowtie2-build Final_RPB1_Phylogeny_headers_one.fasta Final_RPB1_Phylogeny
bowtie2-build Final_RPB2_Phylogeny_headers_one.fasta Final_RPB2_Phylogeny
bowtie2-build Final_TEF1_Phylogeny_headers_one.fasta Final_TEF1_Phylogeny

paste R1_files.txt R2_files.txt names.txt | while read col1 col2 col3 ; do echo bowtie2 -t -x Final_RPB1_Phylogeny -p 60 --very-sensitive-local -1 /mnt/gasablok2/C_plate_genomes/C_plate_reads/C_plate_cleaned_reads/${col1} -2 /mnt/gasablok2/C_plate_genomes/C_plate_reads/C_plate_cleaned_reads/${col2} -S ${col3}.RPB1.sam --no-unal --al-conc ${col3}.RPB1.aligned.fastq; done
paste R1_files.txt R2_files.txt names.txt | while read col1 col2 col3 ; do echo bowtie2 -t -x Final_RPB2_Phylogeny -p 60 --very-sensitive-local -1 /mnt/gasablok2/C_plate_genomes/C_plate_reads/C_plate_cleaned_reads/${col1} -2 /mnt/gasablok2/C_plate_genomes/C_plate_reads/C_plate_cleaned_reads/${col2} -S ${col3}.RPB2.sam --no-unal --al-conc ${col3}.RPB2.aligned.fastq; done
paste R1_files.txt R2_files.txt names.txt | while read col1 col2 col3 ; do echo bowtie2 -t -x Final_TEF1_Phylogeny -p 60 --very-sensitive-local -1 /mnt/gasablok2/C_plate_genomes/C_plate_reads/C_plate_cleaned_reads/${col1} -2 /mnt/gasablok2/C_plate_genomes/C_plate_reads/C_plate_cleaned_reads/${col2} -S ${col3}.TEF1.sam --no-unal --al-conc ${col3}.TEF1.aligned.fastq; done

paste R1_files.txt R2_files.txt names.txt | while read col1 col2 col3 ; do echo megahit -1 ${col1} -2 ${col2} --k-list 79,99,119 --memory 0.7 --num-cpu-threads 70 --out-dir "${col3}"_CDS_megahit_assembly --tmp-dir /mnt/gasablok2/C_plate_genomes/C_plate_reads/C_plate_mitochondrial_genomes/mitochondrial_aligned_fastq_files/CDS_aligned_fastq_files/tmp_CDS_directory; done
paste R1_files.txt R2_files.txt names.txt | while read col1 col2 col3 ; do echo megahit -1 ${col1} -2 ${col2} --k-list 79,99,119 --memory 0.7 --num-cpu-threads 70 --out-dir "${col3}"_genomes_megahit_assembly --tmp-dir /mnt/gasablok2/C_plate_genomes/C_plate_reads/C_plate_mitochondrial_genomes/mitochondrial_aligned_fastq_files/genomes_aligned_fastq_files/tmp_mito_genome; done
paste R1_files.txt R2_files.txt names.txt | while read col1 col2 col3 ; do echo megahit -1 ${col1} -2 ${col2} --k-list 79,99,119 --memory 0.7 --num-cpu-threads 70 --out-dir "${col3}"_RPB1_megahit_assembly --tmp-dir /mnt/gasablok2/C_plate_genomes/C_plate_reads/C_plate_marker_genes/marker_gene_aligned_files/RPB1_aligned_fastq/RPB1_assembly; done
paste R1_files.txt R2_files.txt names.txt | while read col1 col2 col3 ; do echo megahit -1 ${col1} -2 ${col2} --k-list 79,99,119 --memory 0.7 --num-cpu-threads 70 --out-dir "${col3}"_RPB2_megahit_assembly --tmp-dir /mnt/gasablok2/C_plate_genomes/C_plate_reads/C_plate_marker_genes/marker_gene_aligned_files/RPB2_aligned_fastq/RPB2_assembly; done
paste R1_files.txt R2_files.txt names.txt | while read col1 col2 col3 ; do echo megahit -1 ${col1} -2 ${col2} --k-list 79,99,119 --memory 0.7 --num-cpu-threads 70 --out-dir "${col3}"_TEF_megahit_assembly --tmp-dir /mnt/gasablok2/C_plate_genomes/C_plate_reads/C_plate_marker_genes/marker_gene_aligned_files/TEF1_aligned_fastq/TEF1_assembly; done
paste R1_files.txt R2_files.txt names.txt | while read col1 col2 col3 ; do echo megahit -1 ${col1} -2 ${col2} --k-list 79,99,119 --memory 0.7 --num-cpu-threads 70 --out-dir "${col3}"_ITS_C_plate_megahit_assembly --tmp-dir /mnt/gasablok2/C_plate_genomes/C_plate_reads/C_plate_ITS_specific_assembly/ITS_mapping_files/tmp; done




Steps: 
1. create a Blobdir from the assembly file
./blobtools create --fasta examples/assembly.fasta assembly_example
2. Run the blast or the diamond of the genome assembly against the nt database or the uniprot database
./blastn -subject nt -query assembly.fasta -out blast.out -outfmt 6
3. Add the BLAST hits 
./blobtools add --hits blast.out
4. Filter a BlobDir:
./blobtools filter --param length--Min=5000 --output assembly_filter_5000 assembly_example



# 1. Create a new BlobDir from a FASTA file:
    ./blobtools create --fasta examples/assembly.fasta BlobDir

    # 2. Create a new BlobDir from a BlobDB:
    ./blobtools create --blobdb examples/blobDB.json BlobDir

    # 3. Add Coverage data from a BAM file:
    ./blobtools add --cov examples/assembly.reads.bam BlobDir

    # 4. Assign taxonomy from BLAST hits:
    ./blobtools add add --hits examples/blast.out --taxdump ../taxdump BlobDir

    # 5. Add BUSCO results:
    ./blobtools add --busco examples/busco.tsv BlobDir

    # 6. Host an interactive viewer:
    ./blobtools host BlobDir

Filter a BlobDir:
./blobtools filter --param length--Min=5000 --output BlobDir_len_gt_5000 BlobDir


for dir in *_assembly; do echo cp -r $dir/final.contigs.fa $dir.fasta; done >> mitoCDS_assembly_move.sh

kpic: keeps parismony informative and constant sites
kpi: keep only parsimony informative sites


/mnt/gasablok2/fungal_genome_completeness_w2rap_assemblies/K96_assemblies_merging/177_Helicogloea_septifera_96.w2rap_assembly.fasta
/mnt/gasablok2/fungal_genome_completeness_w2rap_assemblies/K96_assemblies_merging/178_Helicogloea_aquilonia_96.w2rap_assembly.fasta
/mnt/gasablok2/fungal_genome_completeness_w2rap_assemblies/K96_assemblies_merging/179_Helicogloea_septifera_96.w2rap_assembly.fasta



blastn -db all.genome -query femsjoenia_filtered.contigs.fasta -outfmt "6 qseqid staxids bitscore std" -max_target_seqs 10 -max_hsps 1 -evalue 1e-25 -num_threads 40 -out femsjoenia.blast.out


177_Helicogloea_septifera_Achroomyces_assembly_msk.msk.fasta   2_Trichaptum_subchartaceum_S2_assembly_msk.msk.fasta           99_Flaviporellus_splitgerberi_S98_assembly_msk.msk.fasta
178_Helicogloea_aquilonia_Holotype_assembly_msk.msk.fasta      30_Rogersella_griselinae_S30_assembly_msk.msk.fasta            9_Fulvifomes_cf_rhytiphloeus_assembly_msk.msk.fasta
179_Helicogloea_septifera_Holotype_assembly_msk.msk.fasta



Parsing taxdump
Traceback (most recent call last):
  File "/mnt/gasablok2/blobtoolkit/blobtools2/lib/add.py", line 194, in <module>
    main()
  File "/mnt/gasablok2/blobtoolkit/blobtools2/lib/add.py", line 141, in main
    taxdump = fetch_taxdump(args["--taxdump"])
  File "/mnt/gasablok2/blobtoolkit/blobtools2/lib/fetch.py", line 93, in fetch_taxdump
    taxdump = Taxdump(path_to_taxdump)
  File "/mnt/gasablok2/blobtoolkit/blobtools2/lib/taxdump.py", line 28, in __init__
    self.load_ranks()
  File "/mnt/gasablok2/blobtoolkit/blobtools2/lib/taxdump.py", line 39, in load_ranks
    filename = os.path.abspath(os.path.join(self.directory, 'nodes.dmp'))
  File "/home/cloud-user/anaconda3/envs/btk_env/lib/python3.6/posixpath.py", line 80, in join
    a = os.fspath(a)
TypeError: expected str, bytes or os.PathLike object, not NoneType

paste R1_files.txt R2_files.txt names.txt | while read col1 col2 col3 ; do echo bowtie2 -t -x SSU -p 60 --very-sensitive-local -1 ${col1} -2 ${col2} -S ${col3}.SSU.sam --no-unal --al-conc ${col3}.aligned.SSU.fastq; done
blastn -db /mnt/gasablok2/blobtoolkit/nt -query femsjoenia_filtered.contigs.fasta  -outfmt "6 qseqid staxids bitscore std" -max_target_seqs 10 -max_hsps 1 -evalue 1e-25 -num_threads 75 -out femsjoenia_filtered.contigs.fasta.blast.out
blastn -db /mnt/gasablok2/blobtoolkit/nt -query All.Femzonia.Pacbio.fasta  -outfmt "6 qseqid staxids bitscore std" -max_target_seqs 10 -max_hsps 1 -evalue 1e-25 -num_threads 75 -out All.Femzonia.Pacbio.fasta.blast.out
for f in *.fasta; do echo python funannotate.py predict -i $f --force -s "$f"_AN --name "$f"_annotations_AN -s "$f"_annotation_AN --protein_evidence All.proteins_agaromycotina.fasta --augustus_species yarrowia_lipolytica --cpus 40; done >> auricularles_annotations.sh
cat names.txt  while read line; do echo bowtie-build "$line".scaffolds.fasta $line; done  
paste R1_files.txt R2_files.txt names.txt | while read col1 col2 col3 ; do echo bowtie2 -t -x "${col3}" -p 60 --very-sensitive-local -1 ${col1} -2 ${col2} -S "${col3}".sam --no-unal --al-conc ${col3}.aligned.fastq; done

#!/bin/bash -l
#SBATCH -J polytrichum_genome_assembly_version_5
#SBATCH --constraint="snb|hsw"
#SBATCH -p hugemem
#SBATCH -n 1
#SBATCH -c 40
#SBATCH --mem=1200000
#SBATCH --workdir=/proj/project_hy1057/polytrichum_genome/polytrichum_w2rap_contigger_assembly
#SBATCH --mail-type=END
#SBATCH --mail-user=sablokg@gmail.com

module load intel/18.0.1

export PATH=/proj/project_hy1057/chaetospiridum_nextseq/chaeto_w2rap_contigger_assembly/w2rap-contigger/bin:$PATH

#w2rap-contigger -t 40 -m 1000 -r /proj/project_hy1057/polytrichum_genome/PC3_Final.R1.fastq.cleaned.R1.fastq,/proj/project_hy1057/polytrichum_genome/PC3_Final.R2.fastq.cleaned.R2.fastq -o /proj/project_hy1057/polytrichum_genome/polytrichum_w2rap_contigger_assembly/poly_w2rap_contigger_assembly_80 -p poly_assembly_verion_5_80K -K 80 --tmp_dir /proj/project_hy1057/polytrichum_genome/polytrichum_w2rap_contigger_assembly/tmp_80K -s 2
#w2rap-contigger -t 40 -m 1000 -r /proj/project_hy1057/polytrichum_genome/PC3_Final.R1.fastq.cleaned.R1.fastq,/proj/project_hy1057/polytrichum_genome/PC3_Final.R2.fastq.cleaned.R2.fastq  -o /proj/project_hy1057/polytrichum_genome/polytrichum_w2rap_contigger_assembly/poly_w2rap_contigger_assembly_88 -p poly_assembly_verion_5_88K -K 88 --tmp_dir /proj/project_hy1057/polytrichum_genome/polytrichum_w2rap_contigger_assembly/tmp_88K -s 2
w2rap-contigger -t 40 -m 1000 -r /proj/project_hy1057/polytrichum_genome/PC3_Final.R1.fastq.cleaned.R1.fastq,/proj/project_hy1057/polytrichum_genome/PC3_Final.R2.fastq.cleaned.R2.fastq -o /proj/project_hy1057/polytrichum_genome/polytrichum_w2rap_contigger_assembly/poly_w2rap_contigger_assembly_96 -p poly_assembly_verion_5_96K -K 96 --tmp_dir /proj/project_hy1057/polytrichum_genome/polytrichum_w2rap_contigger_assembly/tmp_96K -s 2
mafft --auto --thread 10 All_version_1_mitochondrial_genome.fasta > All_version_1_mitochondrial_genome_alignment.fasta
module load gcc/7.3.0
export PATH=/wrk/gasablok/binaries:$PATH
export PATH=/wrk/gasablok/binaries/MECAT-canu/Linux-amd64/bin:$PATH

#mecat2ref -d /proj/project_hy1057/blasia_genome/blasia_organelle_final_assembly/blasia_pacbio_organelle/Combined_All_Runs_PacBio.fasta -r /proj/project_hy1057/blasia_genome/blasia_organelle_final_assembly/blasia_mitochondrial_genome_illumina_mapping/Mito_List_Upload.fasta -o blasia.all.pacbio.MECAT2REF.mitochondria.sam -w blasia.all.pacbio.mitochondria -t 24 -n 20 -n 50 -m 2 -x 0
#lordfast --index /proj/project_hy1057/blasia_genome/blasia_organelle_final_assembly/blasia_mitochondrial_genome_illumina_mapping/Mito_List_Upload.fasta
lordfast --search /proj/project_hy1057/blasia_genome/blasia_organelle_final_assembly/blasia_mitochondrial_genome_illumina_mapping/Mito_List_Upload.fasta --seq /proj/project_hy1057/blasia_genome/blasia_organelle_final_assembly/blasia_pacbio_organelle/Combined_All_Runs_PacBio.fasta --out blasia.all.mitochondria.pacbio.lordfast.sam --threads 24 --numMap 10
lordec-correct -T 10 -S lordec.stats.out -i blasia.embryophyte_mapping.bam.fasta -2 All.mitochondrial.mapping.fasta -k 21 -o blasia.embryophyte_mapping.bam_lordec_error_corrected.fasta -s 3

srun canu gridOptions="--time=24:00:00" -corMhapSensitivity=high -p blasia_mitochondrial_mapping_all_lordec_canu -d /proj/project_hy1057/blasia_genome/blasia_organelle_final_assembly/blasia_mitochondrial_pacbio_mapping/mitochondrial_mapping_all_lordec_canu genomeSize=200000 -minReadLength=500 -merylMemory=64 -gnuplotImageFormat=png -ovsThreads=$SLURM_CPUS_PER_TASK -ovbThreads=$SLURM_CPUS_PER_TASK -pacbio-raw All_mitochondrial_mapping_pacbio_lordec_error_corrected.fasta
lordfast --search /proj/project_hy1057/blasia_genome/blasia_organelle_final_assembly/blasia_mitochondrial_genome_illumina_mapping/Mito_List_Upload.fasta --seq /proj/project_hy1057/blasia_genome/blasia_organelle_final_assembly/blasia_pacbio_organelle/Combined_All_Runs_PacBio.fasta --out blasia.all.mitochondria.pacbio.lordfast.sam --threads 24 --numMap 10
lordec-correct -T 10 -S lordec.stats.out -i blasia.embryophyte_mapping.bam.fasta -2 All.mitochondrial.mapping.fasta -k 21 -o blasia.embryophyte_mapping.bam_lordec_error_corrected.fasta -s 3

export PATH=/wrk/gasablok/binaries:$PATH
export PATH=/wrk/gasablok/binaries/SPAdes-3.10.1-Linux/bin/:$PATH

FILE1=/proj/project_hy1057/blasia_genome/blasia_organelle_final_assembly/blasia_mitochondrial_genome_illumina_mapping/mitochondrial_mapped.R1.fastq
FILE2=/proj/project_hy1057/blasia_genome/blasia_organelle_final_assembly/blasia_mitochondrial_genome_illumina_mapping/mitochondrial_mapped.R2.fastq
FILE3=/proj/project_hy1057/blasia_genome/blasia_organelle_final_assembly/blasia_mitochondrial_hybrid_spades/blasia_mitochondrial_mapping_all_lordec_canu.correctedReads.fasta

spades.py -1 $FILE1 -2 $FILE2 --pacbio $FILE3 --careful --threads 24 --memory 256 --tmp-dir /proj/project_hy1057/blasia_genome/blasia_organelle_final_assembly/blasia_mitochondrial_hybrid_spades -k 77,97,121 -o /proj/project_hy1057/blasia_genome/blasia_organelle_final_assembly/blasia_mitochondrial_hybrid_spades/mitochondrial_hybrid_assembly

perl GapFiller.pl -l lib.txt -s Blasia_Mitochondrial_Assembly_Map.fasta -m 40 -o 5 -r 0.8 -n 10 -t 0 -i 10 -T 20 -b Blasia_mitochondrial_gap_filling
spades.py --pe1-1 $FILE1 --pe1-2 $FILE3 --pe1-s $FILE2 --pe1-s $FILE4 --pacbio $FILE5 --careful --threads 24 --memory 256 --tmp-dir /proj/project_hy1057/blasia_genome/blasia_organelle_final_assembly/blasia_organelle_hybrid_assembly -k 77,97,121 -o /proj/project_hy1057/blasia_genome/blasia_organelle_final_assembly/blasia_organelle_hybrid_assembly
lordec-correct -T 10 -S lordec.stats.out -i blasia.embryophyte_mapping.bam.fasta -2 All.mitochondrial.mapping.fasta -k 21 -o blasia.embryophyte_mapping.bam_lordec_error_corrected.fasta -s 3
FILE1=/proj/project_hy1057/blasia_genome/blasia_organelle_final_assembly/blasia_pacbio_organelle/Combined_All_Runs_PacBio.fasta
FILE2=Embryophyte.marchantia.fasta

lastz $FILE2[multiple] $FILE1 --identity=60 --coverage=50 --format=sam --output=blasia.embryophyte_mapping.sam


lastz ATP6.fas[multiple] Agahy1_GeneCatalog_CDS_20140307.fasta.final.fasta --output=Agahy1_GeneCatalog_CDS_20140307.alignment.ATP6.txt --ambiguous=iupac --coverage=70 --format=general
for f in *.fasta; do echo java -jar -Xmx100g macse_v2.03.jar -prog alignSequences -gc_def 12 -seq $f -out_AA $f.AA -out_NT $f.NT @ $f.MASCSE.run.txt; done
for f in *.fasta; do echo muscle -in $f -out $f.alignment.fasta; done
for f in *.fasta; do echo mafft --thread 40 --threadtb 30 --threadit 0 --reorder --auto $f @ $f.MAFFT.fasta; done
for f in *.fasta; do echo prank -d=$f -o=$f.PRANK.alignment.fasta; done

ava -jar -Xmx100g macse_v2.03.jar -prog alignSequences -gc_def 12 -seq Revised_Final_RPB1_Phylogeny_headers.fasta -out_AA Revised_Final_RPB1_Phylogeny_headers.fasta.AA -out_NT Revised_Final_RPB1_Phylogeny_headers.fasta.NT > Revised_Final_RPB1_Phylogeny_headers.fasta.MASCSE.run.txt

trimal -in Final_ITS_Phylogeny_headers.fasta.PRANK.alignment.fasta.best.fas -out Final_ITS_Phylogeny_headers.fasta.PRANK.alignment.fasta.best.fas.trimal.trimmed.fasta -strict
trimal -in Final_RPB2_Phylogeny_headers.fasta.PRANK.alignment.fasta.best.fas -out Final_RPB2_Phylogeny_headers.fasta.PRANK.alignment.fasta.best.fas.trimal.trimmed.fasta -strict
trimal -in Final_TEF1_Phylogeny_headers.fasta.PRANK.alignment.fasta.best.fas -out Final_TEF1_Phylogeny_headers.fasta.PRANK.alignment.fasta.best.fas.trimal.trimmed.fasta -strict
trimal -in Revised_Final_RPB1_Phylogeny_headers.fasta.PRANK.alignment.fasta.best.fas -out Revised_Final_RPB1_Phylogeny_headers.fasta.PRANK.alignment.fasta.best.fas.trimal.trimmed.fasta -strict

statal -in Final_ITS_Phylogeny_headers.fasta.MAFFT.fasta -sgt > Final_ITS_Phylogeny_headers.fasta.MAFFT.fasta.alignment.stats.txt
statal -in Final_RPB2_Phylogeny_headers.fasta.MAFFT.fasta -sgt > Final_RPB2_Phylogeny_headers.fasta.MAFFT.fasta.alignment.stats.txt
statal -in Final_TEF1_Phylogeny_headers.fasta.MAFFT.fasta -sgt > Final_TEF1_Phylogeny_headers.fasta.MAFFT.fasta.alignment.stats.txt
statal -in Revised_Final_RPB1_Phylogeny_headers.fasta.MAFFT.fasta -sgt > Revised_Final_RPB1_Phylogeny_headers.fasta.MAFFT.fasta.alignment.stats.txt
statal -in Final_ITS_Phylogeny_headers.fasta.MAFFT.fasta -scolidentt > Final_ITS_Phylogeny_headers.fasta.MAFFT.fasta.alignment.column.identity.txt
statal -in Final_RPB2_Phylogeny_headers.fasta.MAFFT.fasta -scolidentt > Final_RPB2_Phylogeny_headers.fasta.MAFFT.fasta.alignment.column.identity.txt
statal -in Final_TEF1_Phylogeny_headers.fasta.MAFFT.fasta -scolidentt > Final_TEF1_Phylogeny_headers.fasta.MAFFT.fasta.alignment.column.identity.txt
statal -in Revised_Final_RPB1_Phylogeny_headers.fasta.MAFFT.fasta -scolidentt > Revised_Final_RPB1_Phylogeny_headers.fasta.MAFFT.fasta.alignment.column.identity.txt

java -jar -Xmx100g macse_v2.03.jar -prog alignSequences -gc_def 12 -seq Final_ITS_Phylogeny_headers.fasta -out_AA Final_ITS_Phylogeny_headers.fasta.AA -out_NT Final_ITS_Phylogeny_headers.fasta.NT > Final_ITS_Phylogeny_headers.fasta.MASCSE.run.txt
java -jar -Xmx100g macse_v2.03.jar -prog alignSequences -gc_def 12 -seq Final_RPB2_Phylogeny_headers.fasta -out_AA Final_RPB2_Phylogeny_headers.fasta.AA -out_NT Final_RPB2_Phylogeny_headers.fasta.NT > Final_RPB2_Phylogeny_headers.fasta.MASCSE.run.txt
java -jar -Xmx100g macse_v2.03.jar -prog alignSequences -gc_def 12 -seq Final_TEF1_Phylogeny_headers.fasta -out_AA Final_TEF1_Phylogeny_headers.fasta.AA -out_NT Final_TEF1_Phylogeny_headers.fasta.NT > Final_TEF1_Phylogeny_headers.fasta.MASCSE.run.txt
java -jar -Xmx100g macse_v2.03.jar -prog alignSequences -gc_def 12 -seq Revised_Final_RPB1_Phylogeny_headers.fasta -out_AA Revised_Final_RPB1_Phylogeny_headers.fasta.AA -out_NT Revised_Final_RPB1_Phylogeny_headers.fasta.NT > Revised_Final_RPB1_Phylogeny_headers.fasta.MASCSE.run.txt
muscle -in Final_ITS_Phylogeny_headers.fasta -out Final_ITS_Phylogeny_headers.fasta.alignment.fasta
muscle -in Final_RPB2_Phylogeny_headers.fasta -out Final_RPB2_Phylogeny_headers.fasta.alignment.fasta
muscle -in Final_TEF1_Phylogeny_headers.fasta -out Final_TEF1_Phylogeny_headers.fasta.alignment.fasta
muscle -in Revised_Final_RPB1_Phylogeny_headers.fasta -out Revised_Final_RPB1_Phylogeny_headers.fasta.alignment.fasta
mafft --thread 40 --threadtb 30 --threadit 0 --reorder --auto Final_ITS_Phylogeny_headers.fasta > Final_ITS_Phylogeny_headers.fasta.MAFFT.fasta
mafft --thread 40 --threadtb 30 --threadit 0 --reorder --auto Final_RPB2_Phylogeny_headers.fasta > Final_RPB2_Phylogeny_headers.fasta.MAFFT.fasta
mafft --thread 40 --threadtb 30 --threadit 0 --reorder --auto Final_TEF1_Phylogeny_headers.fasta > Final_TEF1_Phylogeny_headers.fasta.MAFFT.fasta
mafft --thread 40 --threadtb 30 --threadit 0 --reorder --auto Revised_Final_RPB1_Phylogeny_headers.fasta > Revised_Final_RPB1_Phylogeny_headers.fasta.MAFFT.fasta
prank -d=Final_ITS_Phylogeny_headers.fasta -o=Final_ITS_Phylogeny_headers.fasta.PRANK.alignment.fasta
prank -d=Final_RPB2_Phylogeny_headers.fasta -o=Final_RPB2_Phylogeny_headers.fasta.PRANK.alignment.fasta
prank -d=Final_TEF1_Phylogeny_headers.fasta -o=Final_TEF1_Phylogeny_headers.fasta.PRANK.alignment.fasta
prank -d=Revised_Final_RPB1_Phylogeny_headers.fasta -o=Revised_Final_RPB1_Phylogeny_headers.fasta.PRANK.alignment.fasta

prank -d=Final_ITS_Phylogeny_headers.fasta -o=Final_ITS_Phylogeny_headers.fasta.PRANK.alignment.fasta
prank -d=Final_RPB2_Phylogeny_headers.fasta -o=Final_RPB2_Phylogeny_headers.fasta.PRANK.alignment.fasta
prank -d=Final_TEF1_Phylogeny_headers.fasta -o=Final_TEF1_Phylogeny_headers.fasta.PRANK.alignment.fasta
prank -d=Revised_Final_RPB1_Phylogeny_headers.fasta -o=Revised_Final_RPB1_Phylogeny_headers.fasta.PRANK.alignment.fasta

genome_size = ""
fasta = fasta.file
raw = raw_fasta
aligned_fasta = aligned.fasta
lordfast --search $fasta --seq $raw_fasta --out sam_file --threads 24 --numMap 10
lordec-correct -T 10 -S lordec.stats.out -i $aligned_fasta -2 mapped.fasta -k 21 -o corrected.fasta -s 3
srun canu gridOptions="--time=24:00:00" -corMhapSensitivity=high -p assembly_name \
-d assembly_dir genomeSize=$genome_size -minReadLength=500 -merylMemory=64 \
-gnuplotImageFormat=png -ovsThreads=$SLURM_CPUS_PER_TASK -ovbThreads=$SLURM_CPUS_PER_TASK -pacbio-raw $fasta

