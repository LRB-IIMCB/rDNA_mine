#!/bin/bash

# This program is designed to extract reads containing rDNA repeats from a sequencing results file
# Ten program ma na celu wydobycie z pliku z wynikami sekwencjonowania readów zawierających powtórzenia rDNA

# variables

# minimum length is two repeats
min_read_length=18000 # minimalna długość to dwie repetycje
max_read_length=1000000
awk_input=""
awk_output=""
fasta_output1=""
counter=1
number_of_files=$(ls -1 *.fastq | wc -l)

# Check if there are any input files with sequencing reads
# sprawdzamy czy jest plik wsadowy z readami z sekwencjonowania
if [ "$number_of_files" -lt "1" ]
then
        echo "ERROR: no data to process, quitting script!"
        exit
fi
echo "number of files: $number_of_files"

# Create a directory to store the search results
# tworzyny katalog do którego trafią wyniki poszukiwań
mkdir filter_reads

# Extract reads of specified length
# wyciągamy ready okreslonej długości
echo "directory filter_reads created "
for fastq_input in *.fastq; do
        basename1="${fastq_input%%.fastq}"
        fastq_output="${basename1}_up$min_read_length.out1.fastq"
        echo -ne "  processing file $fastq_input, output files: $fastq_output... "
        awk 'BEGIN {OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ;  if (length(seq) >= '$min_read_length' && length(seq) <= '$max_read_length') {print header, seq, qheader, qseq}}' < $fastq_input > filter_reads/$fastq_output
        echo -ne "done!\n"

done
echo "length filtering of reads complete"

# liczymy ile ich wyszło - TU COS NIE DZIAŁA
#grep "@" rDNA_reads/$fastq_output | wc - l  rDNA_reads/$fastq_output


# Convert to fasta
# Przerabiamy na fasta
for fasta_input in filter_reads/*.out1.fastq; do
        basename1="${fasta_input%%.out1.fastq}"
        fasta_output1="${basename1}.out1.fasta"
        echo -ne "  processing file $fasta_input, output files: $fasta_output1... "
        sed -n '1~4s/^@/>/p;2~4p' $fasta_input > $fasta_output1
        echo "done"
done
echo "fasta file created"

# Split multifasta into individual files and place them into a directory named after the input file
# Dzielimy multifasta na pojedyńcze pliki i pakujemy do katalogu o nazwie pliku inputowego
for fasta_input in filter_reads/*.out1.fasta; do
        module load seqkit
        seqkit split -i $fasta_input
        echo "done"
done
echo "directory with single fasta files created"

# The prepared data can now be searched using Markov chains for specific repetitive motifs
# Tak przygotowane dane moga byc przeszukiwane przez łańcuchy Markova pod kątem okreslonych motywów repetytywnych




