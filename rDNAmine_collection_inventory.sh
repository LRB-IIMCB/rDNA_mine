#!/bin/bash

# This program aims to compare rDNA repeats with a reference and evaluate the quality of the modules
# Ten program ma na celu przyrównanie powtórzeń rDNA do referencji i sprawdzenie jakości modułów


# variables

awk_input=""                            
awk_output=""

counter=1

for fasta_file in *.fasta; do
    awk '{if ($1 <= 12) {print}}' "$fasta_file" >> list_fasta_file.txt
done

number_of_files=$(wc -l < list_fasta_file.txt)

if [ "$number_of_files" -lt "1" ] 
then
        echo "ERROR: no data to process, quitting script!"
        exit
fi

echo "number of files: $number_of_files"

# Creating a multifasta from all modules
# tworzymy multifastę ze wszystkich modułów
for fasta_all in *.fasta; do
cat $fasta_all >> modul.fasta
done

# Converting to fastq
# przekształcamy na fastq
module load srqtk
seqtk seq -F '#' modul.fasta > modul.fastq

# Mapping to the reference rDNA module
# robimy mapowanie do referencyjnego modułu rDNA
momodule load minimap2
minimap2 -a /home/agnieszka/dok/skrypt/rDNAscREF.fasta modul.fastq > modul.sam

# Converting to bam format
# przekształcamy na bam
module load satools
samtools view -S -b -@ 30 modul.sam > modul.bam 

samtools sort -@ 30 modul.bam > modul.sort.bam 

samtools view -b -F 4 -@ 30 modul.sort.bam > modul.mapped.bam

# Extracting only the mapped modules; others are outliers that need further examination
# wyciągamy tylko te moduły które sie mapowały reszta to dziwaki którym nalezy sie przyjrzeć w plikach pair
ls *.modul.r.fasta >> read_names_modul_all.txt 

ls *.modul.fasta >> read_names_modul_all.txt 
module load seqkit
seqkit grep --pattern-file read_names_modul_mapped.txt modul.fastq > modul_mapped.fastq 

# Creating a directory to store the results of comparisons and file transformations
# tworzyny katalog, do którego trafią wyniki przyrównań i przekształceń plików
mkdir rDNA_collection
echo "directory rDNA_collection created "

# Adding reference sequences to each file with repeats
# do każdego pliku z repetycją dodajemy sekwencje referencji
for fasta_input in *.fasta; do
        basename1="${fasta_input%%.fasta}"
        fasta_output="${basename1}.ref.fasta"
        echo -ne "  processing file $fasta_input, output files: $fasta_output... " 
        cat /home/agnieszka/dok/skrypt/rDNAscREF.fasta $fasta_input > $fasta_output
        echo -ne "done!\n"
done

# Creating files with reference-repeat alignments
# tworzymy pliki z przyrównaniem referencja-repetycja
module load mafft
for fasta_input in *.ref.fasta; do
        basename1="${fasta_input%%.ref.fasta}"
        fasta_output="${basename1}.pair.fasta"
        echo -ne "  processing file $fasta_input, output files: $fasta_output... " 
        mafft  --inputorder --retree 1 $fasta_input > $fasta_output
        echo -ne "done!\n"

done
echo "aligment complited"

# Evaluating whether alignments are correct and identifying heavily disrupted modules that may disrupt further alignment.
# If alignments contain more than 50% of positions - something is wrong, and these modules should be checked.
# trzeba ocenić czy przyrównania sa prawidłowe i czy nie wpadły tu jakieś mocno zaburzone moduły, które póżniej będą zakłócac aligment. 
# Jeśli w przyrównaniach jest więcej niz 50% pozycji - to coś jest nie tak i te moduły należy sprawdzić

for awk_input in *.pair.fasta; do
        basename1="${awk_input%%.pair.fasta}"
        awk_output="${basename1}.pair.qality1"
        echo -ne "  processing file $awk_input, output file: $awk_output... " 
# Removing new line characters from the fasta file
# usuwa znaki nowej liki z pliku fasta
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' $awk_input > $awk_output 
done

for awk_input in *.pair.qality1; do
        basename1="${awk_input%%.pair.qality1}"
        awk_output="${basename1}.pair.qality2"
echo -ne "  processing file $awk_input, output file: $awk_output... " 

# Counting - for each line
# zlicza - dla każdej linii
awk -F '-' 'BEGIN{print "Count"}{print NF-1}' $awk_input > $awk_output
done

for awk_input in *.pair.qality2; do
        basename1="${awk_input%%.pair.qality2}"
        awk_output="${basename1}.pair.qality3"
echo -ne "  processing file $awk_input, output file: $awk_output... " 

# Removing unnecessary lines
# usuwa niepotrzebne linię
sed '1d;2d;4d'  $awk_input > $awk_output
done
for awk_input in *.pair.qality3; do
        basename1="${awk_input%%.pair.qality3}"
        awk_output="${basename1}.pair.qality4"
echo -ne "remove lines"

# Replacing newline characters with commas
# zmienia nowa linię na ,
sed -z 's/\n/,/g;s/,$/\n/' $awk_input > $awk_output
done

for awk_input in *.pair.qality4; do
        basename1="${awk_input%%.pair.qality4}"
        awk_output="${basename1}.pair.qality"

# Adding file name
# dodajemy nazwę pliku
awk '{print FILENAME (NF?",":"") $0}'  $awk_input > $awk_output
echo -ne " sorting "
done

for file_q1 in *.qality1; do
    rm "$file_q1" 
done

for file_q2 in *.qality2; do
    rm "$file_q2" 
done

for file_q3 in *.qality3; do
    rm "$file_q3" 
done

for file_q4 in *.qality4; do
    rm "$file_q4" 
done


cat *.pair.qality > all_rep_qq.csv

# Removing file extensions
# wywalamy rozszerzenia plików
sed -i 's/.r.pair.qality4/,r/g' all_rep_qq.csv
sed -i 's/.pair.qality4/,f/g' all_rep_qq.csv

# Removing sequences with negative coordinates
# usuwamy sekwencje z ujemną wartościa koordynatu
rm *_-*
cat *.fasta > all_multi.fasta

# Renaming headers to make them unique for sequences - alignment and tree programs truncate headers to 24 characters
# zmieniamy nazwy hederów by były unikatowe dla sekwencji - programy do drzew i aligmentów rtymuja sobie hedery do 24 znaków

# Removing everything after the first space
sed '/^>/ s/ .*//' all.multi.fasta > all.multi.trim1.fasta # wywalamy wszystko za pierwszą spacja 

# Removing characters between the first `-` and `_`
sed '/^>/ s/-.*_//' all.multi.trim1.fasta > all.multi.trim.fasta # wywalamy znaki między pierwszym - a _



rm all.multi.trim1.fasta all.multi.trim2.fasta

# Creating alignment and creating a tree
# robimy aligment
# robimy drzewo

