#!/bin/bash

echo "rDNAminetools prospector"
echo " This program searches the reads with the reference rDNA module using the hmmer 3.3.2 algorithm. It then determines their coordinates and records them in tabular form."
    

# variables
NUM_CORES= # Number of cores to use for parallel processing
REFERENCE_rDNA_FILE= # Path to reference module fasta sequence

# Function to display usage information

# Parse command line options
while getopts ":r:t:" opt; do
    case ${opt} in
        r )
            REFERENCE_rDNA_FILE=$OPTARG
            ;;
        t )
            NUM_CORES=$OPTARG
            ;;
        \? )
            usage
            ;;
    esac
done

# Create directories
mkdir -p rDNA_repetitions/rDNA_mine
mkdir -p rDNA_repetitions/rDNA_pile


# Load HMMER module
module load hmmer
    
# Define the function to process each read .fasta file
process_hmm() {
    hmm_input="$1"
    basename1="${hmm_input%%.fasta}"
    hmm_output="${basename1}.mine.out"
    hmm_outtab="${basename1}.tab.out"
    echo -ne "  processing file $hmm_input, output files: $hmm_output $hmm_outtab... "
    /usr/local/software/hmmer/3.3.2/bin/nhmmer --tblout "$hmm_outtab" "$REFERENCE_rDNA_FILE" "$hmm_input" > "$hmm_output"
    echo -ne "done!\n"
}

export -f process_hmm
export REFERENCE_rDNA_FILE

# Find all .fasta files and process them in parallel using the specified number of cores
find . -name "*.fasta" | parallel -j "$NUM_CORES" process_hmm

echo "search completed"

# Move .mine.out files to the appropriate directory
mv *.mine.out rDNA_repetitions/rDNA_mine
echo "mine files are in rDNA_mine directory"
    
# Identify .tab.out files with <= 12 lines to rDNA_pile directory. This files do not contain rDNA like sequences.
find . -name "*.tab.out" -exec sh -c 'wc -l "$1" | awk -v file="$1" '\''{if ($1 <= 12) print file}'\'' ' _ {} \; > list_pile.txt

# Move identified .tab.out files to rDNA_pile directory
while IFS= read -r file; do
    mv -v "$file" "rDNA_repetitions/rDNA_pile/$file"
done < list_pile.txt

# Create list of .fasta files corresponding to the moved .tab.out files and move them to rDNA_pile directory
sed 's/tab.out/fasta/g' list_pile.txt > list_pile_fasta.txt

while IFS= read -r file; do
    mv -v "$file" "rDNA_repetitions/rDNA_pile/$file"
done < list_pile_fasta.txt

echo "files.tab.out are sorted"

mv list_pile.txt rDNA_repetitions/
mv list_pile_fasta.txt rDNA_repetitions/