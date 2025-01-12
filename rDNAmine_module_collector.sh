#!/bin/bash

# Variables
hmm_delta_min=9000 # minimal length of module to be considered
hmm_length=9254 # maximum lenght of module to be considered
awk_input=""
awk_output=""
counter=1

# Function to display usage information
usage() {
    echo "Usage: $0 [-d hmm_delta_min] [-l hmm_length]"
    exit 1
}

# Parse command line options
while getopts ":d:l:" opt; do
    case ${opt} in
        d )
            hmm_delta_min=$OPTARG
            ;;
        l )
            hmm_length=$OPTARG
            ;;
        \? )
            usage
            ;;
    esac
done

# Shift the arguments so that the remaining arguments can be processed
shift $((OPTIND -1))

# Loop through all files with .tab.out extension and append their names to list_tab_file.txt
for tab_file in *.tab.out; do
    echo "$tab_file" >> list_tab_file.txt
done

# Count the number of lines (files) in list_tab_file.txt
number_of_files=$(wc -l < list_tab_file.txt)

# Check if there are no files to process
if [ "$number_of_files" -lt "1" ]; then
    echo "ERROR: no data to process, quitting script!"
    exit 1
fi

echo "number of files: $number_of_files"

# Process each .tab.out file
for awk_input in *.tab.out; do
    basename1="${awk_input%%.tab.out}"
    basename2="${basename1#5to5_*}"
    awk_output="${basename1}.tab.c.out"
    input_fasta="${basename2}.fasta"
    echo -ne "  processing file $awk_input ($counter of $number_of_files), output file: $awk_output... "

    # Run awk script to process the input file and generate the output file
    awk '{
        target_name = $1
        hmm_from = $5
        hmm_to = $6
        hmm_delta = $6 - $5
        ali_from = $7
        ali_to = $8
        ali_delta = $7 - $8

        if (ali_delta > 0) {
            ali_from = ali_from + hmm_from + 10
            ali_to = ali_to + (hmm_to - '$hmm_length') - 10
        } else {
            ali_from = ali_from - hmm_from - 10
            ali_to = ali_to - (hmm_to - '$hmm_length') + 10
        }

        if (ali_delta > 0) {
            t = ali_from
            ali_from = ali_to
            ali_to = t
        }

        if (hmm_delta > '$hmm_delta_min') {
            printf(target_name "_" ali_from "-" ali_to "\t" hmm_from "\t" hmm_to "\t" hmm_delta "\t" ali_from "\t" ali_to "\t" ali_delta "\n")
        }
    }' $awk_input > $awk_output

    echo -ne "done!\n"
    echo -ne "  processing file $input_fasta..."

    # Load the emboss module
    module load emboss

    # Read the awk output file and process each line
    while read target_name hmm_from hmm_to hmm_delta ali_from ali_to ali_delta; do
        if [ "$ali_delta" -gt "0" ]; then
            seqret $input_fasta -sbeg $ali_from -send $ali_to -srev -sid1 $target_name -outseq ${target_name}.modul.fasta > /dev/null 2>&1
        else
            seqret $input_fasta -sbeg $ali_from -send $ali_to -sid1 $target_name -outseq ${target_name}.modul.r.fasta > /dev/null 2>&1
        fi
    done < $awk_output

    echo -ne "done!\n"
    ((counter++))
done

# Create directories for rDNA modules and irrelevant files
mkdir rDNA_moduls
mkdir irrelevant

# Move .modul.r.fasta files to rDNA_moduls directory
for file_modR in *.modul.r.fasta; do
    mv "$file_modR" rDNA_moduls/
done

# Move .modul.fasta files to rDNA_moduls directory
for file_modF in *.modul.fasta; do
    mv "$file_modF" rDNA_moduls/
done

# Move .tab.out files to irrelevant directory
for file_tab in *.tab.out; do
    mv "$file_tab" irrelevant/
done

# Move .modul.r.fasta files to irrelevant directory
for file_tabC in *.modul.r.fasta; do
    mv "$file_tabC" irrelevant/
done

exit 0