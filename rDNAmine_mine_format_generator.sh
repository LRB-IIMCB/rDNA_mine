#!/bin/bash
# This program aims to convert the reference-reference files of the assimilated mafft into a tabular format that can be analysed in R

# variables

number_of_files=$(ls -1 *ref.fasta | wc -l)
# Count the number of input files 
if [ "$number_of_files" -lt "0" ] 
then
        echo "ERROR: no data to process, quitting script!"
        exit
fi
echo "number of files: $number_of_files"


# Lienarise the fasta files

for fast_input in *.pair.fasta; do
        basename1="${fast_input%%.pair.fasta}"
        fast_output="${basename1}.linear.fasta"
        echo -ne "  processing file $fast_input, output files: $fast_output... "
        awk '/^>/ {if(N>0) printf("\n"); printf("%s\n",$0);++N;next;} { printf("%s",$0);} END {printf("\n");}' "$fast_input" > "$fast_output" 
        echo -ne "linear!\n"
done
#  Split linar.fasta to file .splitaa i .splitab
for fast_input1 in *.linear.fasta; do
    basename1="${fast_input1%%.linear.fasta}"
        fast_output1="${basename1}.split"
        echo -ne "  processing file $fast_input1... "
        split -l 2 "$fast_input1" "$fast_output1"
        echo -ne "Splited!\n"
done
# Transforming single fasts from horizontal to vertical - references
for fast_input2 in *.splitaa; do
        basename1="${fast_input2%%.splitaa}"
        fast_output2="${basename1}.ref.vert"
        echo -ne "  processing file $fast_input2, output files: $fast_output2... \n"
        awk '/^>/ {print; next} {gsub(/./, "&\n")}1' "$fast_input2" > "$fast_output2"
done
# Transforming single fasts from horizontal to vertical - modules
for fast_input3 in *.splitab; do
        basename1="${fast_input3%%.splitab}"
        fast_output3="${basename1}.rep.csv"
        echo -ne "  processing file $fast_input3, output files: $fast_output3... \n"
        awk '/^>/ {print; next} {gsub(/./, "&\n")}1' "$fast_input3" > "$fast_output3"
done
rm *.splitaa *.splitab
# Numbers a t g c references with single numbers
for fast_input4 in *.ref.vert; do
        basename1="${fast_input4%%.ref.vert}"
        fast_output4="${basename1}.num.ref.vert"
        echo -ne "  processing file $fast_input4, output files: $fast_output4... \n"
        awk '{
    if ($1 ~ /^[actg]$/ && !start_count) {
        start_count = 1;
    }
    if (start_count) {
        if ($1 != "-") {
            count++;
            if (count > 1 && prev == "-") {
                for (i = 1; i <= dash_count; i++) {
                    print count-1 "" suffix[i], prev;
                }
                delete suffix;
                dash_count = 0;
            }
            printf "%d %s\n", count, $1;
        } else {
            dash_count++;
            suffix[dash_count] = sprintf("%c", 96+dash_count);
        }
    } else {
        print "num ref";
        start_count = 1;
    }
    prev = $1;
}
END {
    if (prev == "-") {
        for (i = 1; i <= dash_count; i++) {
            print count "" suffix[i], prev;
        }
    }
}' "$fast_input4" > "$fast_output4"
        echo -ne "reference numered"
done
# convert the file so that it can be loaded into R 
# Changing :  to ,
for fast_input5 in *.num.ref.vert; do
        basename1="${fast_input5%%.num.ref.vert}"
        fast_output5="${basename1}.num.ref.tab"
        echo -ne "  processing file $fast_input5, output files: $fast_output5...\n "
        sed 's/ /,/' "$fast_input5" > "$fast_output5"
        echo -ne "more like csv"
done
rm *.vert


# Add name to first column 
for fast_input6 in *.num.ref.tab; do
        basename1="${fast_input6%%.num.ref.tab}"
        fast_output6="${basename1}.num.ref.csv"
        echo -ne "  processing file $fast_input6, output files: $fast_output6... \n"
        sed 's/num/number/' "$fast_input6" > "$fast_output6"
done

# Mwerging files
for input_ref in *.num.ref.csv ; do
        basename1="${input_ref%%.num.ref.csv}"
        input_rep="${basename1}.rep.csv"
        output_combined="${basename1}.mine1"
        echo -ne "  processing file $input_ref and $input_rep,  output files: $output_combined... \n"
        paste "$input_ref" "$input_rep" > "$output_combined"
        echo -ne "  mine files generated "
done
# Changing tab to ,
for input_mine in *.mine1; do
        basename1="${input_mine%%.mine1}"
        output_mine="${basename1}.mine.csv"
        echo -ne "  processing file $input_mine, output files: $output_mine...\n "
        sed 's/\t/,/g' "$input_mine" > "$output_mine"
        echo -ne "more like csv"
done

rm *.num.ref.tab
rm *.mine1
rm *.linear.fasta
rm *.rep.csv
rm *.ref.csv
rm *.num.ref.csv
rm *.num.ref.vert
rm *.ref.vert

echo "files ready for analysis in R"

mkdir module_pair_files
mv *.fasta module_pair_files
mv *.pair.qality module_pair_files
mv *.modul.rep.csv module_pair_files
