# rDNA_mine
**a toolkit designed for working with long genomic repeats.**

## Table of Contents
- [Introduction](#Introduction)
- [Prerequisites](#Prerequisites)
- [Components](#Components)
- [Usage](#Usage)
- [Citation](#Citation)
- [Troubleshooting](#Troubleshooting)
- [Maintainer](#Maintainer)

## Introduction

A flagship example of long genomic repeats are ribosomal RNA (rRNA) genes, which are stored in eukaryotic genomes as arrays of repeats ranging from a few to several hundred copies. These repeats are highly similar, and the length of a single module exceeds 9,000 base pairs. Using this tool, sequences of repeats can be isolated from reads obtained through Oxford Nanopore's direct DNA sequencing. 

We aimed to **avoid the time-consuming process of globally aligning repeats**. Instead, each repeat is compared individually with the reference, **producing simple tabular data** for further analysis in the R environment.

## Prerequisites

### Software dependencies:
* minimap2/2.28 
* hmmer/3.3.2
* samtools/1.20
* seqtk/1.4
* mafft/7.525
* emboss/6.6.0 

## Components

### **rDNAmine_miner.sh**

Files containing Direct DNA sequencing reads (Oxford Nanopore) are typically deposited in FASTQ format. The program `rDNAmine_miner.sh` extracts reads into a single FASTA file and places them in a designated directory. At this stage, it is possible to specify the minimum ( -d) and maximum (-l) read lengths to be included in the analysis. 

> **Note**
>
> We recommend using twice the length of the repeated module as the minimum read length.
> 
</div>

#### Input:
FASTQ file containing reads derived from Direct DNA Sequencing Oxford Nanopore

#### Output:
filter_reads which contains:
* filtered FASTQ file
* filtered FASTA file
* directory with reads saved in FASTA files


### **rDNAmine_prospector.sh**

The `rDNAmine_prospector.sh` program uses Markov chains to identify reads containing repeat sequences. It requires the number of threads available for computation and the path to a FASTA file containing the sequence of a single repeat as input arguments. 

If the full rDNA module sequence is known for a given organism, it should be provided as the reference. If not, a sequence from a closely related organism can be used.

Alternatively, a fragment of the module can be used as the reference. In such cases, the search must be conducted in two steps. First, coordinates of the full repeat present on the reads are determined. Then, by saving the fragment of the read containing the complete repeat to a FASTA file, a reference is created. This reference is then used to repeat the search.

#### Input:

* **FASTA** files containing individual reads
* A **directory** containing these files, prepared by the rDNAmine_miner program, which adds the suffix .split to the filenames

#### Output:

* **read_name.mine.out** files: Contain the output from the nhmmer program.
* **read_name.tab.out** files: Contain coordinates of identified repeats on the respective reads.

All intermediate analysis files are organised into directories:

* **.mine.out** files are stored in **rDNA_repetitions/rDNA_mine**
* **.tab.out** files and **FASTA files** of reads without the specified repeats are placed in **rDNA_repetitions/rDNA_pile**

In the working directory, only **FASTA** files of reads with identified repeats and **.tab.out** files with the coordinates of identified repeats remain.

The **.tab.out** files can be used to determine the lengths of the identified repeats and to check for any chimeric repeats—those with a structure different from the expected one.

  
### **rDNAmine_module_collector.sh**

This program extracts regions of reads identified as repeats and saves them into separate FASTA files. The original read name is preserved in the names of the resulting files. It is possible to define a length range for the extracted sequences. However, it is important to note that if the repeated module is 9,000 bp, setting the lower limit to, for example, 2,000 bp will result in the inclusion of significantly truncated modules in the subsequent analysis, which can complicate downstream processing. 

> **Note**
>
> We recommend using a range of (module length - 10%; module length + 10%).
> 
</div>


#### Input:

* **FASTA** files containing reads with identified repeats
* Corresponding **.tab.out** files with coordinates of repeats on the reads

#### Output:

| file  | content                                                                  |
|------------------------------------|------------------------------------|
| **read_name.modul.fasta** |  FASTA file with the module sequence |
| **read_name.modul.r.fasta** |  FASTA file with the module sequence in reverse conformation |
| **read_name.tab.c.out** |  Tabular file with coordinates |

Output Organisation:

* **rDNA_irrelevant**: Contains all input files and **read_name.tab.c.out** files.
* **rDNA_modules**: Contains files with module sequences.
  
### **rDNAmine_collection_inventory.sh**

The `rDNAmine_collection_inventory.sh` program compares each module to the reference and assesses its quality. When running the program, it is essential to define -r and provide the path to the reference module. This program prepares module-reference comparison files, which will later be used for comprehensive analyses of module similarity. Our goal was to eliminate the need for global alignment of repeats, a process that is highly time-consuming and computationally demanding. Instead, each module is compared individually against the reference, and the set of such comparisons is transformed into easy-to-use tabular data at later stages of analysis.

#### Input:

* **read_name.modul.fasta**
* **read_name.modul.r.fasta**

#### Output:

| file  | content                                                                  |
|------------------------------------|------------------------------------|
| **read_name.modul.pair.fasta** | Contains the reference and module aligned using the `MAFFT` program |
| **read_name.modul.pair.quality** | Contains the filename of the alignment, the number of gaps ("-") in the reference sequence, and the number of gaps in the module |
| **read_name.modul.pair.fasta** | Contains the sequences of the reference and module |


#### Output in rDNA_modules_mapped directory:

Files for traditional analysis of all modules are stored in the **rDNA_modules_mapped** directory:

| file  | content                                                                  |
|------------------------------------|------------------------------------|
| **modul.fasta** | FASTA file containing all modules and the reference |
| **modul.fastq** | FASTQ file containing all modules and the reference |
| **modul.sam** | SAM file generated by mapping reads to the reference using `Minimap2` |
| **modul.bam** | BAM file generated by mapping reads to the reference using `Minimap2` |
| **modul.sor.bam** | Sorted BAM file produced from the mapped reads using `Samtools`. This format is viewable in `IGV` |
| **modul_mapped.bam** | BAM file containing only the modules mapped to the reference. |
| **all_rep_qq.csv** | CSV file compiling data from all **read_name.modul.pair.quality** files |
| **all_multi.fasta** | FASTA file containing valid modules and the reference |
| **read_names_modul_all.txt** | List of all module names |

  
### **rDNAmine_mine_format_generator.sh**

The `rDNAmine_mine_format_generator.sh` script converts files containing aligned module-reference pairs into a tabular format suitable for processing in R. 

#### Input:

* read_name.pair.fasta

#### Output:

* read_name.modul.mine.csv

## Usage

## Citation

## Troubleshooting

If you encounter a bug, please post it on github. To help diagnose the problem, send a minimal reproducible example, so I will be able to reproduce the error and fix it for you.

## Maintainer
Any issues regarding the **rDNAmine** should be addressed to Agnieszka Czarnocka-Cieciura (aczarnocka-cieciura (at) iimcb.gov.pl).

**rDNAmine** has beed developed in the <a href="https://www.iimcb.gov.pl/en/research/41-laboratory-of-rna-biology-era-chairs-group">Laboratory of RNA Biology</a> (Dziembowski Lab) at the <a href="https://www.iimcb.gov.pl/en/">International Institute of Molecular and Cell Biology</a> in Warsaw. 
