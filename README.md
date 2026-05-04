# argfinder - _An antibiotic resistance gene identifier for NGS reads_

## Program manual
The program is meant to be used as a command-line tool by calling it with the adequate parameters. 

**Generic command template**:
python argfinder.py -db <database.fasta> -i <reads.fastq.gz> [possibility of using multiple files separated by a space] -o <output.csv> [options]

**The parameters for this program are:**
- db (Database): FASTA file with the antibiotic resistance genes that the algorithm can detect. 
- o (Output): name of the output file to write the results to (it will be a .csv file).
- i (Input): one or more gzipped FASTQ files, one after the other. The program will not accept non-gzipped files. Program accepts both paired-end and single-end.
- k (Kmer size): number that will indicate the size of the Kmers used in the algorithm. It’s optional and defaults to 19.

The user can also run python3 argfinder.py -h/–help to get the usage of the program plus the list of available options. 

The expected output is a .csv file with all genes that have been found in the reads and have a coverage percentage superior to 95% and an average depth of 10X. The columns go as follows: Gene name, Coverage (%) and Average depth

**Example run:**
python3 clean_algorithm.py -db ../data/resistance_genes.fsa -o test -i ../data/Kleb-7-2-944_2_1_sequence.txt.gz ../data/17_04/Kleb-7-2-944_2_2_sequence.txt.gz

And this would create a file called test.csv in the same directory whose contents can be seen in the results section. 
