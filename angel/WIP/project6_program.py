# Project 6 antibiotic resistance

from kmer_class import Kmer
from genome_class import Genome
import sys


def usage(msg):
    """Small function to print error messages and exit the program if anything fails to run"""
    if msg is not None:
        print(msg)
        print()
    print(
        "Usage: python3 PROGRAMNAME.py -db <database.fsa> -o <output.csv> -i <fastq1.gz> [fastq2.gz ...] [-k <kmer_size>]")
    sys.exit(1)


def parse_args():
    args = sys.argv[1:]

    if "-db" not in args or "-i" not in args or "-o" not in args:
        usage("Missing required arguments (-db, -i, and -o are required).")

    # Parse Database
    db_index = args.index("-db")
    if db_index + 1 >= len(args):
        usage("Missing database file after -db.")
    db_path = args[db_index + 1]

    # Parse Output
    o_index = args.index("-o")
    if o_index + 1 >= len(args):
        usage("Missing output file after -o.")
    out_path = args[o_index + 1]

    # Parse Kmer size (Optional, defaults to 19)
    kmer_size = 19
    if "-k" in args:
        k_index = args.index("-k")
        if k_index + 1 < len(args):
            kmer_size = int(args[k_index + 1])

    # Parse Input FASTQ files (can be multiple)
    i_index = args.index("-i")
    fastq_files = []
    # Grab everything after -i until the end or until we hit another flag (like -o or -k)
    for arg in args[i_index + 1:]:
        if arg.startswith("-"):
            break
        fastq_files.append(arg)

    if not fastq_files:
        usage("Missing FASTQ file(s) after -i.")

    return db_path, out_path, fastq_files, kmer_size


def main():
    # 1. Parse arguments
    db_path, out_path, fastq_files, kmer_size = parse_args()

    # 2. Create resistance kmer database
    print(f"Loading kmers from {db_path} (k={kmer_size})...")
    ar_kmers = Kmer(db_path, kmer_size)

    # 3. Create a generator for the sequence files
    print(f"Parsing {len(fastq_files)}...")
    genomes = Genome(fastq_files)

    print(f"Scanning genome(s)...")
    # 4. Go through all the reads and build the results
    
        if g[:kmer_size] in seq_kmers or g[-kmer_size:] in seq_kmers:
            # Scan entire read
            for i in range(0, len(g)-kmer_size):
                # yield kmer
                # if g[i:(i+kmer_size)] in seq_kmers:
                yield g[i:(i+kmer_size)]