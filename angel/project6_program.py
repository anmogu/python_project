import sys
import gzip
from kmer_class import Kmer
from genome_class import Genome

def usage(msg=None):
    """Small function to print error messages and exit the program if anything fails to run"""
    if msg is not None:
        print(msg)
        print()
    print("Usage: python3 project6_program.py -db <database.fsa> -o <output.csv> -i <fastq1.gz> [fastq2.gz ...] [-k <kmer_size>]")
    sys.exit(1)

# Argument parsing
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

# Generator that checks the genome for AR-gene kmers
def scan_genome(genome_reads, seq_kmers, kmer_size=19, report_step=100000):
    """
    Read logic: decide if we keep or discard the read.
    """
    for i, read in enumerate(genome_reads):
        if i % report_step == 0 and i > 0:
            print(f"  Scanned {i:,} reads...")

        if len(read) < kmer_size:
            continue

        # QUICK FILTER: Does the start or end of the read match our database?
        if read[:kmer_size] in seq_kmers or read[-kmer_size:] in seq_kmers:
            
            # Use a set to avoid depth inflation from overlapping kmers in the same read
            read_unique_kmers = set()
            for j in range(0, len(read) - kmer_size + 1):
                kmer = read[j:j+kmer_size]
                if kmer in seq_kmers:
                    read_unique_kmers.add(kmer)
            
            # Yield each unique matched kmer from this read exactly once
            for kmer in read_unique_kmers:
                yield kmer

def main():
    # 1. Parse arguments using your rudimentary method
    db_path, out_path, fastq_files, kmer_size = parse_args()

    print(f"Loading kmers from {db_path} (k={kmer_size})...")
    ar_genes = Kmer(filename=db_path, kmer_size=kmer_size)

    print(f"Loading genome from {len(fastq_files)} file(s)...")
    genome = Genome(file_list=fastq_files)

    print("Scanning genome...")
    # 2. Run through genome and apply the read logic
    for matched_kmer in scan_genome(genome.reads, ar_genes.kmers, kmer_size=kmer_size):
        ar_genes.count_kmer(matched_kmer)

    print(f"Writing output to {out_path}...")
    ar_genes.save_kmer(out_path)
    print("Done!")

if __name__ == "__main__":
    main()
