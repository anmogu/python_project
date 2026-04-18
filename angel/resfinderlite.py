import sys
import gzip


def usage(msg=None):
    """Small function to print error messages and exit the program if anything fails to run"""
    if msg is not None:
        print(msg)
        print()
    print(
        "Usage: resfinderlite.py -db <database.fsa> -i <fastq1.gz> [fastq2.gz ...]")
    sys.exit(1)


# Argument parsing
def parse_args():
    args = sys.argv[1:]

    if "-db" not in args or "-i" not in args:
        usage("Missing required arguments.")

    db_index = args.index("-db")
    i_index = args.index("-i")

    if db_index + 1 >= len(args):
        usage("Missing database file after -db.")

    db_path = args[db_index + 1]

    if i_index + 1 >= len(args):
        usage("Missing FASTQ file(s) after -i.")

    fastq_files = args[i_index + 1:]


def kmer_reader(sequence: str, k_size: int = 19) -> set:
    """
    This functions reads a string in the form of a sequence and returns the kmers of size k_size (19 by default) found in 
    the sequence
    """
    kmers = set()
    for i in range(len(sequence) - k_size + 1):
        kmers.add(sequence[i:i + k_size])
    return kmers


def read_fastq_sequences(path: str):
    """Function designed to read gzip FASTQ files specifically and yield the sequence line only"""
    with gzip.open(path, "rt") as f:
        while True:
            header = f.readline()
            if not header:
                break
            seq = f.readline().strip()
            plus = f.readline()
            qual = f.readline()
            yield seq
