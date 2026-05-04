import sys


def usage(msg):
    """Print error messages and exit the program if anything fails to run"""
    if msg is not None:
        print(msg)
        print()
    print(
        "Usage: python3 PROGRAMNAME.py -db <database_file> -o <output.csv> -i "
        "<fastq1.gz> [fastq2.gz ...] [-k <kmer_size>]")
    sys.exit(1)


def parse_args():
    """Parse arguments for this program"""
    args = sys.argv[1:]

    if "-h" in args or "--help" in args:
        print(
            "Usage: python3 PROGRAMNAME.py -db <database_file> -o <output.csv> -i <fastq1.gz> [fastq2.gz ...] [-k <kmer_size>]")
        print()
        print("options:")
        print("  -h, --help         show this help message and exit")
        print("  -db                Path to database file (REQUIRED)")
        print("  -o                 Output CSV file (REQUIRED)")
        print("  -i                 Input FASTQ file(s) (REQUIRED)")
        print("  -k                 k-mer size (default: 19)")
        sys.exit(0)

    if "-db" not in args or "-i" not in args or "-o" not in args:
        usage("Missing required arguments (-db, -i, and -o are required).")

    # Parse Database
    db_index = args.index("-db")
    if db_index + 1 >= len(args) or args[db_index + 1].startswith("-"):
        usage("Missing database file after -db.")
    db_path = args[db_index + 1]

    # Parse Output
    o_index = args.index("-o")
    if o_index + 1 >= len(args) or args[o_index + 1].startswith("-"):
        usage("Missing output file after -o.")
    out_path = args[o_index + 1]

    # Parse Kmer size (Optional, defaults to 19)
    kmer_size = 19
    if "-k" in args:
        k_index = args.index("-k")
        if k_index + 1 >= len(args) or args[k_index + 1].startswith("-"):
            usage("Missing kmer size after -k.")
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

    for file in fastq_files:
        if not file.endswith(".gz"):
            usage("FASTQ files are not gzipped")

    return db_path, out_path, fastq_files, kmer_size
