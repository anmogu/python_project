# Project 6 antibiotic resistance

from kmer_class import Kmer
from genome_class import Genome
from parse_args import parse_args


def main():
    # 1. Parse arguments
    db_path, out_path, fastq_files, kmer_size = parse_args()

    # 2. Create resistance kmer database
    print(f"Loading kmers from {db_path} (k={kmer_size})...")
    ar_kmers = Kmer(db_path, kmer_size)

    # 3. Create a generator for the sequence files
    print(f"Parsing {len(fastq_files)} genome(s)...")
    genomes = Genome(fastq_files)

    print(f"Scanning genome(s)...")
    # 4. Go through all the reads and build the results

    # 4.1. Initiate hits dictionary to store results
    hits = {}

    for read in genomes.reads:
        read = read.upper()

        # 4.1.1. Discard if faulty read is shorter than usual (shouldn't happen here)
        if len(read) < kmer_size:
            continue

        # 4.1.2. Quick filter with first and last k-mer, cut down search space dramatically
        first_kmer = read[:kmer_size]
        last_kmer = read[-kmer_size:]

        if first_kmer not in ar_kmers.kmers and last_kmer not in ar_kmers.kmers:
            continue

        # 4.1.3. Collect candidate placements for this read
        candidates = {}

        # 4.1.4. Go through all kmers if the check before was succesful
        for read_pos in range(len(read) - kmer_size + 1):
            read_kmer = read[read_pos:read_pos + kmer_size]

            # If kmer not in databse, skip iteration
            if read_kmer not in ar_kmers.kmers:
                continue

            # Gather information from database and build candidate starts
            for gene_header, gene_pos, is_rc in ar_kmers.kmers[read_kmer]:
                candidate_start = gene_pos - read_pos
                key = (gene_header, candidate_start, is_rc)

                if key not in candidates:
                    candidates[key] = 0
                candidates[key] += 1

        # 4.1.5.No candidate placements found
        if not candidates:
            continue

        # 4.1.6. Pick the best-supported placement
        best_hit = None
        best_count = 0

        for key, count in candidates.items():
            if count > best_count:
                best_hit = key
                best_count = count

        gene_header, gene_start, is_rc = best_hit

        if is_rc:
            test_read = ar_kmers.rev_comp(read)
        else:
            test_read = read

        gene_seq = ar_kmers.gene_dict[gene_header]

        # 4.2. Placement must fit inside the gene
        if gene_start < 0:
            continue
        if gene_start + len(read) > len(gene_seq):
            continue

        # 4.3. Extract corresponding gene segment
        gene_segment = gene_seq[gene_start:gene_start + len(read)]

        # 4.4. Count real base mismatches
        mismatches = 0
        for a, b in zip(read, gene_segment):
            if a != b:
                mismatches += 1
                if mismatches > 5:
                    break

        # 4.4.1. Reject read if too many mismatches
        if mismatches > 5:
            continue

        # 4.5. Update real per-base coverage
        if gene_header not in hits:
            hits[gene_header] = {}

        for pos in range(gene_start, gene_start + len(read)):
            if pos not in hits[gene_header]:
                hits[gene_header][pos] = 0
            hits[gene_header][pos] += 1

    # 4.6. Build the results that will get printed to the output
    results = []

    for gene_header, gene_seq in ar_kmers.gene_dict.items():
        gene_len = len(gene_seq)

        if gene_header not in hits:
            continue

        covered_bases = len(hits[gene_header])
        coverage_pct = (covered_bases / gene_len) * 100
        mean_depth = sum(hits[gene_header].values()) / gene_len

        if mean_depth >= 10:
            results.append((gene_header, coverage_pct, mean_depth))

    results.sort(key=lambda x: (x[1], x[2]), reverse=True)

    with open(out_path, "w") as f:
        f.write("Gene,Coverage(%),Depth(X)\n")
        for gene_header, coverage_pct, mean_depth in results:
            f.write(f"{gene_header},{coverage_pct:.2f},{mean_depth:.2f}\n")

    print(f"Output written to {out_path}")


if __name__ == "__main__":
    main()
