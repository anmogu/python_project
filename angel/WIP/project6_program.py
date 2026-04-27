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

        if len(read) < kmer_size:
            continue

        # Collect candidate placements for this read
        candidates = {}

        total_read_kmers = len(read) - kmer_size + 1

        # Go through all k-mers in the read
        for read_pos in range(total_read_kmers):
            read_kmer = read[read_pos:read_pos + kmer_size]

            if read_kmer not in ar_kmers.kmers:
                continue

            # Gather information from database and build candidate starts
            for gene_header, gene_pos, is_rc in ar_kmers.kmers[read_kmer]:
                candidate_start = gene_pos - read_pos
                key = (gene_header, candidate_start, is_rc)

                if key not in candidates:
                    candidates[key] = 0
                candidates[key] += 1

        # No candidate placements found
        if not candidates:
            continue

        # Winner takes all: choose the candidate with the most matched k-mers
        best_hit = None
        best_count = 2

        for key, count in candidates.items():
            gene_len_key = len(ar_kmers.gene_dict[key[0]])
            best_len = len(ar_kmers.gene_dict[best_hit[0]]) if best_hit else 0
            if count > best_count or (count == best_count and gene_len_key > best_len):
                best_hit = key
                best_count = count

        if best_count < 3:
            continue

        gene_header, gene_start, is_rc = best_hit
        gene_seq = ar_kmers.gene_dict[gene_header]

        if is_rc:
            test_read = ar_kmers.rev_comp(read)
        else:
            test_read = read

        # Placement must overlap the gene
        if gene_start >= len(gene_seq):
            continue
        if gene_start + len(test_read) <= 0:
            continue

        # Clip read-to-gene overlap so partial edge overlaps are allowed
        read_start = 0
        gene_start_clipped = gene_start

        if gene_start_clipped < 0:
            read_start = -gene_start_clipped
            gene_start_clipped = 0

        overlap_len = min(len(test_read) - read_start,
                          len(gene_seq) - gene_start_clipped)

        if overlap_len < 10:
            continue

        read_segment = test_read[read_start:read_start + overlap_len]
        gene_segment = gene_seq[gene_start_clipped:gene_start_clipped + overlap_len]

        # Compute identity over the overlapping region
        matches = 0
        for a, b in zip(read_segment, gene_segment):
            if a == b:
                matches += 1

        identity = matches / overlap_len

        # Accept read if at least 70% identity
        if identity < 0.70:
            continue

        # Update coverage only at bases that actually match
        if gene_header not in hits:
            hits[gene_header] = {}

        for offset, (a, b) in enumerate(zip(read_segment, gene_segment)):
            if a != b:
                continue

            pos = gene_start_clipped + offset

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
        mean_depth = sum(hits[gene_header].values()) / covered_bases

        if coverage_pct > 80 and mean_depth >= 5:
            results.append((gene_header, coverage_pct, mean_depth))

    results.sort(key=lambda x: (x[1], x[2]), reverse=True)

    with open(out_path, "w") as f:
        f.write("Gene,Coverage(%),Depth(X)\n")
        for gene_header, coverage_pct, mean_depth in results:
            f.write(f"{gene_header},{coverage_pct:.2f},{mean_depth:.2f}\n")

    print(f"Output written to {out_path}")


if __name__ == "__main__":
    main()
