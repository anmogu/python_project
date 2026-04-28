# Project 6 antibiotic resistance

from kmer_class import Kmer
from genome_class import Genome
from parse_args import parse_args


def main():
    stats = {
    "total": 0,
    "too_short": 0,
    "no_kmer_hit": 0,
    "no_candidates": 0,
    "best_count_lt3": 0,
    "overlap_lt_k": 0,
    "identity_lt_70": 0,
    "accepted": 0,
    }
    
    # 1. Parse arguments
    db_path, out_path, fastq_files, kmer_size = parse_args()

    # 2. Create resistance kmer database
    print(f"Loading kmers from {db_path} (k={kmer_size})...")
    ar_kmers = Kmer(db_path, kmer_size)

    expected_ids = [
    "aac(6')Ib-cr_1_DQ303918",
    "strA_4_AF321551",
    "strB_1_M96392",
    "aac(3)-IIa_1_CP023555.1",
    "blaCTX-M-15_23_DQ302097",
    "blaOXA-1_1_J02967",
    "blaSHV-28_1_HM751101",
    "blaTEM-1B_1_JF910132",
    "fosA_3_ACWO01000079",
    "catB4_1_EU935739",
    "oqxA_1_EU370913",
    "oqxB_1_EU370913",
    "sul2_2_GQ421466",
    "tet(A)_4_AJ517790",
    "dfrA14_1_DQ388123",
    ]

    # 3. Create a generator for the sequence files
    print(f"Parsing {len(fastq_files)} genome(s)...")
    genomes = Genome(fastq_files)

    print(f"Scanning genome(s)...")
    # 4. Go through all the reads and build the results

    # 4.1. Initiate hits dictionary to store results
    hits = {}

    for read in genomes.reads:
        stats["total"] += 1
        read = read.upper()

        if len(read) < kmer_size:
            stats["too_short"] += 1
            continue

        # Collect candidate placements for this read
        candidates = {}

        total_read_kmers = len(read) - kmer_size + 1
        any_kmer_hit = False
        # Go through all k-mers in the read
        for read_pos in range(total_read_kmers):
            read_kmer = read[read_pos:read_pos + kmer_size]

            if read_kmer not in ar_kmers.kmers:
                continue
            any_kmer_hit = True

            # Gather information from database and build candidate starts (kmer voting)
            for gene_header, gene_pos, is_rc in ar_kmers.kmers[read_kmer]:
                candidate_start = gene_pos - read_pos
                key = (gene_header, candidate_start, is_rc)

                if key not in candidates:
                    candidates[key] = 0
                candidates[key] += 1

        if not any_kmer_hit:
            stats["no_kmer_hit"] += 1
        # No candidate placements found
        if not candidates:
            stats["no_candidates"] += 1
            continue

        # Choose the candidate with the most matched k-mers
        best_hit = None
        best_count = 0

        for key, count in candidates.items():
            gene_len_key = len(ar_kmers.gene_dict[key[0]])
            best_len = len(ar_kmers.gene_dict[best_hit[0]]) if best_hit else 0
            if count > best_count or (count == best_count and gene_len_key > best_len):
                best_hit = key
                best_count = count

        if best_count < 3:
            stats["best_count_lt3"] += 1
            continue

        gene_header, gene_start, is_rc = best_hit
        gene_seq = ar_kmers.gene_dict[gene_header]

        # Strand adjustment
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

        if overlap_len < kmer_size:
            stats["overlap_lt_k"] += 1
            continue

        read_segment = test_read[read_start:read_start + overlap_len]
        gene_segment = gene_seq[gene_start_clipped:gene_start_clipped + overlap_len]

        # Compute identity over the overlapping region
        matches = 0
        for a, b in zip(read_segment, gene_segment):
            if a == b:
                matches += 1

        identity = matches / overlap_len

        # Accept read if at least 50% identity
        if identity < 0.50:
            stats["identity_lt_70"] += 1
            continue

        stats["accepted"] += 1
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

###################################################################################################
    # Build the results that will get printed to the output
    results = []

    for gene_header, gene_seq in ar_kmers.gene_dict.items():
        gene_len = len(gene_seq)

        if gene_header not in hits:
            continue

        covered_bases = len(hits[gene_header])
        coverage_pct = (covered_bases / gene_len) * 100
        mean_depth = sum(hits[gene_header].values()) / covered_bases

    # Map expected IDs to their stats, regardless of thresholds
    expected_stats = {gid: [] for gid in expected_ids}

    for gene_header, gene_seq in ar_kmers.gene_dict.items():
        gene_len = len(gene_seq)

        if gene_header not in hits:
            # No reads mapped at all
            for gid in expected_ids:
                if gid in gene_header:
                    expected_stats[gid].append((gene_header, 0.0, 0.0))
            continue

        covered_bases = len(hits[gene_header])
        coverage_pct = (covered_bases / gene_len) * 100
        mean_depth = sum(hits[gene_header].values()) / covered_bases

        # Record stats for watch-list genes, even if they fail thresholds
        for gid in expected_ids:
            if gid in gene_header:
                expected_stats[gid].append((gene_header, coverage_pct, mean_depth))

    print("\n=== Expected genes debug ===")
    for gid in expected_ids:
        entries = expected_stats[gid]
        if not entries:
            print(gid, "→ NO MATCHING HEADER IN DB")
        else:
            for header, cov, depth in entries:
                print(f"{gid} matched header {header}: coverage={cov:.2f}%, depth={depth:.2f}x")
####################################################################################
    # Original filter
    if coverage_pct > 80 and mean_depth >= 5:
        results.append((gene_header, coverage_pct, mean_depth))

    results.sort(key=lambda x: (x[1], x[2]), reverse=True)

    with open(out_path, "w") as f:
        f.write("Gene,Coverage(%),Depth(X)\n")
        for gene_header, coverage_pct, mean_depth in results:
            f.write(f"{gene_header},{coverage_pct:.2f},{mean_depth:.2f}\n")

    print(f"Output written to {out_path}")
    print("Read stats:", stats)


if __name__ == "__main__":
    main()
