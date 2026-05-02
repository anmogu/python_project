from kmer_clean import Kmer
from parse_args import parse_args
import gzip

def read_genome(fastq_files):
    comp_trans = str.maketrans("ATCG", "TAGC")
    valid_bases = set("ATGC")
    for file in fastq_files:
        with gzip.open(file, "rt") as f:
            for i, line in enumerate(f):
                if i % 4 == 1:
                    seq = line.strip()
                    seq = seq.upper()
                    if any(base not in valid_bases for base in seq):
                        raise ValueError("Invalid DNA sequence in FASTQ file")
                    yield seq
                    yield seq[::-1].translate(comp_trans)

def recursive_match(gene, position, gene_segment, read_segment, ar_genes):
    """
    Standalone recursive function. ar_genes is passed explicitly.
    """
    if gene_segment == read_segment:
        for p in range(position, position + len(read_segment)):
            ar_genes.coverage[gene][p] += 1
    elif len(read_segment) > 1:
        mid = len(read_segment) // 2
        recursive_match(gene, position, gene_segment[:mid], read_segment[:mid], ar_genes)
        recursive_match(gene, position + mid, gene_segment[mid:], read_segment[mid:], ar_genes)

def scan_genome(ar_genes, fastq_files, kmer_size, report_step=1000000):
    for i, read in enumerate(read_genome(fastq_files)):
        if i % report_step == 0:
            print(f"  read {i:,}")

        # Try anchoring on the first k-mer
        hits = ar_genes.kmer_lookup.get(read[:kmer_size])
        if hits is not None:
            for gene, pos in hits:
                gene_seq = ar_genes.seq_lookup[gene]
                end = min(pos + len(read), len(gene_seq))
                recursive_match(gene, pos, gene_seq[pos:end], read[:end - pos], ar_genes)
            continue

        # Try anchoring on the last k-mer
        hits = ar_genes.kmer_lookup.get(read[-kmer_size:])
        if hits is not None:
            for gene, pos in hits:
                gene_seq = ar_genes.seq_lookup[gene]
                start = max(0, pos + kmer_size - len(read))
                L = pos + kmer_size - start
                recursive_match(gene, start, gene_seq[start:pos + kmer_size], read[-L:], ar_genes)

def main():
    # 1. Parse arguments
    db_path, out_path, fastq_files, kmer_size, indent_size = parse_args()

    # 2. Load resistance gene k-mer database
    print(f"Loading kmers from {db_path} (k={kmer_size})...")
    ar_genes = Kmer(filename=db_path, kmer_size=kmer_size)

    # 3. Scan genome reads and build coverage
    print("Scanning genome...")
    scan_genome(ar_genes, fastq_files, kmer_size)

    # 4. Call genes with >95% positions covered at depth > avg_cutoff
    print("Calling genes...")
    avg_cutoff = 10
    min_cutoff = 3

    genes95 = []
    for h in ar_genes.header:
        cover = ar_genes.coverage[h]
        gene_len = len(cover)
        well_covered = 0
        for c in cover:
            if c > avg_cutoff:
                well_covered += 1
        coverage_pct = well_covered / gene_len * 100
        if coverage_pct > 95:
            genes95.append((h, coverage_pct))

    # 5. Remove subset genes using covered k-mer sets
    print("Removing subset genes...")
    gene_sets = []
    for h, coverage_pct in genes95:
        cover = ar_genes.coverage[h]
        kmer_set = set()
        for k, p in ar_genes.gene_to_kmer[h]:
            if min(cover[p:p + kmer_size]) > min_cutoff:
                kmer_set.add(k)
        gene_sets.append(kmer_set)

    x = 0
    while x < len(gene_sets):
        removed = False
        for y in range(len(gene_sets)):
            if x != y and gene_sets[x].issubset(gene_sets[y]):
                del gene_sets[x]
                del genes95[x]
                x -= 1
                removed = True
                break
        x += 1

    # 6. Write output
    print(f"Writing output to {out_path}...")
    with open(out_path + ".csv", "w") as f:
        f.write("Gene,Coverage (%),Average depth\n")
        for h, coverage_pct in genes95:
            cover = ar_genes.coverage[h]
            avg_depth = sum(cover) / len(cover)
            if avg_depth >= avg_cutoff and coverage_pct >= 95:
                f.write(f"{h},{coverage_pct:.2f},{avg_depth:.2f}\n")

    print("Done.")


if __name__ == "__main__":
    main()