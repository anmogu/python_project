# Project 6: Antibiotic resistance gene detection
from kmer_class import Kmer
from parse_args import parse_args
import gzip


def main():
    # 1. Parse arguments
    db_path, out_path, fastq_files, kmer_size = parse_args()

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

    # Sort results by coverage first, then by average depth (highest first)
    genes95.sort(
        key=lambda x: (
            x[1],  # coverage_pct
            sum(ar_genes.coverage[x[0]]) /
            len(ar_genes.coverage[x[0]])  # avg_depth
        ),
        reverse=True
    )

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


def scan_genome(ar_genes, fastq_files, kmer_size, report_step=1000000):
    """Scan FASTQ reads and update per-base coverage on ar_genes."""

    def read_genome(fastq_files):
        comp_trans = str.maketrans("ATCGN", "TAGCN")
        valid_bases = set("ATGCN")

        for file in fastq_files:
            try:
                with gzip.open(file, "rt") as f:
                    line_num = 0

                    while True:
                        header = f.readline()
                        if not header:
                            break  # normal end of file

                        seq = f.readline()
                        plus = f.readline()
                        qual = f.readline()
                        line_num += 4

                        # Check incomplete FASTQ record
                        if not seq or not plus or not qual:
                            raise ValueError(
                                f"Incomplete FASTQ record in {file} near line {line_num}")

                        header = header.strip()
                        seq = seq.strip().upper()
                        plus = plus.strip()
                        qual = qual.strip()

                        # Check FASTQ format
                        if not header.startswith("@"):
                            raise ValueError(
                                f"Invalid FASTQ header in {file} near line {line_num - 3}"
                            )

                        if not plus.startswith("+"):
                            raise ValueError(
                                f"Missing '+' line in {file} near line {line_num - 1}"
                            )

                        # Check empty sequence
                        if len(seq) == 0:
                            raise ValueError(
                                f"Empty sequence in {file} near line {line_num - 2}"
                            )

                        # Check valid DNA alphabet
                        if any(base not in valid_bases for base in seq):
                            raise ValueError(
                                f"Invalid DNA sequence in FASTQ file {file} near line {line_num - 2}"
                            )

                        # Check sequence-quality length match
                        if len(seq) != len(qual):
                            raise ValueError(
                                f"Sequence/quality length mismatch in {file} near line {line_num}"
                            )

                        yield seq
                        yield seq[::-1].translate(comp_trans)

            except FileNotFoundError:
                raise FileNotFoundError(f"Cannot find FASTQ file: {file}")
            except OSError:
                raise ValueError(f"Could not read gzipped FASTQ file: {file}")

    def recursive_match(gene, position, gene_segment, read_segment):
        """
        Recursively match read_segment to gene_segment.
        If they are identical, increment coverage for all covered positions.
        If not, split in half and recurse on each half independently.
        """
        if gene_segment == read_segment:
            for p in range(position, position + len(read_segment)):
                ar_genes.coverage[gene][p] += 1
        elif len(read_segment) > 1:
            mid = len(read_segment) // 2
            recursive_match(
                gene, position, gene_segment[:mid], read_segment[:mid])
            recursive_match(gene, position + mid,
                            gene_segment[mid:], read_segment[mid:])

    for i, read in enumerate(read_genome(fastq_files)):
        if i % report_step == 0:
            print(f"  read {i:,}")

        # Try anchoring on the first k-mer of the read
        hits = ar_genes.kmer_lookup.get(read[:kmer_size])
        if hits is not None:
            for gene, pos in hits:
                gene_seq = ar_genes.seq_lookup[gene]
                end = min(pos + len(read), len(gene_seq))
                recursive_match(
                    gene, pos,
                    gene_seq[pos:end],
                    read[:end - pos]
                )
            continue  # skip last-kmer check if first kmer already matched

        # Try anchoring on the last k-mer of the read
        hits = ar_genes.kmer_lookup.get(read[-kmer_size:])
        if hits is not None:
            for gene, pos in hits:
                gene_seq = ar_genes.seq_lookup[gene]
                start = max(0, pos + kmer_size - len(read))
                L = pos + kmer_size - start
                recursive_match(
                    gene, start,
                    gene_seq[start:pos + kmer_size],
                    read[-L:]
                )


if __name__ == "__main__":
    main()
