class Kmer:
    """Resistance gene k-mer database for read mapping."""

    def __init__(self, filename=None, kmer_size=19):
        self.kmer_size = kmer_size
        if filename is not None:
            self.load(filename)
            self.split(kmer_size=kmer_size)

    def load(self, filename):
        """Parse a FASTA file and store headers, sequences, coverage arrays, and seq lookup."""
        try:
            with open(filename, "r") as f:
                lines = f.read().splitlines()
        except FileNotFoundError:
            raise FileNotFoundError(f"Cannot find file: {filename}")

        headers = []
        sequences = []
        seq = ""
        first = True

        for line in lines:
            if line.startswith(">"):
                if not first:
                    sequences.append(seq.upper())
                else:
                    first = False
                headers.append(line)
                seq = ""
            else:
                seq += line
        sequences.append(seq.upper())  # Append last sequence

        self.header = headers
        self.sequence = sequences
        self.seq_lookup = {h: s for h, s in zip(headers, sequences)}
        self.coverage = {h: [0] * len(s) for h, s in zip(headers, sequences)}

    def split(self, kmer_size=19):
        """Build kmer_lookup and gene_to_kmer dictionaries from loaded sequences."""
        kmer_lookup = {}
        gene_to_kmer = {}

        for h, s in zip(self.header, self.sequence):
            gene_to_kmer[h] = []
            for i in range(len(s) - kmer_size + 1):  
                kmer = s[i:i + kmer_size]   # Slide a window of kmer_size to extract all

                if kmer not in kmer_lookup:
                    kmer_lookup[kmer] = []
                kmer_lookup[kmer].append((h, i))

                gene_to_kmer[h].append((kmer, i))

        self.kmer_lookup = kmer_lookup
        self.gene_to_kmer = gene_to_kmer

    def __len__(self):
        return len(self.header)

    def __repr__(self):
        return f"Kmer(genes={len(self.header)}, kmer_size={self.kmer_size})"