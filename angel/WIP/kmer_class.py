class Kmer:
    # Initate class
    def __init__(self, filename=None, kmer_size=19) -> None:
        # If given a filename when instanciated, load that file into the instance and run split on it (see below)
        self.kmer_size = kmer_size
        if filename:
            self.load(filename)
            self.split(kmer_size=kmer_size)
            self.remove_subset_genes()
            self.remove_similar_genes(similarity_threshold=0.95)

    def load(self, filename):
        """Reads a fastafile from filename"""
        # Check if filename exists
        try:
            with open(filename, "r") as f:
                lines = f.read().splitlines()

        # If filename does not exist print error message and exit
        except:
            raise FileNotFoundError("Cannot find file or it doesn't exist")

        # initialize name and sequence list
        dna_header = []
        dna_seq = []

        # Initialize loop logic
        seq = ""
        first_sequence = True

        # Go through all lines
        for l in lines:
            # Check if a new header is detected
            if l.startswith(">"):
                # Append seq to list if it is not the first sequence
                if not first_sequence:
                    dna_seq.append(seq)
                else:
                    first_sequence = False

                # Append the new header
                dna_header.append(l)
                # Reset the sequence string to blank
                seq = ""
            else:
                # Append the line to the sequence
                seq += l
        # Append last sequence
        dna_seq.append(seq)

        self.sequence = dna_seq
        self.header = dna_header

    def rev_comp(self, seq):
        """Returns reverse complement strand of DNA only"""
        comp_table = str.maketrans("ATCGN", "TAGCN")
        return seq.translate(comp_table)[::-1]

    def split(self, kmer_size=19):
        """Creates 2 dictionaries with kmers from given sequences linked to header (gene name in this case)"""

        self.kmers = {}
        self.gene_dict = {}

        for h, s in zip(self.header, self.sequence):
            s = s.upper()
            self.gene_dict[h] = s

            # Create all kmers from the sequences
            for i in range(len(s) - kmer_size + 1):

                # Including reverse sequence
                tmp_kmer = s[i:i + kmer_size]
                rev_tmp_kmer = self.rev_comp(tmp_kmer)

                # Create an entry in the dictionary if it didn't exist
                if tmp_kmer not in self.kmers:
                    self.kmers[tmp_kmer] = []
                self.kmers[tmp_kmer].append((h, i, False))

                if rev_tmp_kmer not in self.kmers:
                    self.kmers[rev_tmp_kmer] = []
                self.kmers[rev_tmp_kmer].append((h, i, True))

    def remove_subset_genes(self):
        gene_kmers = {}

        for h, s in self.gene_dict.items():
            kmers = set()
            for i in range(len(s) - self.kmer_size + 1):
                tmp_kmer = s[i:i + self.kmer_size]
                rev_tmp_kmer = self.rev_comp(tmp_kmer)
                kmers.add(tmp_kmer)
                kmers.add(rev_tmp_kmer)
            gene_kmers[h] = kmers

        to_remove = set()

        headers = list(gene_kmers.keys())
        for i, h1 in enumerate(headers):
            if h1 in to_remove:
                continue
            for j, h2 in enumerate(headers):
                if i == j or h2 in to_remove:
                    continue

                kmers1 = gene_kmers[h1]
                kmers2 = gene_kmers[h2]

                if kmers1 < kmers2:
                    to_remove.add(h1)
                elif kmers2 < kmers1:
                    to_remove.add(h2)

        for h in to_remove:
            self.gene_dict.pop(h, None)
            self.sequence = [s for hdr, s in zip(
                self.header, self.sequence) if hdr != h]
            self.header = [hdr for hdr in self.header if hdr != h]

        self.kmers = {}

        for h, s in zip(self.header, self.sequence):
            for i in range(len(s) - self.kmer_size + 1):

                tmp_kmer = s[i:i + self.kmer_size]
                rev_tmp_kmer = self.rev_comp(tmp_kmer)

                if tmp_kmer not in self.kmers:
                    self.kmers[tmp_kmer] = []
                self.kmers[tmp_kmer].append((h, i, False))
                if rev_tmp_kmer not in self.kmers:
                    self.kmers[rev_tmp_kmer] = []
                self.kmers[rev_tmp_kmer].append((h, i, True))

    def remove_similar_genes(self, similarity_threshold=0.95):
        """Remove genes that share more than `similarity_threshold` of their kmers with another gene, keeping the longest."""
        gene_kmers = {}
        for h, s in self.gene_dict.items():
            kmers = set()
            for i in range(len(s) - self.kmer_size + 1):
                tmp_kmer = s[i:i + self.kmer_size]
                kmers.add(tmp_kmer)
                kmers.add(self.rev_comp(tmp_kmer))
            gene_kmers[h] = kmers

        to_remove = set()
        headers = list(gene_kmers.keys())

        for i, h1 in enumerate(headers):
            if h1 in to_remove:
                continue
            for j, h2 in enumerate(headers):
                if i >= j or h2 in to_remove:
                    continue

                kmers1 = gene_kmers[h1]
                kmers2 = gene_kmers[h2]

                intersection = len(kmers1 & kmers2)
                smaller = min(len(kmers1), len(kmers2))

                if intersection / smaller >= similarity_threshold:
                    # Keep the longer gene, discard the shorter
                    if len(self.gene_dict[h1]) >= len(self.gene_dict[h2]):
                        to_remove.add(h2)
                    else:
                        to_remove.add(h1)
                        break  # h1 is gone, stop inner loop

        # Reuse the existing removal + kmer rebuild logic
        for h in to_remove:
            self.gene_dict.pop(h, None)
            self.sequence = [s for hdr, s in zip(
                self.header, self.sequence) if hdr != h]
            self.header = [hdr for hdr in self.header if hdr != h]

        # Rebuild kmer index
        self.kmers = {}
        for h, s in zip(self.header, self.sequence):
            for i in range(len(s) - self.kmer_size + 1):
                tmp_kmer = s[i:i + self.kmer_size]
                rev_tmp_kmer = self.rev_comp(tmp_kmer)
                if tmp_kmer not in self.kmers:
                    self.kmers[tmp_kmer] = []
                self.kmers[tmp_kmer].append((h, i, False))
                if rev_tmp_kmer not in self.kmers:
                    self.kmers[rev_tmp_kmer] = []
                self.kmers[rev_tmp_kmer].append((h, i, True))
