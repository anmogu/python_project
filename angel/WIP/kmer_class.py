class Kmer:
    # Initate class
    def __init__(self, filename=None, kmer_size=19) -> None:
        # If given a filename when instanciated, load that file into the instance and run split on it (see below)
        if filename:
            self.load(filename)
            self.split(kmer_size=kmer_size)
        self.kmer_size = kmer_size

    # function that reads a fastafile from filename

    def load(self, filename):
        # Check if filename exists
        try:
            with open(filename, "r") as f:
                lines = f.read().splitlines()

        # If filename does not exist print error message and exit
        except:
            raise FileNotFoundError("Cannot find filename")

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
                # reset the sequence string to blank
                seq = ""
            else:
                # Append the line to the sequence
                seq += l
        # Append last sequence
        dna_seq.append(seq)

        self.sequence = dna_seq
        self.header = dna_header

    # Small function to return reverse complement strand of DNA
    def rev_comp(self, seq):
        comp_table = str.maketrans("ATCGN", "TAGCN")
        return seq.translate(comp_table)[::-1]

    def split(self, kmer_size=19):

        self.kmers = {}
        self.gene_dict = {}

        for h, s in zip(self.header, self.sequence):
            s = s.upper()
            self.gene_dict[h] = s

            for i in range(len(s) - kmer_size + 1):

                tmp_kmer = s[i:i + kmer_size]
                rev_tmp_kmer = self.rev_comp(tmp_kmer)

                if tmp_kmer not in self.kmers:
                    self.kmers[tmp_kmer] = []
                self.kmers[tmp_kmer].append((h, i, False))

                if rev_tmp_kmer not in self.kmers:
                    self.kmers[rev_tmp_kmer] = []
                self.kmers[rev_tmp_kmer].append((h, i, True))

