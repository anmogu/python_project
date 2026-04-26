class Kmer:
    # Initate class
    def __init__(self, filename=None, kmer_size=19) -> None:
        # If given a filename when instanciated, load that file into the instance and run split on it (see below)
        if filename:
            self.load(filename)
            self.split(kmer_size=kmer_size)
        self.kmer_size = kmer_size

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
