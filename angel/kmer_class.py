class Kmer:
    #Class variables
    alphabets = {"DNA": "ATCG",
                 "RNA": "AUCG",
                 "Protein": "ARNDCEQGHILKMFPOUSTWYV"}

    ##initate class
    def __init__(self, filename = None, kmer_size = 19) -> None:
        #if given a filename when instanciated, load that file into the instance
        if filename:
            self.load(filename)
            #split current header/seqs into kmers with the variables: origin, position, count, seq
            self.split(kmer_size=kmer_size)
            self.kmer_size = kmer_size

    #function that reads a fastafile from filename
    def load(self,filename):
        #check if filename exists
        try:
                with open(filename,"r") as f:
                        lines = f.read().splitlines()

        # if filename does not exist print error message and exit
        except:
            raise FileNotFoundError("Cannot find filename")

        #initialize name and sequence list
        dna_header = []
        dna_seq = []

        #initialize loop logic
        seq = ""
        first_sequence = True

        #go through all lines
        for l in lines:
                #check if a new header is detected
                if l.startswith(">"):
                        #append seq to list if it is not the first sequence
                        if not first_sequence:
                                dna_seq.append(seq)
                        else:
                                first_sequence = False

                        #append the new header
                        dna_header.append(l)
                        #reset the sequence string to blank
                        seq = ""
                else:
                        #apppend the line to the sequence
                        seq += l
        #append last sequence
        dna_seq.append(seq)
        
        self.sequence = dna_seq
        self.header = dna_header

    def split(self, kmer_size = 19):

        origins = []
        positions = []
        seqs = []
        count = []

        #go through each sequence 
        for h,s in zip(self.header,self.sequence):
            # from first position to last - kmer_size, make kmers
            for i in range(0,len(s)-kmer_size + 1):
                #add kmer to list
                seqs.append(s[i:(i+kmer_size)])
                #add origin and position
                origins.append(h)
                positions.append(i)
                count.append(0)


        self.kmers = seqs
        self.positions = positions
        self.origin = origins
        self.count = count

    

    #make functions
    def kmer_match(self):
        #check that object currently being iterated
        if not self.iterationPointer:
            raise TypeError("instance is not being iterated")
        else:
            #add 1 to current kmer count
            self.count[self.iterationPointer-1] += 1

    def scan_genome(self,genome_list, report_step = 1000):
        len_genome = len(genome_list)
        for i,g in enumerate(genome_list):
            if i% report_step == 0:
                print(f"read {i} of {len_genome}")
            for i in range(0,len(g)-self.kmer_size):
                    #yield kmer
                    if g[i:(i+self.kmer_size)] in self.kmers:
                        yield g[i:(i+self.kmer_size)]

    
    def count_kmer(self,kmer):
        for i,s in enumerate(self.kmers):
            if s == kmer:
                self.count[i] += 1
     
    def save_kmer(self,filename):
        #check if the instance has the attribute sequence and header
        if hasattr(self, "kmers") and hasattr(self, "origin") and hasattr(self, "origin") and hasattr(self, "origin"):
            ########UNITTEST HERE========
            #check if headers and sequences are the same length
            if len(self.kmers) != len(self.count):
                ########UNITTEST HERE
                raise ValueError("Lengths do not match")
            else:
                #make a file
                print(f"writing file{filename}")
                with open(filename, "w") as f:
                    #write header
                    f.write("Kmer" + "," + "Gene of origin" + "," + "Position in gene" + ","+ "Count"+ "\n")
                    for k,o,p,c in zip(self.kmers,self.origin, self.positions, self.count):
                        f.write(k + "," + str(o) + "," + str(p) + ","+ str(c)+ "\n")
