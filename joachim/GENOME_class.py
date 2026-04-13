import gzip

class Genome:
    
    ##initate class
    def __init__(self, file_list = None) -> None:
        #if given a filename when instanciated, load that file into the instance
        if file_list:
            reads = []
            for f in file_list:
                reads += type(self).read_genome(f)
        self.reads = reads


    #function that reads a fastafile from filename
    def read_genome(filename, genes = None):

        with gzip.open(filename, "rt") as f:
            lines = f.read().splitlines()

        genome = []
        #for every 4 lines save line 2 to genome
        for i in range(0,int(len(lines)/4)):
            genome.append(lines[i*4+1])
        return genome

    def subset(self, kmers):
        self.reads = [r for r in self.reads if any(kmer in r for kmer in set(kmers)) ]
    
    # split genome into kmers
    def split(self, kmer_size = 19):
        print("Splitting genome")
        seqs = []
        

        #go through each sequence 
        for s in self.reads:
            # from first position to last - kmer_size, make kmers
            for i in range(0,len(s)-kmer_size):
                #add kmer to list
                seqs.append(s[i:(i+kmer_size)])
                #add origin and position


        self.kmers = seqs