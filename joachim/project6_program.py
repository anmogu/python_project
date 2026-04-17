#### Project 6 antibiotic resistance

from kmer_class import Kmer
from GENOME_class import Genome
import gzip


#set kmer size = 19
kmer = 19


print("Loading kmers")
#Load the kmers of the anti-biotics resistance genes through the filename included in the input
ar_genes = Kmer(filename = "short_res_genes.fsa")


print("Loading genome")
#Make a list of all the raw read files
#genomes_to_load = ["Unknown3_raw_reads_1.txt.gz","Unknown3_raw_reads_2.txt.gz"]
genomes_to_load = ["test_genome.txt.gz"]

outputfile = "out.csv"

#make generator that check entire genome for AR-gene kmers and yields reads with kmers in them
def scan_genome(genome_files, seq_kmers, kmer_size = 19, report_step = 1000):
    #function that lazily loads genes
    def read_genome(file_list):
        for file in file_list:
            with gzip.open(file,"rt") as f:
                for i,line in enumerate(f):
                    if i % 4 == 1:
                        yield line.strip()
                    
                    

    for i,g in enumerate(read_genome(genome_files)):
    #    if i% report_step == 0:
    #        print(f"read {i} of {len_genome}")

        #if first or last kmer in gene kmers
        if g[:kmer_size] in seq_kmers or g[-kmer_size:] in seq_kmers:
            #Scan entire read
            for i in range(0,len(g)-kmer_size):
                #yield kmer
                #if g[i:(i+kmer_size)] in seq_kmers:
                yield g[i:(i+kmer_size)]
        
                
                

print("Scanning genome")
#Run through genome and for each gene in AR-genes kmer add 1 to the count list with the matching kmer
for s in scan_genome(genomes_to_load, ar_genes.kmers):
    ar_genes.count_kmer(s)



#We now have the counts for each kmer in the form of a list. 

#write a csv with all the outputs
print("Write output")
ar_genes.save_kmer(outputfile)

