#### Project 6 antibiotic resistance

from kmer_class import Kmer
import gzip


#set kmer size = 19
kmer = 19


print("Loading kmers")
#Load the kmers of the anti-biotics resistance genes through the filename included in the input
#ar_genes = Kmer(filename = "short_res_genes.fsa")
ar_genes = Kmer(filename = "resistance_genes.fsa")


print("Loading genome")
#Make a list of all the raw read files
genomes_to_load = ["Unknown3_raw_reads_1.txt.gz","Unknown3_raw_reads_2.txt.gz"]
#genomes_to_load = ["test_genome.txt.gz"]

outputfile = "out.csv"

#make generator that check entire genome for AR-gene kmers and yields reads with kmers in them
def scan_genome(genome_files, kmer_size = 19, report_step = 1000000):
    #function that lazily loads genes
    def read_genome(file_list):
        comp_trans = str.maketrans("ATCG", "TAGC")
        for file in file_list:
            with gzip.open(file,"rt") as f:
                for i,line in enumerate(f):
                    if i % 4 == 1:
                        yield line.strip()
                        yield line.strip()[::-1].translate(comp_trans)
                


    def recursive_match_read(gene, position,seq_to_match, read):
      
        #check if entire sequence is a match
        if seq_to_match == read:
            
            if position + len(read) > len(ar_genes.coverage[gene]):
                print("error")
                print(f"endpos{position + len(read)}")
                print(len(ar_genes.coverage[gene]))
                print(seq_to_match)
                print(read)
                print(position)

            #if yes, add 1 to all coverage positions
            for p in range(position, position + len(read)):
            #add 1 to all positions
                ar_genes.coverage[gene][p] += 1

            

            #if no, and length > 1
        elif len(read) > 1:
            #split read in the middle and rerun recursive_match_read
            new_pos = len(read) // 2

            #print(position)
            #print(position + new_pos)
            
            #run match on 2 new sequences
            recursive_match_read(gene,position, seq_to_match[:new_pos], read[:new_pos])
            recursive_match_read(gene,position + new_pos, seq_to_match[new_pos:], read[new_pos:])
            
                    
                    
    #for each read
    for i,g in enumerate(read_genome(genome_files)):
        if i% report_step == 0:
            print(f"read {i}")

        #try to look for first kmer in read
        ori_pos = ar_genes.kmer_lookup.get( g[:kmer_size] )
        if ori_pos is not None:
           
           #for each matched gene in ori_pos:
           for ori,pos in ori_pos:

               
               
               #match the parts of the genes that match, starting from the found kmer to kmer position + read length
               seq_to_check = ar_genes.seq_lookup[ori]

               end_position = min(pos + len(g), len(seq_to_check))

               
               potential_gene_match = seq_to_check[pos:end_position]

               #print(f"F, pos: {pos}, endpos: {end_position}, ")
               #print(potential_gene_match)
               #print(g)
               recursive_match_read(gene = ori, position= pos, seq_to_match = potential_gene_match, read = g[:end_position-pos] )
        else:
            #try last kmer
            ori_pos = ar_genes.kmer_lookup.get( g[-kmer_size:] )

            if ori_pos is not None:
                
                for ori,pos in ori_pos:
                    seq_to_check = ar_genes.seq_lookup[ori]
                    start_position = max(0, pos + kmer_size - len(g))
                    #match the parts of the genes that match, starting from the found kmer to kmer position - read length to the kmer at the end of the read
                    potential_gene_match = seq_to_check[start_position:(pos + kmer_size)]


                    #print(f"R, pos: {pos}, startpos: {start_position}, ")
                    #print(potential_gene_match)
                    #print(g)

                    L = pos + kmer_size - start_position
                    read = g[-L:]

                    recursive_match_read(gene = ori, position= start_position, seq_to_match = potential_gene_match, read = read )

        

        #if first or last kmer in gene kmers
        #if g[:kmer_size] in seq_kmers or g[-kmer_size:] in seq_kmers:
         #   #Scan entire read
         #   for i in range(0,len(g)-kmer_size):
         #       #yield kmer
         #       #if g[i:(i+kmer_size)] in seq_kmers:
         #       yield g[i:(i+kmer_size)]
        
                
                

print("Scanning genome")
#Run through genome and for each gene in AR-genes kmer add 1 to the count list with the matching kmer
scan_genome(genomes_to_load)



#We now have the counts for each kmer in the form of a list. 

#write a csv with all the outputs
print("Write output")
#ar_genes.save_kmer(outputfile)

avg_cutoff = 10
min_cutoff = 0
with open(outputfile, "w") as f:
    f.write(f"gene, average depth, minimum depth \n")
    for h in ar_genes.header:
        avg = sum(ar_genes.coverage[h])/len(ar_genes.coverage[h])
        minimum = min(ar_genes.coverage[h])
        if  avg > avg_cutoff and minimum > min_cutoff:
            f.write(f"{h}, {avg}, {minimum} \n")
            #f.write(f"{ar_genes.coverage[h]} \n")
  
            

#95% coverage
with open("percent_out", "w") as f:
    f.write(f"gene, average depth, minimum depth \n")
    for h in ar_genes.header:
        avg = sum(ar_genes.coverage[h])/len(ar_genes.coverage[h])
        minimum = min(ar_genes.coverage[h])
        gene_coverage = 0
        for c in ar_genes.coverage[h]:
            if c > avg_cutoff:
                gene_coverage += 1
        gene_coverage_percent = gene_coverage/len(ar_genes.coverage[h])*100
        if  gene_coverage_percent > 95 and avg > avg_cutoff and minimum > min_cutoff:
            f.write(f"{h}, {gene_coverage_percent}, {avg}, {minimum} \n")
            #f.write(f"{ar_genes.coverage[h]} \n")