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

    ## functions for class

    #function that reads a fastafile from filename
    def load(self,filename):
        #check if filename exists
        try:
                with open(filename,"r") as f:
                        lines = f.read().splitlines()

        # if filename does not exist print error message and exit
        except:
            raise ValueError("cannot find filename")

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

        # add dna_seq and dna_header to instance list header and sequence (here we assume that each instance only contains one fasta file and that we cannot concat them, but instead running the load function again will overwrite the fasta file)
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
            for i in range(0,len(s)-kmer_size):
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
     
    







    def save(self,filename):
        #check if the instance has the attribute sequence and header
        if hasattr(self, "sequence") and hasattr(self, "header"):
            #check if headers and sequences are the same length
            if len(self.header) != len(self.sequence):
                raise ValueError("header and sequence not of same length")
            else:
                #make a file
                with open(filename, "w") as f:
                    for h,s in zip(self.header,self.sequence):
                        f.write(h + "\n" + s + "\n")

    #exercise 5
    def content(self,start = None, end = None):
        #check if instance has a loaded fasta file
        if not hasattr(self,"header") and not hasattr(self,"sequence"):
            raise NameError("Instance does not have a loaded fasta file")

        #helper function to check that input is a whole number
        def WholeNumber(number):
            #check if it is a number
            if isinstance(number, (float,int)):
                #check if it is a whole number
                if not number % 1 == 0:
                    raise ValueError(f"not given a whole number, for number {number}")
            else:
                raise TypeError("not given a number")

        #if both given start and end (meaning they are no longer None)
        if start and end:
            #check that they are both whole numbers (function gives error if they are not both whole numbers)
            WholeNumber(start)
            WholeNumber(end)

            #return betweed start and end
            return (self.header[start:end],self.sequence[start:end])

        #if only given start,
        elif start:
            #check that start is whole number (function gives error if they are not both whole numbers)
            WholeNumber(start)

            #return from start
            return (self.header[start:],self.sequence[start:])
        else:
            #no input is given, return all
            return (self.header,self.sequence)

    #exercise 4 delete function
    def delete(self,start = None, end = None):
        #check if instance has a loaded fasta file
        if not hasattr(self,"header") and not hasattr(self,"sequence"):
            raise NameError("Instance does not have a loaded fasta file")

        #helper function to check that input is a whole number
        def WholeNumber(number):
            #check if it is a number
            if isinstance(number, (float,int)):
                #check if it is a whole number
                if not number % 1 == 0:
                    raise ValueError(f"not given a whole number, for number {number}")
            else:
                raise TypeError("not given a number")

        #if both given start and end (meaning they are no longer None)
        if start and end:
            #check that they are both whole numbers (function gives error if they are not both whole numbers)
            WholeNumber(start)
            WholeNumber(end)

            #save new lists
            self.header = self.header[:start] + self.header[end:]
            self.sequence = self.sequence[:start] + self.sequence[end:]

        #if only given start,
        elif start:
            #check that start is whole number (function gives error if they are not both whole numbers)
            WholeNumber(start)

            #save new lists
            self.header = self.header[:start] + self.header[start+1:]
            self.sequence = self.sequence[:start] + self.sequence[start+1:]
        else:
            #no input is given, meaning delete all
            self.header = None
            self.sequence = None

    #Exercise 6
    def insert(self,header, sequence,position = None):
        #check that given header and sequence are lists
        if not isinstance(header,list) and not isinstance(sequence,list):
            raise TypeError("all given objects are not lists")

        #check if instance has a header and sequnce attatched to it
        if not hasattr(self,"header") and not hasattr(self,"sequence"):
            #if it does not, we add these headers and sequences as the first
            self.header = header
            self.sequence = sequence
            #if a position was given, it is superflous as we do not have any elements in the list to place them between, so we just give a warning and ignore the position
            if position:
                print("Warning: instance does not have any fasta elemets, and the position was therefore ignored")
        else:
            #we will now add new sequences
            if position:
                #if given position add new sequences at position
                self.header = self.header[:position] + header + self.header[position:]
                self.sequence = self.sequence[:position] + sequence + self.sequence[position]

            else:
                #add sequences at the end
                self.header = self.header + header
                self.sequence = self.sequence + sequence

    def verify(self,alphabet = "DNA", start = None, end = None):
        #check if instance has a loaded fasta file
        if not hasattr(self,"header") and not hasattr(self,"sequence"):
            raise NameError("Instance does not have a loaded fasta file")

        #check if given alphabet is in class variables
        if not alphabet in Fasta.alphabets.keys():
            raise ValueError("alphabet not existing")

        #helper function to check that input is a whole number
        def WholeNumber(number):
            #check if it is a number
            if isinstance(number, (float,int)):
                #check if it is a whole number
                if not number % 1 == 0:
                    raise ValueError(f"not given a whole number, for number {number}")
            else:
                raise TypeError("not given a number")
        #Helper function to check string for only having letters from alphabet
        def checkAlphabet(string, alphabet):
            #returns true if all characters ,c in string is in the alphabet
            return all(c in Fasta.alphabets[alphabet] for c in string.upper())


        #if both given start and end (meaning they are no longer None)
        if start and end:
            #check that they are both whole numbers (function gives error if they are not both whole numbers)
            WholeNumber(start)
            WholeNumber(end)

            #check if sequences have legal alphabets
            sequences_to_check = self.sequence[start:end]

        #if only given start,
        elif start:
            #check that start is whole number (function gives error if they are not both whole numbers)
            WholeNumber(start)

            #check if sequences have legal alphabets
            sequences_to_check = self.sequence[start:]
        else:
            #no input is given, check all
           #check if sequences have legal alphabets
            sequences_to_check = self.sequence

        #now check sequences to check for alphabet
        return all(checkAlphabet(seq, alphabet=alphabet) for seq in sequences_to_check)

    #Exercise 8
    def discard(self,alphabet = "DNA", start = None, end = None):
        #check if instance has a loaded fasta file
        if not hasattr(self,"header") and not hasattr(self,"sequence"):
            raise NameError("Instance does not have a loaded fasta file")

        #check if given alphabet is in class variables
        if not alphabet in Fasta.alphabets.keys():
            raise ValueError("alphabet not existing")

        #helper function to check that input is a whole number
        def WholeNumber(number):
            #check if it is a number
            if isinstance(number, (float,int)):
                #check if it is a whole number
                if not number % 1 == 0:
                    raise ValueError(f"not given a whole number, for number {number}")
            else:
                raise TypeError("not given a number")

        ##identify the sequences to check
        #if both given start and end (meaning they are no longer None)
        if start and end:
            #check that they are both whole numbers (function gives error if they are not both whole numbers)
            WholeNumber(start)
            WholeNumber(end)

            #check if sequences have legal alphabets
            sequences_to_check = self.sequence[start:end]
            headers_to_check = self.header[start:end]

        #if only given start,
        elif start:
            #check that start is whole number (function gives error if they are not both whole numbers)
            WholeNumber(start)

            #check if sequences have legal alphabets
            sequences_to_check = self.sequence[start:]
            headers_to_check = self.header[start:]
        else:
            #no input is given, check all
           #check if sequences have legal alphabets
            sequences_to_check = self.sequence
            headers_to_check = self.header

        #Helper function to check string for only having letters from alphabet
        def checkAlphabet(string, alphabet):
            #returns true if all characters ,c in string is in the alphabet
            return all(c in Fasta.alphabets[alphabet] for c in string.upper())

        #now go through sequences to check, and only save sequences that are verified
        verified_pairs = [
            (h,s) for h,s in zip(headers_to_check,sequences_to_check)
            if checkAlphabet(s,alphabet)
            ]
        #split the tuple into two lists
        verified_sequences = [s for h,s in verified_pairs]
        verified_headers = [h for h,s in verified_pairs]

        #graft original sequence list with the verified sequences
        if start and end:
            self.sequence = self.sequence[:start] + verified_sequences + self.sequence[end:]
            self.header = self.header[:start] + verified_headers + self.header[end:]

        #if only given start,
        elif start:
            self.sequence = self.sequence[:start] + verified_sequences
            self.header = self.header[:start] + verified_headers

        else:
            #no input is given, check all
            self.sequence = verified_sequences
            self.header = verified_headers

########################  Exercises starts here   ####################################

    ##Exercise 1
    #make a function for length
    def __len__(self):
        return len(self.header)

    # make instance iterable
    #def __iter__(self):
    #    return iter(zip(self.header, self.sequence))

    #Exercise 2
    

    #make functions
    def deletethis(self):
        #check that object currently being iterated
        if not self.iterationPointer:
            raise TypeError("instance is not being iterated")
        else:
            #delete the current object(pointer is pointing to the next object)
            del self.header[self.iterationPointer-1]
            del self.sequence[self.iterationPointer-1]
            # The list is one object shorter so the pointer should point to the element one short
            self.iterationPointer -= 1

    def insertthis(self,header,sequence):
        #check that object currently being iterated
        if not self.iterationPointer:
            raise TypeError("instance is not being iterated")

        #check that header and sequence are strings
        elif not isinstance(header, str) and isinstance(sequence,str):
            raise TypeError("input are not strings")

        else:
            #add the string and header as the next elements
            self.header = self.header[:self.iterationPointer] + [header] + self.header[self.iterationPointer:]
            self.sequence = self.sequence[:self.sequence] + [sequence] + self.sequence[self.iterationPointer:]
            # add 1 to pointer so that the new element is not iterated over
            self.iterationPointer += 1

    def verifythis(self,alphabet):
        #use existing verify function on current sequence
        return type(self).verify(self,alphabet,self.iterationPointer-1, self.iterationPointer)

    def discardthis(self, alphabet):
        type(self).discard(alphabet= alphabet, start = self.iterationPointer-1, end = self.iterationPointer)

    ### Exericse 3
    def __iadd__(self, other):

        #check that other is the same type as instance
        if type(self) != type(other):
            raise TypeError("the added objects are not of same type")

        #check if other has the header and sequence attribute
        elif not hasattr(other,"header") and not hasattr(other,"sequence"):
            raise NameError("The added object does not have header and sequence")

        #check if instance has header and sequence
        elif not hasattr(self,"header") and not hasattr(self,"sequence"):
            #make header and sequence of other the current header and sequence
            self.header = other.header
            self.sequence = other.sequence

        else:
            self.header = self.header + other.header
            self.sequence = self.sequence + other.sequence
        return self

    #Exercise 4
    def __add__(self,other):
        #check that other is the same type as instance
        if type(self) != type(other):
            raise TypeError("the added objects are not of same type")

        #check if other has the header and sequence attribute
        elif not hasattr(other,"header") and not hasattr(other,"sequence"):
            #make header and sequence of self the current header and sequence
            new_instance = type(self)()
            new_instance.insert(self.header, self.sequence )

        #check if instance has header and sequence
        elif not hasattr(self,"header") and not hasattr(self,"sequence"):
            #make header and sequence of other the current header and sequence
            new_instance = type(self)()
            new_instance.insert(other.header, other.sequence)

        else:
            new_instance = type(self)()
            new_instance.insert(self.header + other.header, self.sequence + other.sequence)
        return new_instance
