import warnings


class Fasta:

    # This is from last week, ignore everything in this block
    ##################################################################################################

    # Initialization, create empty lists and if filename is provided, load it
    def __init__(self, filename=None):
        self.headers, self.sequences = [], []
        if filename is not None:
            self.load(filename)

    # Adapted version of fastaread function to work with self.
    def load(self, filename):
        with open(filename, "r") as infile:
            temp_seq = ""
            for line in infile:
                line = line.strip()
                if line.startswith(">"):
                    self.headers.append(line)
                    if temp_seq:
                        self.sequences.append(temp_seq)
                        temp_seq = ""
                else:
                    temp_seq += line
            if temp_seq:
                self.sequences.append(temp_seq)

        # Check to make sure each header is unique
        if len(self.headers) > len(set(self.headers)):
            raise ValueError("There are repeated headers! Can't be")

    # This is a function to check that start and none are the correct kind of input.
    # i.e. they are an integer and above 0 for start or below the length of the lists for end
    def _check_input(self, start=None, end=None):

        if start is not None:
            if start <= 0:
                raise IndexError("Index must be positive")

            if not isinstance(start, int):
                raise ValueError("Start type is invalid, must be an integer")

            else:
                start = start - 1

        if end is not None:
            if end > len(self.sequences):
                raise IndexError("End index is too large, out of range")

            if not isinstance(end, int):
                raise ValueError("End type is invalid, must be an integer")

        return start, end

    # Return a copy of the headers and sequences on the selected range
    def content(self, start=None, end=None):

        start, end = self._check_input(start, end)

        # Here I define the 4 different cases depending on start and or end
        # and manipulate the lists accordingly
        if start is None and end is None:
            return self.headers.copy(), self.sequences.copy()

        elif start is not None and end is None:
            return [self.headers[start]].copy(), [self.sequences[start]].copy()

        elif start is None and end is not None:
            raise ValueError("Are you stupid why are you just using end, no")

        elif start is not None and end is not None:
            return self.headers[start:end].copy(), self.sequences[start:end].copy()

    # Modified version of fastawrite to work with self (peep the 60 character blocks)
    def save(self, filename):
        with open(filename, "w") as output:
            for i in range(len(self.headers)):
                output.write(f"{self.headers[i]}\n")
                for j in range(0, len(self.sequences[i]), 60):
                    output.write(self.sequences[i][j:j+60] + "\n")

    # Delete the selected sequences with slicing according to input of method
    def delete(self, start=None, end=None):

        start, end = self._check_input(start, end)

        if start is not None and end is not None:
            del self.sequences[start:end]
            del self.headers[start:end]

        elif start is not None and end is None:
            del self.sequences[start]
            del self.headers[start]

    # Introduce header(s) and sequence(s) on the position provided. If not provided, will be appended to end of lists
    def insert(self, header, sequence, position=None):
        # Check that both headers and sequences are the same type, we want to mantain a strong control on both lists being harmonic
        if not ((isinstance(header, str) and isinstance(sequence, str)) or (isinstance(header, list) and isinstance(sequence, list))):
            raise TypeError(
                "Arguments are different classes, that's not allowed in this method")

        # Second check to make sure they are the same length (Control!!!!)
        if isinstance(header, list) and len(header) != len(sequence):
            raise ValueError(
                "Length of headers and sequences is different, no!")

        # Handling of position argument
        if position is None:
            position = len(self.sequences) + 1

        if position is not None:
            if position <= 0:
                raise IndexError("Must insert in a positive integer")

            if position > len(self.sequences) + 1:
                raise IndexError(
                    "Trying to insert outside range, check length of lists")

            # Make it readable for python
            position = position - 1

            # Finally, handle inputs differently if they are strings or lists
            if isinstance(header, str):
                self.sequences.insert(position, sequence)
                self.headers.insert(position, header)

            elif isinstance(header, list):
                for i in reversed(range(len(header))):
                    self.headers.insert(position, header[i])
                    self.sequences.insert(position, sequence[i])

        if len(self.headers) > len(set(self.headers)):
            raise ValueError("There are repeated headers! Can't be")

    # Mini check function for the alphabet since it's used more than once

    def _is_valid(self, seq, alpha_set):
        return set(seq) <= alpha_set

    # Function to make sure the entries in a specific start and/or end position belong to the alphabet provided
    # Did not include pro version where alphabet is its own class

    def verify(self, alphabet, start=None, end=None):
        if not isinstance(alphabet, str):
            raise TypeError(
                "Alphabet is not a string of characters, that's not allowed")

        start, end = self._check_input(start, end)

        # As before, make sure each case is covered according to start/end input
        if start is not None and end is not None:
            unchecked_sequences = self.sequences[start:end]

        elif start is not None and end is None:
            unchecked_sequences = [self.sequences[start]]

        elif start is None and end is not None:
            raise ValueError("Are you stupid why are you just using end, no")

        elif start is None and end is None:
            unchecked_sequences = self.sequences

        # Create alphabet set (faster lookup speed) and check each sequence
        # Return true if it follow the alphabet and false if not.
        alpha_set = set(alphabet)
        for seq in unchecked_sequences:
            if not self._is_valid(seq, alpha_set):
                return False
        return True

    # This method checks specific entries (or all) and discards any that don't follow the provided alphabet,
    # which must be a string. Returns the list without the "faulty" entries.
    def discard(self, alphabet, start=None, end=None):
        if not isinstance(alphabet, str):
            raise TypeError(
                "Alphabet is not a string of characters, that's not allowed")

        start, end = self._check_input(start, end)

        # The idea in the following three blocks is that we break the entries into:
        #   before
        #   middle (sequences of interest to inspect)
        #   after
        # to be able to rebuild the sequences after
        if start is not None and end is not None:
            before_seq = self.sequences[:start]
            middle_seq = self.sequences[start:end]
            after_seq = self.sequences[end:]

            before_head = self.headers[:start]
            middle_head = self.headers[start:end]
            after_head = self.headers[end:]

        elif start is not None and end is None:
            before_seq = self.sequences[:start]
            middle_seq = [self.sequences[start]]
            after_seq = self.sequences[start + 1:]

            before_head = self.headers[:start]
            middle_head = [self.headers[start]]
            after_head = self.headers[start + 1:]

        elif start is None and end is None:
            before_seq = []
            middle_seq = self.sequences
            after_seq = []

            before_head = []
            middle_head = self.headers
            after_head = []

        # Check that sequences follow the alphabet and store in new variable
        new_seq = []
        new_head = []
        alpha_set = set(alphabet)

        for seq, header in zip(middle_seq, middle_head):
            if self._is_valid(seq, alpha_set):
                new_seq.append(seq)
                new_head.append(header)

        # Build the sequences and headers from the check before
        self.sequences = before_seq + new_seq + after_seq
        self.headers = before_head + new_head + after_head

######################################################################

############ Exercise 1 ############
    def __iter__(self):
        self.index = 0
        return self

    def __next__(self):
        if self.index == len(self.headers):
            raise StopIteration
        self.index += 1
        header, seq = self.headers[self.index-1], self.sequences[self.index-1]
        return header, seq

    def __len__(self):
        return len(self.sequences)

############ Exercise 2 ############

    def deletethis(self):
        del self.headers[self.index - 1]
        del self.sequences[self.index - 1]
        self.index -= 1

    def verifythis(self, alphabet):
        alpha_set = set(alphabet)
        if set(self.sequences[self.index - 1]) <= alpha_set:
            return True
        else:
            return False

    def insertthis(self, header, sequence):
        if not ((isinstance(header, str) and isinstance(sequence, str)) or (isinstance(header, list) and isinstance(sequence, list))):
            raise TypeError(
                "Arguments are different classes, that's not allowed in this method")

        if isinstance(header, list) and len(header) != len(sequence):
            raise ValueError(
                "Lists are different lengths, not possible in this class")

        if isinstance(header, str):
            self.headers.insert(self.index - 1, header)
            self.sequences.insert(self.index - 1, sequence)
            self.index += 1

        elif isinstance(header, list):
            for i, j in zip(reversed(header), reversed(sequence)):
                self.headers.insert(self.index - 1, i)
                self.sequences.insert(self.index - 1, j)
            self.index += len(header)

    def discardthis(self, alphabet):
        if not self.verifythis(alphabet):
            self.deletethis()

############ Exercise 3 ############

    def __iadd__(self, other):
        if not isinstance(other, Fasta):
            raise TypeError(
                "Trying to add a non-Fasta class object to a fasta object, not possible")

        self.headers.extend(other.headers)
        self.sequences.extend(other.sequences)

        return self

############ Exercise 4 ############

    def __add__(self, other):
        if not isinstance(other, Fasta):
            raise TypeError(
                "Trying to add a non-Fasta class object to a fasta object, not possible")

        new = Fasta()

        new.headers = self.headers.copy()
        new.sequences = self.sequences.copy()

        new.headers.extend(other.headers)
        new.sequences.extend(other.sequences)

        return new
