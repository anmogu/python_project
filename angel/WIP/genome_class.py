import gzip

# Initiate class


class Genome:
    def __init__(self, file_list=list) -> None:
        self.file_list = file_list

    @property
    def reads(self):
        """Generator that yields reads one by one"""
        for filename in self.file_list:
            with gzip.open(filename, "rt") as infile:
                for i, line in enumerate(infile):
                    # In FASTQ, the sequence is on line 2 (index 1), line 6 (index 5), etc.
                    if i % 4 == 1:
                        yield line.strip()
