import gzip

class Genome:
    def __init__(self, file_list=None) -> None:
        self.file_list = file_list if file_list else []

    @property
    def reads(self):
        """
        Generator that lazily loads reads one by one.
        This solves 'How to handle the sequences (can’t be stored in a list)'
        """
        for filename in self.file_list:
            with gzip.open(filename, "rt") as f:
                for i, line in enumerate(f):
                    # In FASTQ, the sequence is on line 2 (index 1), line 6 (index 5), etc.
                    if i % 4 == 1:
                        yield line.strip()
