import pytest
import sys
import gzip

sys.path.append("../../final_code")

from kmer_clean import Kmer
from alternative_algorithm_for_testing import recursive_match
from alternative_algorithm_for_testing import read_genome
from alternative_algorithm_for_testing import scan_genome 
from alternative_algorithm_for_testing import main 


class MockArGenes:
    def __init__(self):
        # Define a small coverage array to test against
        self.coverage = {">GeneA": [0, 0, 0, 0, 0, 0]}

@pytest.fixture
def test_setup():
    # This setup returns the objects needed for any algorithm test
    return MockArGenes()

@pytest.fixture
def gzipped_fastq(tmp_path):
    # Since the function uses gzip, we need to also create it in this format
    raw_content = "@read1\nAAAA\n+\n!!!!\n@read2\nGGGG\n+\n!!!!\n"
    
    file_path = tmp_path / "reads.fastq.gz"
    
    with gzip.open(file_path, "wt") as f:
        f.write(raw_content)
        
    return str(file_path)

# Copied from test_arg_pars.py
@pytest.fixture
def kmer_db(tmp_path):
    db_file = tmp_path / "genes.fasta"
    db_file.write_text(">GeneA\nAAAAGGGG\n>GeneB\nACGTACGT\n")
    
    # Return an initialized Kmer object
    return Kmer(filename=str(db_file), kmer_size=4)

@pytest.fixture
def kmer_db_file(tmp_path):
    # This is for the main pipeline test
    db_file = tmp_path / "genes.fasta"
    db_file.write_text(">GeneA\nAAAAGGGG\n>GeneB\nACGTACGT\n")
    return str(db_file)

#############################################################
##################### Sunny day testing #####################
#############################################################

def test_recursive_match_correct(test_setup):
    ar_genes = test_setup
    recursive_match(">GeneA", 0, "AAAA", "AAAA", ar_genes)
    assert ar_genes.coverage[">GeneA"] == [1, 1, 1, 1, 0, 0]

def test_read_genome_correct(gzipped_fastq):
    file_list = [gzipped_fastq]
    read_count = 0
    all_reads = []
    for read in read_genome(file_list):
        read_count += 1
        all_reads.append(read)
    assert read_count == 4
    assert all_reads[0] == "AAAA"
    assert all_reads[1] == "TTTT"
    assert all_reads[2] == "GGGG"
    assert all_reads[3] == "CCCC"

def test_scan_genome_correct(gzipped_fastq, kmer_db):
    ar_db = kmer_db
    file_list = [gzipped_fastq]
    scan_genome(ar_db, file_list, 4)
    assert ar_db.coverage[">GeneA"][0] == 1
    assert ar_db.coverage[">GeneA"][-1] == 1
    assert ar_db.coverage[">GeneB"][0] == 0

def test_main_pipeline_correct(tmp_path, gzipped_fastq, kmer_db_file):
    out_base = tmp_path / "result"

    sys.argv = [
        "program",
        "-db", str(kmer_db_file),
        "-o", str(out_base),
        "-i", str(gzipped_fastq),
        "-k", "4"
    ]

    main()
    output_file = tmp_path / "result.csv"
    assert output_file.exists()

    content = output_file.read_text()
    assert "Gene,Coverage (%),Average depth" in content

#############################################################
##################### Rainy day testing #####################
#############################################################

@pytest.fixture
def gzipped_fastq_lowercase(tmp_path):
    # Since the function uses gzip, we need to also create it in this format
    raw_content = "@read1\naaaa\n+\n!!!!\n@read2\nGGGG\n+\n!!!!\n"
    
    file_path = tmp_path / "reads.fastq.gz"
    
    with gzip.open(file_path, "wt") as f:
        f.write(raw_content)
        
    return str(file_path)

def test_recursive_match_with_mismatch(test_setup):
    ar_genes = test_setup
    recursive_match(">GeneA", 0, "AAAAAA", "AAAGGA", ar_genes)
    assert ar_genes.coverage[">GeneA"] == [1, 1, 1, 0, 0, 1]

def test_read_genome_lowercase(gzipped_fastq_lowercase):
    file_list = [gzipped_fastq_lowercase]
    read_count = 0
    all_reads = []
    for read in read_genome(file_list):
        read_count += 1
        all_reads.append(read)
    assert read_count == 4
    assert all_reads[0] == "AAAA"
    assert all_reads[1] == "TTTT"
    assert all_reads[2] == "GGGG"
    assert all_reads[3] == "CCCC"
