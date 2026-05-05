from kmer_class import Kmer
import gzip
import sys

import pytest

sys.path.append("../src/")

from alternative_algorithm_for_testing import (
    main,
    read_genome,
    recursive_match,
    scan_genome,
)


class MockArGenes:
    def __init__(self):
        self.coverage = {">GeneA": [0, 0, 0, 0, 0, 0]}


#############################################################
######################## Fixtures ###########################
#############################################################

@pytest.fixture
def test_setup():
    return MockArGenes()


@pytest.fixture
def kmer_db(tmp_path):
    db_file = tmp_path / "genes.fasta"
    db_file.write_text(">GeneA\nAAAAGGGG\n>GeneB\nACGTACGT\n")
    return Kmer(filename=str(db_file), kmer_size=4)


@pytest.fixture
def kmer_db_file(tmp_path):
    db_file = tmp_path / "genes.fasta"
    db_file.write_text(">GeneA\nAAAAGGGG\n>GeneB\nACGTACGT\n")
    return str(db_file)


@pytest.fixture
def gzipped_fastq(tmp_path):
    raw_content = "@read1\nAAAA\n+\n!!!!\n@read2\nGGGG\n+\n!!!!\n"
    file_path = tmp_path / "reads.fastq.gz"

    with gzip.open(file_path, "wt") as f:
        f.write(raw_content)

    return str(file_path)


@pytest.fixture
def gzipped_fastq_lowercase(tmp_path):
    raw_content = "@read1\naaaa\n+\n!!!!\n@read2\nGGGG\n+\n!!!!\n"
    file_path = tmp_path / "reads.fastq.gz"

    with gzip.open(file_path, "wt") as f:
        f.write(raw_content)

    return str(file_path)


@pytest.fixture
def gzipped_fastq_no_hits(tmp_path):
    raw_content = "@read1\nATGC\n+\n!!!!\n@read2\nGCTA\n+\n!!!!\n"
    file_path = tmp_path / "reads.fastq.gz"

    with gzip.open(file_path, "wt") as f:
        f.write(raw_content)

    return str(file_path)


@pytest.fixture
def gzipped_fastq_invalid(tmp_path):
    raw_content = "@read1\nACGTX\n+\n!!!!!\n"
    file_path = tmp_path / "invalid.fastq.gz"

    with gzip.open(file_path, "wt") as f:
        f.write(raw_content)

    return str(file_path)


#############################################################
##################### recursive_match #######################
#############################################################

def test_recursive_match_correct(test_setup):
    ar_genes = test_setup
    recursive_match(">GeneA", 0, "AAAA", "AAAA", ar_genes)
    assert ar_genes.coverage[">GeneA"] == [1, 1, 1, 1, 0, 0]


def test_recursive_match_with_mismatch(test_setup):
    ar_genes = test_setup
    recursive_match(">GeneA", 0, "AAAAAA", "AAAGGA", ar_genes)
    assert ar_genes.coverage[">GeneA"] == [1, 1, 1, 0, 0, 1]


#############################################################
######################## read_genome ########################
#############################################################

def test_read_genome_correct(gzipped_fastq):
    file_list = [gzipped_fastq]
    all_reads = list(read_genome(file_list))

    assert len(all_reads) == 4
    assert all_reads[0] == "AAAA"
    assert all_reads[1] == "TTTT"
    assert all_reads[2] == "GGGG"
    assert all_reads[3] == "CCCC"


def test_read_genome_lowercase(gzipped_fastq_lowercase):
    file_list = [gzipped_fastq_lowercase]
    all_reads = list(read_genome(file_list))

    assert len(all_reads) == 4
    assert all_reads[0] == "AAAA"
    assert all_reads[1] == "TTTT"
    assert all_reads[2] == "GGGG"
    assert all_reads[3] == "CCCC"


def test_read_genome_raises_error_for_invalid_dna(gzipped_fastq_invalid):
    file_list = [gzipped_fastq_invalid]

    with pytest.raises(ValueError):
        list(read_genome(file_list))

def test_read_genome_incomplete_record(tmp_path):
    raw_content = "@read1\nACGT\n+\n"
    file_path = tmp_path / "bad.fastq.gz"
    with gzip.open(file_path, "wt") as f:
        f.write(raw_content)

    with pytest.raises(ValueError):
        list(read_genome([str(file_path)]))


def test_read_genome_bad_header(tmp_path):
    raw_content = "read1\nACGT\n+\n!!!!\n"
    file_path = tmp_path / "bad.fastq.gz"
    with gzip.open(file_path, "wt") as f:
        f.write(raw_content)

    with pytest.raises(ValueError):
        list(read_genome([str(file_path)]))


def test_read_genome_missing_plus(tmp_path):
    raw_content = "@read1\nACGT\n-\n!!!!\n"
    file_path = tmp_path / "bad.fastq.gz"
    with gzip.open(file_path, "wt") as f:
        f.write(raw_content)

    with pytest.raises(ValueError):
        list(read_genome([str(file_path)]))


def test_read_genome_quality_length_mismatch(tmp_path):
    raw_content = "@read1\nACGT\n+\n!!!\n"
    file_path = tmp_path / "bad.fastq.gz"
    with gzip.open(file_path, "wt") as f:
        f.write(raw_content)

    with pytest.raises(ValueError):
        list(read_genome([str(file_path)]))



#############################################################
######################## scan_genome ########################
#############################################################

def test_scan_genome_correct(gzipped_fastq, kmer_db):
    ar_db = kmer_db
    file_list = [gzipped_fastq]

    scan_genome(ar_db, file_list, 4)

    assert ar_db.coverage[">GeneA"][0] == 1
    assert ar_db.coverage[">GeneA"][-1] == 1
    assert ar_db.coverage[">GeneB"][0] == 0


def test_scan_genome_no_hits(gzipped_fastq_no_hits, kmer_db):
    ar_db = kmer_db
    file_list = [gzipped_fastq_no_hits]

    scan_genome(ar_db, file_list, 4)

    assert ar_db.coverage[">GeneA"] == [0] * len(ar_db.coverage[">GeneA"])


#############################################################
########################### main ############################
#############################################################

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

