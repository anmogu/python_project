import pytest
import sys

sys.path.append("../../final_code")

from kmer_clean import Kmer

@pytest.fixture
def kmer_db(tmp_path):
    # 1. Create a dummy fasta
    db_file = tmp_path / "genes.fasta"
    db_file.write_text(">GeneA\nAAAAGGGG\n>GeneB\nACGTACGT\n")
    
    # 2. Return an initialized Kmer object
    return Kmer(filename=str(db_file), kmer_size=4)

def test_kmer_initialization(kmer_db):
    # Check that it loaded both genes
    assert len(kmer_db) == 2
    assert ">GeneA" in kmer_db.seq_lookup
    assert ">GeneB" in kmer_db.seq_lookup

def test_kmer_splitting(kmer_db):
    # Check that k-mers were generated
    # GeneA (8 chars, k=4) should have (8-4+1) = 5 kmers
    assert len(kmer_db.gene_to_kmer[">GeneA"]) == 5

    # Verify known kmers are in the lookup
    assert "AAAA" in kmer_db.kmer_lookup
    assert "AAAG" in kmer_db.kmer_lookup
    assert "AAGG" in kmer_db.kmer_lookup
    assert "AGGG" in kmer_db.kmer_lookup
    assert "GGGG" in kmer_db.kmer_lookup

def test_coverage_initialization(kmer_db):
    # Coverage should be a list of 0s, same length as the sequence
    gene_a_cov = kmer_db.coverage[">GeneA"]
    assert len(gene_a_cov) == 8
    assert all(c == 0 for c in gene_a_cov)

def test_single_fasta_sequence(tmp_path):
    # Test that it can handle the last sequence (adding sequence after the loop)
    db_file = tmp_path / "genes.fasta"
    db_file.write_text(">GeneA\nAAAAGGGG\n")
    kmer_db = Kmer(filename=str(db_file), kmer_size=4)
    assert ">GeneA" in kmer_db.seq_lookup
    assert "AAAAGGGG" in kmer_db.sequence

def test_repeated_kmers_across_genes(tmp_path):
    db_file = tmp_path / "genes.fasta"

    # Both genes contain 'AAAA' at position 0
    db_file.write_text(">GeneA\nAAAAGGGG\n>GeneB\nAAAAGGGG\n")
    kmer_db = Kmer(filename=str(db_file), kmer_size=4)
    hits = kmer_db.kmer_lookup["AAAA"]
    
    # We expect 2 hits total
    assert len(hits) == 2
    
    # Extract just the header names from the hits (which are tuples like (header, pos))
    found_headers = []
    for hit in hits:
        found_headers.append(hit[0])

    # Assert that both GeneA and GeneB are in the list of hits
    assert ">GeneA" in found_headers
    assert ">GeneB" in found_headers

