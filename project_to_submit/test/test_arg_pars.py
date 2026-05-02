import pytest
import sys

sys.path.append("../src")

from parse_args import parse_args

def test_correct_parsing():
    # Save the real sys.argv
    old_args = sys.argv

    # Assign the made-up values
    sys.argv = ["program", "-db", "db.fasta", "-o", "out.csv",
                "-i", "file1.fastq.gz", "-k", "19", "-n", "20"]

    # Assert!
    db, out, fastq, k, n = parse_args()
    assert db == "db.fasta"
    assert out == "out.csv"
    assert fastq == ["file1.fastq.gz"]
    assert k == 19
    assert n == 20


def test_correct_parsing_without_optionals():
    # Save the real sys.argv
    old_args = sys.argv

    # Assign the made-up values
    sys.argv = ["program", "-db", "db.fasta", "-o", "out.csv",
                "-i", "file1.fastq.gz"]

    # Assert!
    db, out, fastq, k, n = parse_args()
    assert db == "db.fasta"
    assert out == "out.csv"
    assert fastq == ["file1.fastq.gz"]


def test_parse_muliple_fastq():
    # Save the real sys.argv
    old_args = sys.argv

    # Assign the made-up values
    sys.argv = ["program", "-db", "db.fasta", "-o",
                "out.csv", "-i", "file1.fastq.gz", "file2.fastq.gz"]

    dv, out, fastq, k, n = parse_args()
    assert len(fastq) == 2


def test_parse_args_missing_dashes():
    # Simulate a user typing "db file.fasta" instead of "-db file.fasta"
    sys.argv = ["program", "db", "file.fasta", "-o", "out.csv", "-i", "in.fq.gz"]

    with pytest.raises(SystemExit):
        parse_args()


def test_missing_arguments():
    # Simulate input with missing "-o" option
    sys.argv = ["program", "-db", "file.fasta", "out.csv", "-i", "in.fq.gz"]

    with pytest.raises(SystemExit):
        parse_args()


def test_parse_args_incorrect_spacing():
    # Incorrect syntax (no space after -db)
    sys.argv = ["program", "-dbfile.fasta", "-o", "out.csv", "-i", "in.fq.gz"]

    with pytest.raises(SystemExit):
        parse_args()


def test_parse_args_error_message(capsys):
    sys.argv = ["program", "wrong", "args"]
    with pytest.raises(SystemExit):
        parse_args()

    # Capture the print statements
    captured = capsys.readouterr()
    assert "Missing required arguments" in captured.out


def test_empty_program():
    sys.argv = []
    with pytest.raises(SystemExit):
        parse_args()


def test_main_flag_without_value():
    sys.argv = ["program", "-db", "-o", "out.csv",
                "-i", "file1.fastq.gz", "-k", "19", "-n", "20"]
    
    with pytest.raises(SystemExit):
        parse_args()

def test_optional_flag_without_value():
    sys.argv = ["program", "-db", "hola", "-o", "out.csv",
                "-i", "file1.fastq.gz", "-k", "19", "-n"]
    
    with pytest.raises(SystemExit):
        parse_args()

def test_input_is_gzipped():
    sys.argv = ["program", "-db", "hola", "-o", "out.csv",
                "-i", "file1.fastq", "-k", "19", "-n"]
    
    with pytest.raises(SystemExit):
        parse_args()