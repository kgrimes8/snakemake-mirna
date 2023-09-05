import pytest
import sys

# setting path
sys.path.append('../workflow/')

from scripts.extract_precursors import parse_gff, parse_fasta, get_seqs, read_line
from scripts.convert_to_tsv import convert_to_tsv, remove_brackets, parse_file

"""
Tests for extract_precursors.py
"""
@pytest.mark.parametrize(
    "filepath, expected_output",
    [
        (
            "c-elegans.gff3.gz",
            {'chr': 'I', 'start': '1738635', 'end': '1738733', 'strand': '-', 'attributes': 'ID=transcript:Y71G12B.29;Parent=gene:WBGene00003278;biotype=pre_miRNA;tag=Ensembl_canonical;transcript_id=Y71G12B.29'},
        ),
    ],
)
def test_parse_gff(filepath, expected_output):
    assert parse_gff(filepath)[0] == expected_output


@pytest.mark.parametrize(
    "input_data, expected_output",
    [
        (
            "1\t2\t3\t4\t5\t6\t7\t8\t9",
        {
            "chr": '1',
            "start": '4',
            "end": '5',
            "strand": '7',
            "attributes": '9',
        }),
        (
            "this\tis\ta\ttest\tstring\tfor\tthis\tfunction\t!",
        {
            "chr": 'this',
            "start": 'test',
            "end": 'string',
            "strand": 'this',
            "attributes": '!',
        }),
    ],
)
def test_read_line(input_data, expected_output):
    assert read_line(input_data) == expected_output


def test_parse_fasta():
    fasta_dict = parse_fasta('c-elegans.fa.gz')
    assert list(fasta_dict.keys()) == ['I', 'II', 'III']


def test_get_seqs():
    fasta = parse_fasta('c-elegans.fa.gz')
    gff = parse_gff('c-elegans.gff3.gz')
    seqs = get_seqs(gff, fasta)
    assert len(seqs) == 4 # 2 mirna with 2 pieces
    assert seqs[0] == '>Y71G12B.29'
    assert seqs[2] == '>Y23H5B.14'


"""
Tests for convert_to_tsv.py
"""

def test_convert_to_tsv():
    assert convert_to_tsv(['1', '2', '3', '4', '5']) == "1\t2\t3\t4\t5"


@pytest.mark.parametrize(
    "input_data, expected_output",
    [
        (
            "(no brackets here)",
            "no brackets here",
        ),
        (
            "no brackets here",
            "no brackets here",
        ),
    ],
)
def test_remove_brackets(input_data, expected_output):
    assert remove_brackets(input_data) == expected_output


def test_parse_file():
    output_list = parse_file('precursor-test.fa.gz')
    assert output_list[0] == "Y71G12B.29\tCUGCCCGCCGGCCGCUGAUAUGUCUGGUAUUCUUGGGUUUGAACUUCCAGCGUUGAACCCGCAUAUUAGACGUAUCGACGGCCGGCGGGGCAGGUAAUG\t((((((((((((((.(((((((((((((((...(((((((.(((.......)))))))))).))))))))))))))).)))))))).))))))......\t-59.80"
    assert output_list[1] == "Y23H5B.14	CGCGGGACUAGUCAAGUGUCGGCUGCAACCAGAGCAGCCGACACUUUCCGGUUUUGCGAA	(((((((((.(..(((((((((((((.......))))))))))))).).)))))))))..	-35.60"