import os, sys
import pytest

sys.path.append(os.path.join(sys.path[0], '../src'))

from alignment import Alignment, Base, CigarUnavailableError
from test_alignment_cases import *

@pytest.mark.parametrize('string,tokens', bowtie2_cigars + hisat2_cigars)
def test_parse_cigar(string, tokens):
    assert Alignment._parse_cigar(string) == tokens

def test_parse_cigar_raises_cigar_unavailable_error():
    with pytest.raises(CigarUnavailableError):
        Alignment._parse_cigar('*')

@pytest.mark.parametrize('string,tokens', bowtie2_mds)
def test_parse_md(string, tokens):
    assert Alignment._parse_md(string) == tokens

@pytest.mark.parametrize('cigar,md', bowtie2_fields + hisat2_fields)
def test_alignment_get_cigar(cigar, md):
    alignment = Alignment(1, cigar, md)
    assert alignment.get_cigar() == cigar

@pytest.mark.parametrize('cigar,md', bowtie2_fields + hisat2_fields)
def test_alignment_get_md(cigar, md):
    alignment = Alignment(1, cigar, md)
    assert alignment.get_md() == md

@pytest.mark.parametrize('pre_clip,to_clip,post_clip', clip_test_set)
def test_alignment_soft_clip(pre_clip, to_clip, post_clip):
    alignment_original = Alignment(*pre_clip)
    alignment_original.soft_clip(*to_clip)
    assert alignment_original == Alignment(*post_clip)
