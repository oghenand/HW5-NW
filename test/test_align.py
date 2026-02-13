# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np

def test_nw_alignment():
    """
    TODO: Write your unit test for NW alignment
    using test_seq1.fa and test_seq2.fa by
    asserting that you have correctly filled out
    the your 3 alignment matrices.
    Use the BLOSUM62 matrix and a gap open penalty
    of -10 and a gap extension penalty of -1.
    """
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")
    
    #NOTE: only doing test cases for general alignment matrix (as prof. said it's sufficient to do linear penalty which I employ with gap open penalty)
    nw = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", -10, -1)
    nw_alignment_score, seqA_aligned, seqB_aligned = nw.align(seq1, seq2)

    # first check that matrix top left corner is 0,0
    assert nw._align_matrix[0,0] == 0, 'Top left corner must be of score 0!'

    # now check values of first row and first col (should be multiples of gap open)
    # first check top row and compare to expected (array of multiples of gap open)
    expected_top_row = np.arange(nw._align_matrix.shape[1]) *-10
    assert np.allclose(nw._align_matrix[0,:], expected_top_row), 'Top row must be multiples of gap open!'

    # now check leftmost col and compare to expected (array of multiples of gap open)
    expected_leftmost_col = np.arange(nw._align_matrix.shape[0]) *-10
    assert np.allclose(nw._align_matrix[:,0], expected_leftmost_col), 'First col must be multiples of gap open!'

def test_nw_backtrace():
    """
    TODO: Write your unit test for NW backtracing
    using test_seq3.fa and test_seq4.fa by
    asserting that the backtrace is correct.
    Use the BLOSUM62 matrix. Use a gap open
    penalty of -10 and a gap extension penalty of -1.
    """
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")
    
    #NOTE: only doing test cases for general alignment matrix (as prof. said it's sufficient to do linear penalty which I employ with gap open penalty)
    #NOTE: using gap penalty of -4 (as indicated by prof. on slack)
    nw = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", -4, -1)
    nw_alignment_score, seqA_aligned, seqB_aligned = nw.align(seq3, seq4)

    #NOTE: aligned sequences not provided for new test case for linear penalty -- only alignment score provided on slack!
    assert nw_alignment_score==18, "got wrong alignment score!"
    # assert seqA_aligned = "MAVHQLIRRP", "got wrong sequence A aligned!"
    # assert seqB_aligned = "M---QLIRHP", "got wrong sequence B algined!"



