import pytest
from align import algs

# initial unit tests
# probably will need more
# almost any time add a function, should add unit test

# also make sure to add requirements
# probably don't change the yaml file though

# PART 1
@pytest.fixture
def some_relevant_data():
	return np.ones(10)

#import PairwiseAligner as algs

### Unit testing
def test_fasta_io():
    # testing I/O of protein fasta sequence
    fasta = algs.read_fasta("./sequences/prot-0004.fa")
    assert(fasta[0] == ">d1flp__ 1.1.1.1.2")
    assert(fasta[1][:10] == "SLEAAQKSNV")
    assert(fasta[1][30:40] == "ALFDAHDDVF")
    assert(fasta[1][-10:] == "ALMGEIEPDM")

    fasta = algs.read_fasta("./sequences/prot-0075.fa")
    assert(fasta[0] == ">d1sra__ 1.31.1.3.1")
    assert(fasta[1][:10] == "PPCLDSELTE")
    assert(fasta[1][30:40] == "EDNNLLTEKQ")
    assert(fasta[1][-10:] == "QKDIDKDLVI")

def test_scoring_matrix_io():

    ### testing I/O of scoring matrix using PAM100
    align = algs.PairwiseAligner("./scoring_matrices/PAM100.mat")

    assert(align.scoring_mat[0][0] == 4) # test corner
    assert(align.scoring_mat[-1][-1] == 1) # test corner
    assert(align.scoring_mat[2][5] == -1) # test middle
    
    assert(align.alphabet[0] == "A") # test first
    assert(align.alphabet[-1] == "*") # test last

    # test scoring dictionary
    assert(align.scores_dict[("A", "A")] == 4)
    assert(align.scores_dict[("A","T")] == 1)
    assert(align.scores_dict[("*","C")] == -9)
    
    ### testing I/O of scoring matrix using BLOSUM62
    align = algs.PairwiseAligner("./scoring_matrices/BLOSUM62.mat")
    
    # testing scoring matrix
    assert(align.scoring_mat[0][0] == 4) # test corner
    assert(align.scoring_mat[-1][-1] == 1) # test corner
    assert(align.scoring_mat[4][10] == -1) # test middle
    
    # testing alphabet
    assert(align.alphabet[0] == "A") # test first
    assert(align.alphabet[-1] == "*") # test last
    
    # test scoring dictionary
    assert(align.scores_dict[("A", "A")] == 4)
    assert(align.scores_dict[("A","T")] == 0)
    assert(align.scores_dict[("*","C")] == -4)

    ### testing I/O of scoring matrix for SW
    align = algs.SmithWaterman("./scoring_matrices/BLOSUM62.mat")
    
    # testing scoring matrix
    assert(align.scoring_mat[0][0] == 4) # test corner
    assert(align.scoring_mat[-1][-1] == 1) # test corner
    assert(align.scoring_mat[4][10] == -1) # test middle
    
    # testing alphabet
    assert(align.alphabet[0] == "A") # test first
    assert(align.alphabet[-1] == "*") # test last
    
    # test scoring dictionary
    assert(align.scores_dict[("A", "A")] == 4)
    assert(align.scores_dict[("A","T")] == 0)
    assert(align.scores_dict[("*","C")] == -4)

    ### testing I/O of scoring matrix for NW
    align = algs.SmithWaterman("./scoring_matrices/BLOSUM62.mat")
    
    # testing scoring matrix
    assert(align.scoring_mat[0][0] == 4) # test corner
    assert(align.scoring_mat[-1][-1] == 1) # test corner
    assert(align.scoring_mat[4][10] == -1) # test middle
    
    # testing alphabet
    assert(align.alphabet[0] == "A") # test first
    assert(align.alphabet[-1] == "*") # test last
    
    # test scoring dictionary
    assert(align.scores_dict[("A", "A")] == 4)
    assert(align.scores_dict[("A","T")] == 0)
    assert(align.scores_dict[("*","C")] == -4)


def test_identical():
    # test all algorithms for ability to align all identical sequences

    # test nw
    nw = algs.NeedlemanWunsch("./scoring_matrices/BLOSUM62.mat")
    assert (nw.align("AAAAA", "AAAAA") == ('AAAAA', 'AAAAA', 20))
    assert (nw.align("AAAAA", "GGGGG") == ('AAAAA', 'GGGGG', 0))
    
    # test sw
    sw = algs.SmithWaterman("./scoring_matrices/BLOSUM62.mat")
    assert (sw.align("AAAAA", "AAAAA") == ('AAAAA', 'AAAAA', 20))
    assert (sw.align("AAAAA", "GGGGG") == ('', '', 0))

def test_alignment_score():
    # test all algorithms for ability to align all sequences
    nw = algs.NeedlemanWunsch("./scoring_matrices/BLOSUM62.mat")
    assert (nw.align("CCTT", "CCCCTTTT") == ('CC----TT', 'CCCCTTTT', 5)) or ('C----CTT', 'CCCCTTTT', 5)
    assert (nw.align("ABCDEFG", "HIKLMNP") == ('ABCDEFG', 'HIKLMNP', -19))

    # test sw
    sw = algs.SmithWaterman("./scoring_matrices/BLOSUM62.mat")
    assert (sw.align("CCTT", "CCCCTTTT") == ('CCTT', 'CCTT', 28))
    assert (sw.align("ABCDEFG", "HIKLMNP") == ('B', 'N', 2))







