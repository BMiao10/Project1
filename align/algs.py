
import numpy as np

class PairwiseAligner:

    """
    Parent class for :class:`align.algs.SmithWaterman` and :class:`align.algs.NeedlemanWunsch` that 
    contains shared alignment functions and variables

    :param gapB_mat: Matrix to keep track of scores for adding gaps to seqB 
    :type gapB_mat: np.ndarray 
    :param back: Matrix that keeps track of pointers for backtracking through match matrix
    :type back: np.ndarray
    :param back_A: Matrix that keeps track of pointers for backtracking through gap A matrix
    :type back_A: np.ndarray
    :param back_B: Matrix that keeps track of pointers for backtracking through gap B matrix
    :type back_B: np.ndarray
    :param opt_score: Variable to keep track of final alignment score
    :type opt_score: int
    :param seqA: Query sequence A
    :type seqA: str
    :param seqB: Query sequence B
    :type seqB: str
    :param seqA_align: Aligned version of query sequence A
    :type seqA_align: str
    :param seqB_align: Aligned version of query sequence B
    :type seqB_align: str
    :param D_open: Opening gap penalty
    :type D_open: int
    :param D_extend: Extension gap penalty
    :type D_extend: int
        
    """
    
    def __init__( self, scoring_mat ):
        
        """
        Constructor

        :param scoring_mat: Pathlike str to .mat file with substitution matrix 
        :type scoring_mat: str
            
        """

        # init alignment and gap matrices
        self.align_mat = None
        self.gapA_mat = None
        self.gapB_mat = None
        
        # init backtrace matrices
        self.back = None
        self.back_A = None
        self.back_B = None
        
        # optimal score
        self.opt_score = 0
        
        # sequences
        self.seqA = ""
        self.seqB = ""
        
        # alignments
        self.seqA_align = ""
        self.seqB_align = ""
        
        # penalties
        self.D_open = -12
        self.D_extend = -2
        
        # read in scoring matrix and create scoring dictionary
        if isinstance(scoring_mat, str): 
            self.scoring_mat, self.alphabet = self.read_scoring_file(scoring_mat)
            self.scores_dict = self.init_scoring_dict()
        elif isinstance(scoring_mat, dict):
            self.scores_dict = scoring_mat
        else: 
            self.scoring_mat = None
            self.alphabet = None
            print("set scoring matrix and alphabet manually")
            
    
    def read_scoring_file(self, scores_file):
        
        """
        Reads in file containing scoring matrix
        
        :param scores_file: Pathlike string to .mat file with substitution matrix 
        :type scores_file: str

        :return: tuple containing n x n substitution matrix and a list of n amino acids

        """
        
        # read in file for scoring matrix
        with open(scores_file, 'r') as f:
            scores_mat = [line for line in f.read().splitlines() if "#" not in line]
            
        # extract out header characters
        alphabet = [char for char in scores_mat[0] if " " not in char]
        scores_mat = scores_mat[1:]
        
        # create list with transition matrix values
        for i in range(len(scores_mat)):
            scores_mat[i] = [int(score) for score in scores_mat[i].split(" ") if len(score) > 0]
        
        # convert to nparray
        scores_mat = np.array(scores_mat)
        
        return scores_mat, alphabet

    def init_scoring_dict(self, scoring_mat=None, alphabet=None):
        """
        Read in scoring file and creating dictionary 
        
        :param scoring_mat: n x n substitution matrix
        :type scoring_mat: np.ndarray
        :param alphabet: list of n amino acids
        :type alphabet: list

        :return: dictionary that maps tuple of amino acid (aa1, aa2) to score

        """
        
        # option to update scoring matrix
        if scoring_mat is not None: 
            self.scoring_mat = scoring_mat
        
        # option to update alphabet
        if alphabet is not None: 
            self.alphabet = alphabet
        
        # check scoring matrix type and create scoring dictionary
        scores_dict = {}
        if isinstance(self.scoring_mat, np.ndarray):
            
            # step through each character and add each pair to dict
            for i in range(len(self.alphabet)):
                for j in range(len(self.alphabet)): 
                    scores_dict[self.alphabet[i],self.alphabet[j]] = self.scoring_mat[i][j]    
                    
            return scores_dict

        return "Scoring matrix invalid"
    
    def clean_sequence(self, seq, remove_unknown=False):
        """
        Makes sequence uppercase and replaces unknown characters with * (or removes them if remove_unknown=True)
        
        :param seq: sequence to clean
        :type seq: str 
        :param remove_unknown: Whether to replace or remove unknown characters in query sequence
        :type remove_unknown: bool, default=False 

        :return: cleaned sequence with all uppercase letters

        """

        # check for alphabet
        if self.alphabet is None: 
            return "Error: add an alphabet"
        
        # check for alphabet
        seq = seq.upper()
        seq = "".join([i if i in self.alphabet else "!" for i in seq])
        
        # check for alphabet
        if "!" in seq:
            if remove_unknown: seq = seq.replace("!", "*")
            else: seq = seq.replace("!", "")
        
        return seq
    
    def align(self, seqA, seqB, sw=False):
        """
        Populates scoring and backtracing matrices with alignment scores and pointers
        
        :param seqA: First query sequence
        :type seqA: str
        :param seqB: Second query sequence
        :type seqB: str
        :param sw: whether to perform Smith-Waterman alignment (default, values clipped to 0)
        :type sw: bool, default = False

        """
        # clean and store sequences
        self.seqA = self.clean_sequence(seqA)
        self.seqB = self.clean_sequence(seqB)
        
        # reset alignment
        self.seqA_align = ""
        self.seqB_align = ""
        
        # initialize first row of gapA and backtrace matrices
        # since once i is 0, only gaps can be added to A
        for j in range(self.gapA_mat.shape[1]): 
            self.gapA_mat[0][j] = self.D_open + j * self.D_extend
            self.back_A[0][j] = 1
            self.back[0][j] = 1 

        # initalize first col of gapB and back
        # since once j is 0, only gaps can be added to B
        for i in range(self.gapB_mat.shape[0]): 
            self.gapB_mat[i][0] = self.D_open + i * self.D_extend
            self.back_B[i][0] = 2
            self.back[i][0] = 2

        # initalize corner of align_mat
        self.align_mat[0][0] = 0 
        
        # begin heuristic with backtrace values stored as 
        # 0 = diag, 1 = left, 2 = up, 3 = end (sw only)
        for i in range(1, len(self.seqA) + 1): 
            for j in range(1, len(self.seqB) + 1):
                
                # update I(A) and I(A) backtrace matrix values to 
                # keep track of whether gaps should be added to seqA
                _currScores = [self.align_mat[i, j-1] + self.D_open, # open gap in seqA
                               self.gapA_mat[i, j-1] + self.D_extend] # extend gap in seqA
                    
                self.gapA_mat[i, j] = max(_currScores) #I(A) 
                self.back_A[i, j] = np.argmax(_currScores) #backtrace

                # update I(B) and I(B) backtrace matrix values to 
                # keep track of whether gaps should be added to seqB
                _currScores = [self.align_mat[i-1, j] + self.D_open, # open gap in seqB
                               -np.inf, # ignored
                               self.gapB_mat[i-1, j] + self.D_extend] # extend gap in seqB
        
                self.gapB_mat[i,j] = max(_currScores) #I(B) 
                self.back_B[i, j] = np.argmax(_currScores) #backtrace 
        
                # update M and M backtrace matrix values
                s = self.scores_dict[self.seqA[i-1], self.seqB[j-1]]
                
                _currScores = [self.align_mat[i-1, j-1] + s, #match (diag)
                               self.gapA_mat[i, j],  #gap seqA (left)
                               self.gapB_mat[i, j]]  #gap seqB (up)
                    
                self.align_mat[i,j] = max(_currScores) #M
                
                if sw and (max(_currScores) <= 0): 
                    self.back[i, j] = 3 #add marker for end for sw
                else:
                    self.back[i, j] = np.argmax(_currScores) #backtrace


    def backtrace(self, i, j, curr_back):
        """
        Backtracking steps shared by SW and NW algorithms
        
        :param i: current row value
        :type i: int
        :param j: current col value
        :type j: int
        :param curr_back: current pointer matrix to consider
        :type curr_back: np.ndarray
        
        :return: tuple (i, j) containing updated row and column values
        """

        # update the sequence
        # if pointing to M, match sequence and move diagonally
        if curr_back is self.back: 
            self.seqA_align = self.seqA[i-1] + self.seqA_align
            self.seqB_align = self.seqB[j-1] + self.seqB_align
            i -= 1
            j -= 1

        # if pointing to gapA, add a gap to sequence A,  
        # get current character in B, and move left
        elif curr_back is self.back_A:
            self.seqA_align = "-" + self.seqA_align
            self.seqB_align = self.seqB[j-1] + self.seqB_align
            j -= 1

        # if pointing to gapB, add a gap to sequence B and move up
        elif curr_back is self.back_B: 
            self.seqB_align = "-" + self.seqB_align
            self.seqA_align = self.seqA[i-1] + self.seqA_align
            i -= 1
            
        return i, j

class NeedlemanWunsch(PairwiseAligner):
    """
    Needleman-Wunsch global alignment algorithm
    """

    def __init__( self, scoring_mat):
        """
        Constructor

        :param scoring_mat: Pathlike str to .mat file with substitution matrix passed to parent :class:`align.algs.PairwiseAligner`
        :type str or matrix
        """

        PairwiseAligner.__init__(self, scoring_mat) 

    def align(self, seqA, seqB, score_only=False):

        """
        Perform Needleman-Wunch global alignment with affine gap scoring.
        Initializes matrices and calls parent align method for scoring
        
        :param seqA: First query sequence
        :type seqA: str
        :param seqB: Second query sequence
        :type seqA: str
        :param score_only: Whether to return the aligned sequences with the score (default) or just the score
        :type score_only: bool, default=False

        :return: Gapped version of seqA, gapped version of seqB, alignment score; if score_only = True is passed, only the alignment score is returned
            
        """

        # create matrices for alignment scores and gaps 
        self.align_mat = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf
        self.gapA_mat = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf
        self.gapB_mat = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf

        # create matrices for backtracing pointers
        self.back = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf
        self.back_A = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf
        self.back_B = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf
        
        # initialize other variables with parent class align
        super(NeedlemanWunsch, self).align(seqA, seqB, sw=False)
                
        return self.backtrace(score_only=score_only)
    
    def backtrace(self, score_only=False):

        """
        Backtracing through scoring matrices using parameters for Needleman-Wunch global alignment parameters
        
        :param score_only: Whether to return the aligned sequences with the score (default) or just the score
        :type score_only: bool, default=False

        :return: Gapped version of seqA, gapped version of seqB, alignment score; if score_only = True is passed, only the alignment score is returned

        """

        # Use list to keep track of what the value of the backtrace matrix points to, with
        # the index corresponding to the pointers used (0=diag=M, 1=left=IA, 2=up=IB). 
        all_mat = [self.align_mat, self.gapA_mat, self.gapB_mat]
        all_back = [self.back, self.back_A, self.back_B]
        
        # get pointer to matrix with largest optimal value
        mat_ind = np.argmax([self.align_mat[-1][-1], self.gapA_mat[-1][-1], self.gapB_mat[-1][-1]])
        
        # get matrix containing maximum score
        curr_mat = all_mat[mat_ind]
        curr_back = all_back[mat_ind]
        
        # get starting index
        i, j = curr_mat.shape[0] - 1, curr_mat.shape[1] - 1
        
        # update optimal score and return if not backtracking
        self.opt_score = int(curr_mat[i][j])
        if score_only: return self.opt_score
        
        # begin backtracing
        while (i > 0 or j > 0):
            pointer = curr_back[i, j] # store current pointer
            i, j = super(NeedlemanWunsch, self).backtrace(i, j, curr_back) # use parent backtrack function
            curr_back = all_back[int(pointer)] # use pointer to select next backtrace matrix
        
        # return the alignments and score
        return (self.seqA_align, self.seqB_align, self.opt_score)
    
class SmithWaterman(PairwiseAligner):
    """
    Smith-Waterman local alignment algorithm

    """

    def __init__( self, scoring_mat):
        """
        Constructor
        
        :param scoring_mat: Pathlike str to .mat file with substitution matrix passed to parent :class:`align.algs.PairwiseAligner`
        :type str or matrix
        """
        PairwiseAligner.__init__(self, scoring_mat) 

    def align(self, seqA, seqB, score_only=False):

        """
        Perform Smith-Waterman local alignment with affine gap scoring.
        Initializes matrices and calls parent align method for scoring
        
        :param seqA: First query sequence
        :type seqA: str
        :param seqB: Second query sequence
        :type seqA: str
        :param score_only: Whether to return the aligned sequences with the score (default) or just the score
        :type score_only: bool, default=False

        :return: Gapped version of seqA, gapped version of seqB, alignment score; if score_only = True is passed, only the alignment score is returned

        """

        # create matrices for alignment scores and gaps 
        self.align_mat = np.zeros((len(seqA) + 1, len(seqB) + 1))
        self.gapA_mat = np.zeros((len(seqA) + 1, len(seqB) + 1))
        self.gapB_mat = np.zeros((len(seqA) + 1, len(seqB) + 1))

        # create matrices for backtracing pointers
        self.back = np.zeros((len(seqA) + 1, len(seqB) + 1))
        self.back_A = np.zeros((len(seqA) + 1, len(seqB) + 1))
        self.back_B = np.zeros((len(seqA) + 1, len(seqB) + 1))
        
        # initialize other variables with parent class align
        super(SmithWaterman, self).align(seqA, seqB, sw=True)
        
        # backtracking
        return self.backtrace(score_only)
    
    def backtrace(self, score_only=False):
        """
        Backtracing through scoring matrices using parameters for Smith-Waterman local alignment parameters
        
        :param score_only: Whether to return the aligned sequences with the score (default) or just the score
        :type score_only: bool, default=False

        :return: Gapped version of seqA, gapped version of seqB, alignment score; if score_only = True is passed, only the alignment score is returned

        """

        # Use list to keep track of what the value of the backtrace matrix points to, with
        # the index corresponding to the pointers used (0=diag=M, 1=left=IA, 2=up=IB). 
        all_mat = [self.align_mat, self.gapA_mat, self.gapB_mat]
        all_back = [self.back, self.back_A, self.back_B]
        
        # get pointer to matrix with largest optimal value
        mat_ind = np.argmax([np.max(self.align_mat), np.max(self.gapA_mat), np.max(self.gapB_mat)])
        
        # get matrix containing maximum score
        curr_mat = all_mat[mat_ind]
        curr_back = all_back[mat_ind]
        
        # get index of largest value
        max_ind = np.where(curr_mat == np.amax(curr_mat))
        i, j = max_ind[0][0], max_ind[1][0]
        
        # update optimal score and return if not backtracking
        self.opt_score = int(curr_mat[i][j])
        if score_only: return self.opt_score
        
        # begin backtracing
        while (i > 0 and j > 0) and (curr_back[i,j] < 3):
            pointer = curr_back[i, j] # store current pointer
            i, j = super(SmithWaterman, self).backtrace(i, j, curr_back) # use parent backtrack function
            curr_back = all_back[int(pointer)] # use pointer to select next backtrace matrix
        
        # return the aligned sequences and score
        return (self.seqA_align, self.seqB_align, self.opt_score)

# input fasta files of alignments
def read_fasta(input_file):
    """
    Parses fasta file to retrieve fasta sequence
        
    :param input_file: Path to fasta file
    :type input_file: str

    :return: Gapped version of seqA, gapped version of seqB, alignment score; if score_only = True is passed, only the alignment score is returned

    """

    # read 
    with open(input_file, 'r') as f:
        fasta = [line for line in f.read().splitlines()]
        return fasta[0], "".join(fasta[1:])
    





