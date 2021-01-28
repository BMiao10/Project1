
import numpy as np

class PairwiseAligner:
    """
    Main parent class containing shared alignment functions and variables
    for Smith-Waterman and Needleman-Wunsch flavored alginments
    
    Attributes
    ----------
        align_mat : np.ndarray
            Matrix to keep track of scores for matching
        gapA_mat : np.ndarray
            Matrix to keep track of scores for adding gaps to seqA
        gapB_mat : np.ndarray
            Matrix to keep track of scores for adding gaps to seqB
        back : np.ndarray
            Matrix that keeps track of pointers for backtracking through match matrix
        back_A  : np.ndarray
            Matrix that keeps track of pointers for backtracking through gap A matrix
        back_B  : np.ndarray
            Matrix that keeps track of pointers for backtracking through gap B matrix
        opt_score : int
            Variable to keep track of final alignment score
        seqA : str
            Query sequence A
        seqB : str
            Query sequence B
        seqA_align : str
            Aligned version of query sequence A
        seqB_align : str
            Aligned version of query sequence B
        D_open : int
            Opening gap penalty
        D_extend : int
            Extension gap penalty

    Methods
    -------
        read_scoring_file(scores_file)
            Reads in scoring matrix from file
        init_scoring_dict(scoring_mat=None, alphabet=None)
            Creates a dictionary of substitution values from scoring matrix & alphabet
        clean_sequence(seq, remove_unknown=False)
            Makes sequence uppercase and removes unknown characters (or replaces with *)
        align(seqA, seqB, sw=False)
            Fills scoring and backtracing matrices 
        backtrace(i, j, curr_back)
            Backtracking steps shared by NW and SW
        
    """
    
    def __init__( self, scoring_mat ):
        
        """
        Parameters
        ----------
        scoring_mat : str
            Pathlike str to .mat file with substitution matrix 
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
        
        Params:
          scores_file (str): pathlike to .mat file with substitution matrix 

        Returns:
          scores_mat (np.ndarray) : n x n substitution matrix
          alphabet (list) : list of n amino acids

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
        
        Params:
          scoring_mat (np.ndarray) : n x n substitution matrix
          alphabet (list) : list of n amino acids

        Returns:
          scores_dict (dict) : dict that maps tuple of amino acids to score ( (aa1, aa2) -> score )

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
        Sequence -> uppercase and unknown characters replaced with * (or removed if remove_unknown=True)
        
        Params:
          seq (str) : sequence to clean
          remove_unknown (bool) : whether to replace (default=True) or remove unknown characters in query sequence

        Returns:
          seq (str) : cleaned sequence 

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
        
        Params:
          seqA (string): first sequence
          seqB (string): second sequence
          sw (bool) : whether to perform Smith-Waterman alignment (default, values clipped to 0)  

        Returns:
          N/A

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
        
        Params:
          i (int) : current row value
          j (int) : current col value
          curr_back (ndarray) : current backtracking matrix to consider

        Returns:
          i, j (int, int) : updated row and column values

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
    
    Attributes
    ----------
        See PairwiseAligner for inherited attributes

    Methods
    -------
        align(seqA, seqB, sw=False)
            Initializes matrices and calls super.align()
        backtrace(i, j, curr_back)
            Backtracking steps unique to Needleman-Wunsch

        See PairwiseAligner for additional, inherited methods 
    """

    def __init__( self, scoring_mat):
        """
        Parameters
        ----------
        scoring_mat : str
            Pathlike str to .mat file with substitution matrix 
            Passed to parent __init__
        """

        PairwiseAligner.__init__(self, scoring_mat) 

    def align(self, seqA, seqB, score_only=False):

        """
        Perform Needleman-Wunch global alignment with affine gap scoring.
        Initializes matrices and calls parent align method for scoring
        
        Params:
            seqA (str): first sequence
            seqB (str): second sequence
            score_only (bool) : whether to return the aligned sequences with the score (default) or just the score

        Returns:
            (str, str, int): gapped version of seqA, gapped version of seqB, alignment score
            if score_only = True is passed, only the alignment score is returned

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
        
        Params:
            score_only (bool) : whether to return the aligned sequences with the score (default) or just the score

        Returns:
            (str, str, int): gapped version of seqA, gapped version of seqB, alignment score
            if score_only = True is passed, only the alignment score is returned

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
    
    Attributes
    ----------
        See PairwiseAligner for inherited attributes

    Methods
    -------
        align(seqA, seqB, sw=False)
            Initializes matrices and calls super.align()
        backtrace(i, j, curr_back)
            Backtracking steps unique to Smith-Waterman

        See PairwiseAligner for additional, inherited methods 
    """

    def __init__( self, scoring_mat):
        """
        Parameters
        ----------
        scoring_mat : str
            Pathlike str to .mat file with substitution matrix 
            Passed to parent __init__
        """
        PairwiseAligner.__init__(self, scoring_mat) 

    def align(self, seqA, seqB, score_only=False):

        """
        Perform Smith-Waterman local alignment with affine gap scoring.
        Initializes matrices and calls parent align method for scoring
        
        Params:
            seqA (str): first sequence
            seqB (str): second sequence
            score_only (bool) : whether to return the aligned sequences with the score (default) or just the score

        Returns:
            (str, str, int): gapped version of seqA, gapped version of seqB, alignment score
            if score_only = True is passed, only the alignment score is returned

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
        
        Params:
            score_only (bool) : whether to return the aligned sequences with the score (default) or just the score

        Returns:
            (str, str, int): gapped version of seqA, gapped version of seqB, alignment score
            if score_only = True is passed, only the alignment score is returned

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
        
    Params:
        input_file (str) : path to fasta file

    Returns:
        (str, str) : (fasta header, fasta sequence)
    """

    # read 
    with open(input_file, 'r') as f:
        fasta = [line for line in f.read().splitlines()]
        return fasta[0], "".join(fasta[1:])
    





