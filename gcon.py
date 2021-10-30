import operator
import sys

import numpy as np
from Bio import SeqIO
from Bio.SubsMat.MatrixInfo import blosum62


def global_alignment_affine_gap_scores(seq1, seq2, scoring_function, gap_start, gap_cont):
    """Global sequence alignment using the Needleman–Wunsch algorithm, but we'll
    penalize gap continuation less.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable
    gap_start: float
        The penalty for introducing a new gap.
    gap_cont: float
        The penalty for continuing a gap.

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    """
    assert gap_start <= 0, "`gap_start` must be negative!"
    assert gap_cont <= 0, "`gap_cont` must be negative!"

    seq1, seq2 = "-" + seq1, "-" + seq2

    M_gap_s1 = np.empty(shape=(len(seq1), len(seq2)), dtype=np.float32)
    M_gap_s2 = np.empty(shape=(len(seq1), len(seq2)), dtype=np.float32)
    M = np.empty(shape=(len(seq1), len(seq2)), dtype=np.float32)

    M_gap_s1[:, 0] = M[:, 0] = gap_start + np.arange(len(seq1)) * gap_cont
    M_gap_s2[0, :] = M[0, :] = gap_start + np.arange(len(seq2)) * gap_cont
    M[0, 0] = 0

    traceback = np.empty(shape=(len(seq1), len(seq2), 2), dtype=np.int64)
    traceback[:, 0] = [[i - 1, 0] for i in range(len(seq1))]
    traceback[0, :] = [[0, j - 1] for j in range(len(seq2))]

    # Fill in dynamic programming table
    for i in range(1, len(seq1)):
        for j in range(1, len(seq2)):
            M_gap_s1[i, j] = max(
                M_gap_s1[i, j - 1] + gap_cont,
                M[i, j - 1] + gap_start,
            )
            M_gap_s2[i, j] = max(
                M_gap_s2[i - 1, j] + gap_cont,
                M[i - 1, j] + gap_start,
            )
            M[i, j], traceback[i, j] = max(
                (M_gap_s1[i, j], (i, j - 1)),
                (M_gap_s2[i, j], (i - 1, j)),
                (M[i - 1, j - 1] + scoring_function(seq1[i], seq2[j]), (i - 1, j - 1)),
                key=operator.itemgetter(0),
            )

    # Trace-back and construct aligned sequences
    i, j, alignment1, alignment2 = len(seq1) - 1, len(seq2) - 1, [], []
    while i > 0 or j > 0:
        ni, nj = traceback[i, j]
        alignment1.insert(0, seq1[i] if i != ni else "-")
        alignment2.insert(0, seq2[j] if j != nj else "-")
        i, j = ni, nj

    alignment1 = "".join(alignment1)
    alignment2 = "".join(alignment2)
    score = M[-1, -1]

    return alignment1, alignment2, score


def score_function(x, y):
    return blosum62.get((x, y), blosum62.get((y, x)))


if __name__ == "__main__":
    records = list(SeqIO.parse(sys.stdin, format="fasta"))
    s1, s2 = records[0].seq, records[1].seq

    print(global_alignment_affine_gap_scores(str(s1), str(s2), score_function, -5, -0))
