import operator
import sys
from collections import defaultdict

import numpy as np
from Bio import SeqIO
from Bio.Align import substitution_matrices

BLOSUM62 = substitution_matrices.load("BLOSUM62")


def global_alignment_affine_gap_scores(seq1, seq2, scoring_function, gap_start, gap_cont):
    assert gap_start <= 0, "`gap_start` must be negative!"
    assert gap_cont <= 0, "`gap_cont` must be negative!"

    seq1, seq2 = "-" + seq1, "-" + seq2

    M_gap_s1 = np.zeros(shape=(len(seq1), len(seq2)), dtype=np.float32)
    M_gap_s2 = np.zeros(shape=(len(seq1), len(seq2)), dtype=np.float32)
    M = np.zeros(shape=(len(seq1), len(seq2)), dtype=np.float32)

    M_gap_s1[:, 0] = M[:, 0] = gap_start + np.arange(len(seq1)) * gap_cont
    M_gap_s2[0, :] = M[0, :] = gap_start + np.arange(len(seq2)) * gap_cont
    M[0, 0] = 0

    MT = defaultdict(lambda: (-1, -1))
    MT_gap_s1 = defaultdict(lambda: (-1, -1))
    MT_gap_s2 = defaultdict(lambda: (-1, -1))

    # Fill in dynamic programming table
    for i in range(1, len(seq1)):
        for j in range(1, len(seq2)):
            M_gap_s1[i, j], MT_gap_s1[i, j] = max(
                ((M_gap_s1[i, j - 1] + gap_cont), (MT_gap_s1, (i, j - 1))),
                ((M[i, j - 1] + gap_start), (MT, (i, j - 1))),
                key=operator.itemgetter(0),
            )
            M_gap_s2[i, j], MT_gap_s2[i, j] = max(
                ((M_gap_s2[i - 1, j] + gap_cont), (MT_gap_s2, (i - 1, j))),
                ((M[i - 1, j] + gap_start), (MT, (i - 1, j))),
                key=operator.itemgetter(0),
            )
            M[i, j], MT[i, j] = max(
                (M_gap_s1[i, j], (MT_gap_s1, (i, j))),
                (M_gap_s2[i, j], (MT_gap_s2, (i, j))),
                (M[i - 1, j - 1] + scoring_function(seq1[i], seq2[j]), (MT, (i - 1, j - 1))),
                key=operator.itemgetter(0),
            )

    # Trace-back and construct aligned sequences
    i, j, trace_matrix, alignment1, alignment2 = len(seq1) - 1, len(seq2) - 1, MT, [], []
    while i > 0 or j > 0:
        trace_matrix, (ni, nj) = trace_matrix[i, j]
        # If ni and nj are the same, we just switched matrices -> no output
        if i == ni and j == nj:
            continue
        alignment1.insert(0, seq1[i] if i != ni else "-")
        alignment2.insert(0, seq2[j] if j != nj else "-")
        i, j = ni, nj

    alignment1 = "".join(alignment1)
    alignment2 = "".join(alignment2)
    score = M[-1, -1]

    return alignment1, alignment2, score


def score_function(x, y):
    return BLOSUM62.get((x, y), BLOSUM62.get((y, x)))


if __name__ == "__main__":
    records = list(SeqIO.parse(sys.stdin, format="fasta"))
    s1, s2 = records[0].seq, records[1].seq

    align1, align2, score = global_alignment_affine_gap_scores(
        str(s1), str(s2), score_function, -11, -1
    )
    print(int(score), align1, align2, sep="\n")
