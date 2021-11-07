import sys
from collections import defaultdict
from operator import itemgetter

import numpy as np
from Bio import SeqIO


def semiglobal_alignment(seq1, seq2, scoring_function):
    # Prepend indel to strings to pad table
    seq1, seq2 = "-" + seq1, "-" + seq2

    M = np.empty(shape=(len(seq1), len(seq2)))
    M[:, 0] = M[0, :] = 0
    MT = defaultdict(lambda: (-1, -1))
    for i in range(1, len(seq1)):
        MT[i, 0] = (i - 1, 0)
    for i in range(1, len(seq2)):
        MT[0, i] = (0, i - 1)

    for i in range(1, len(seq1)):
        for j in range(1, len(seq2)):
            M[i, j], MT[i, j] = max(
                (M[i - 1, j] + scoring_function(seq1[i], "-"), (i - 1, j)),
                (M[i, j - 1] + scoring_function("-", seq2[j]), (i, j - 1)),
                (M[i - 1, j - 1] + scoring_function(seq1[i], seq2[j]), (i - 1, j - 1)),
                key=itemgetter(0),
            )

    # Find the maximum entry in either the last row or last column
    s1_max, s2_max = np.argmax(M[:, -1]), np.argmax(M[-1, :])
    s1_max_val, s2_max_val = M[s1_max, -1], M[-1, s2_max]
    score, (i_max, j_max) = max(
        (s1_max_val, (s1_max, len(seq2) - 1)),
        (s2_max_val, (len(seq1) - 1, s2_max)),
        key=itemgetter(0),
    )

    # Trace-back and construct aligned sequences
    i, j, alignment1, alignment2 = len(seq1) - 1, len(seq2) - 1, [], []
    while i > i_max:
        alignment1.insert(0, seq1[i])
        alignment2.insert(0, "-")
        i -= 1
    while j > j_max:
        alignment1.insert(0, "-")
        alignment2.insert(0, seq2[j])
        j -= 1
    while i > 0 or j > 0:
        ni, nj = MT[i, j]
        alignment1.insert(0, seq1[i] if i != ni else "-")
        alignment2.insert(0, seq2[j] if j != nj else "-")
        i, j = ni, nj

    alignment1 = "".join(alignment1)
    alignment2 = "".join(alignment2)

    return alignment1, alignment2, score


def scoring_function(x, y):
    if x == "-" or y == "-":
        return -1
    return [-1, 1][x == y]


if __name__ == "__main__":
    records = list(SeqIO.parse(sys.stdin, format="fasta"))
    s1, s2 = records[0].seq, records[1].seq

    align1, align2, score = semiglobal_alignment(str(s1), str(s2), scoring_function)
    print(int(score), align1, align2, sep="\n")
