import sys

import numpy as np
from Bio import SeqIO
from Bio.SubsMat.MatrixInfo import pam250


def local_alignment(seq1, seq2, scoring_function):
    # Prepend indel to strings to pad table
    seq1, seq2 = "-" + seq1, "-" + seq2

    table = np.empty(shape=(len(seq1), len(seq2)), dtype=np.float32)
    table[:, 0] = table[0, :] = 0

    traceback = np.empty(shape=(len(seq1), len(seq2), 2), dtype=np.int64)
    traceback[:, 0] = traceback[0, :] = 0

    # Fill in dynamic programming table
    for i in range(1, len(seq1)):
        for j in range(1, len(seq2)):
            table[i, j], traceback[i, j] = max(
                (table[i - 1, j] + scoring_function(seq1[i], "-"), (i - 1, j)),
                (table[i, j - 1] + scoring_function("-", seq2[j]), (i, j - 1)),
                (table[i - 1, j - 1] + scoring_function(seq1[i], seq2[j]), (i - 1, j - 1),),
                (0, (-1, -1)),
                key=lambda x: x[0],
            )

    i, j = np.unravel_index(np.argmax(table), table.shape)
    score = table[i, j]
    alignment1, alignment2 = [], []
    while table[i, j] > 0:
        ni, nj = traceback[i, j]
        alignment1.insert(0, seq1[i] if i != ni else "-")
        alignment2.insert(0, seq2[j] if j != nj else "-")
        i, j = ni, nj

    alignment1 = "".join(alignment1)
    alignment2 = "".join(alignment2)

    return alignment1, alignment2, score


def score_function(x, y):
    if x == "-" or y == "-":
        return -5
    else:
        return pam250.get((x, y), pam250.get((y, x)))


if __name__ == "__main__":
    records = list(SeqIO.parse(sys.stdin, format="fasta"))
    s1, s2 = records[0].seq, records[1].seq
    seq1, seq2, score = local_alignment(str(s1), str(s2), score_function)
    print(int(score), seq1.replace("-", ""), seq2.replace("-", ""), sep="\n")
