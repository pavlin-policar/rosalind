from collections import defaultdict

import numpy as np


def shortest_common_supersequence(seq1, seq2):
    """
    We define the algorithm as follows:

    We want to find the shortest path through the DPT, so we can do a +1 for
    matches and -1 for any time we take a letter from one string, but not the
    other. Since mismatches are not allwed, we can only move diagonally if the
    characters match. We can do this by adding a -inf penalty on diagnoal
    mismatch movements. The rest is pretty much the same as global alignment.

    """
    # Prepend indel to strings to pad table
    seq1, seq2 = "-" + seq1, "-" + seq2

    M = np.zeros(shape=(len(seq1), len(seq2)))
    M[:, 0] = -np.arange(len(seq1))
    M[0, :] = -np.arange(len(seq2))
    MT = defaultdict(lambda: (-1 , -1))
    for i in range(1, len(seq1)):
        MT[i, 0] = (i - 1, 0)
    for j in range(1, len(seq2)):
        MT[0, j] = (0, j - 1)

    for i in range(1, len(seq1)):
        for j in range(1, len(seq2)):
            M[i, j], MT[i, j] = max(
                (M[i - 1, j] - 1, (i - 1, j)),
                (M[i, j - 1] - 1, (i, j - 1)),
                (M[i - 1, j - 1] + [-np.inf, 1][seq1[i] == seq2[j]], (i - 1, j - 1))
            )

    # Trace-back and construct aligned sequences
    i, j, superstring = len(seq1) - 1, len(seq2) - 1, []
    while i > 0 or j > 0:
        ni, nj = MT[i, j]
        if i == ni:  # right move
            superstring.insert(0, seq2[j])
        else:
            superstring.insert(0, seq1[i])
        i, j = ni, nj

    superstring = "".join(superstring)

    return superstring


if __name__ == "__main__":
    s1 = input()
    s2 = input()

    print(shortest_common_supersequence(s1, s2))
