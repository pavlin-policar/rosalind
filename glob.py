import sys
import math
from collections import defaultdict
from itertools import product

from Bio import SeqIO
from Bio.SubsMat.MatrixInfo import blosum62


def global_alignment(seq1, seq2, score):
    """Calculate global alignment with linear gap penalty using a given matrix.
    """
    dpt = defaultdict(lambda: -math.inf)
    for i, j in product(range(-1, len(seq1)), range(-1, len(seq2))):
        dpt[i, j] = 0 if (i, j) == (-1, -1) else max(
            dpt[i - 1, j] + score(seq1[i], '-'),
            dpt[i, j - 1] + score('-', seq2[j]),
            dpt[i - 1, j - 1] + score(seq1[i], seq2[j]))
    return dpt


def score_function(x, y):
    if x == '-' or y == '-':
        return -5
    else:
        return blosum62.get((x, y), blosum62.get((y, x)))


if __name__ == '__main__':
    records = list(SeqIO.parse(sys.stdin, format='fasta'))
    s1, s2 = records[0].seq, records[1].seq
    M = global_alignment(s1, s2, score_function)
    print(M[len(s1) - 1, len(s2) - 1])
