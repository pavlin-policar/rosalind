import math
import operator
import sys
from collections import defaultdict
from collections import deque
from itertools import product

from Bio import SeqIO
from Bio.SubsMat.MatrixInfo import blosum62


def global_alignment(seq1, seq2, open_penalty, extend_penalty, score):
    # Initialization
    inf = lambda: -math.inf
    M, X, Y = defaultdict(inf), defaultdict(inf), defaultdict(inf)
    M[-1, -1] = 0
    X[-1, -1] = -open_penalty
    Y[-1, -1] = -open_penalty
    for i in range(len(seq1)):
        X[i, -1] = X[i - 1, -1] - extend_penalty
    for j in range(len(seq2)):
        Y[-1, j] = Y[-1, j - 1] - extend_penalty

    trace = lambda: (-1, -1)
    MT, XT, YT = defaultdict(trace), defaultdict(trace), defaultdict(trace)
    # Computation
    for (i, ci), (j, cj) in product(enumerate(seq1), enumerate(seq2)):
        M[i, j], MT[i, j] = max(
            (M[i - 1, j - 1] + score(ci, cj), ('M', (i - 1, j - 1))),
            (X[i - 1, j - 1] + score(ci, cj), ('X', (i - 1, j - 1))),
            (Y[i - 1, j - 1] + score(ci, cj), ('Y', (i - 1, j - 1))),
            key=operator.itemgetter(0))
        X[i, j], XT[i, j] = max(
            (X[i - 1, j] - extend_penalty, ('X', (i - 1, j))),
            (M[i - 1, j] - open_penalty, ('M', (i - 1, j))),
            key=operator.itemgetter(0))
        Y[i, j], YT[i, j] = max(
            (Y[i, j - 1] - extend_penalty, ('Y', (i, j - 1))),
            (M[i, j - 1] - open_penalty, ('M', (i, j - 1))),
            key=operator.itemgetter(0))
    return M, MT, X, XT, Y, YT


def score_function(x, y):
    return blosum62.get((x, y), blosum62.get((y, x)))


def traceback(seq1, seq2, mi, MT, XT, YT):
    subs1, subs2 = deque(), deque()
    x, y = len(seq1) - 1, len(seq2) - 1
    while x >= 0 and y >= 0:
        # Determine the active matrix
        if mi == 'M':
            am = MT
        elif mi == 'X':
            am = XT
        else:
            am = YT
        # Append appropriate characters to the string
        subs1.appendleft(seq1[x] if mi in ('M', 'X') else '-')
        subs2.appendleft(seq2[y] if mi in ('M', 'Y') else '-')
        # Update the matrix index,
        mi, (x, y) = am[x, y]
    return ''.join(subs1), ''.join(subs2)


if __name__ == '__main__':
    records = list(SeqIO.parse(sys.stdin, format='fasta'))
    s1, s2 = records[0].seq, records[1].seq

    M, MT, X, XT, Y, YT = global_alignment(
        s1, s2, score=score_function, open_penalty=11, extend_penalty=1)

    pos = (len(s1) - 1, len(s2) - 1)
    best_score, starting_t = max((M[pos], 'M'), (X[pos], 'X'), (Y[pos], 'Y'))
    print(best_score, *traceback(s1, s2, starting_t, MT, XT, YT), sep='\n')
