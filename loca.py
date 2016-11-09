import operator
import sys
from collections import defaultdict, deque
from itertools import product

from Bio import SeqIO
from Bio.SubsMat.MatrixInfo import pam250


# TODO This doesn't work... Figure out why...


def print_dpt(table, x, y):
    for i in range(x):
        for j in range(y):
            print('%2s' % table[i, j], end=' ')
        print()


def print_tbt(table, s1, s2):
    print('     ', end=' ')
    for i in range(len(s1)):
        print('%7s ' % i, end='')
    print()
    for i in range(len(s1)):
        print('%4s' % i, end='  ')
        for j in range(len(s2)):
            print('(%2s,%2s)' % table[i, j], end=' ')
        print()


def print_table(dpt, tbt, s1, s2):
    print('     ', end=' ')
    for i in range(len(s2)):
        print('%8s%3s ' % (s2[i], i), end='')
    print()
    for i in range(len(s1)):
        print('%2s %2s' % (s1[i], i), end=' ')
        for j in range(len(s2)):
            print('[%2s](%2s,%2s)' % (dpt[i, j], *tbt[i, j]), end=' ')
        print()


def local_alignment(seq1, seq2, score):
    M, T = defaultdict(int), defaultdict(lambda: (-1, -1))
    for (i, ci), (j, cj) in product(enumerate(seq1), enumerate(seq2)):
        M[i, j], T[i, j] = max(
            (0, (-1, -1)),
            (M[i - 1, j] + score(ci, '-'), (i - 1, j)),
            (M[i, j - 1] + score('-', cj), (i, j - 1)),
            (M[i - 1, j - 1] + score(ci, cj), (i - 1, j - 1)),
            key=operator.itemgetter(0))
    # print_table(M, T, seq1, seq2)
    return M, T


def traceback(seq1, seq2, T, pos):
    s1, s2 = deque(), deque()
    (x, y), (px, py) = pos, T[pos]
    while px >= 0 and py >= 0:
        if x is not px:
            s1.appendleft(seq1[x])
        if y is not py:
            s2.appendleft(seq2[y])
        x, y = px, py
        px, py = T[x, y]
    return ''.join(s1), ''.join(s2)


def score_function(x, y):
    if x == '-' or y == '-':
        return -5
    else:
        return pam250.get((x, y), pam250.get((y, x)))


if __name__ == '__main__':
    records = list(SeqIO.parse(sys.stdin, format='fasta'))
    s1, s2 = records[0].seq, records[1].seq
    # Perform local alignment
    M, T = local_alignment(s1, s2, score=score_function)
    # Find the largest alignment score in matrix
    k, v = max(M.items(), key=operator.itemgetter(1))
    # Print out the results and the traceback
    print(v, *traceback(s1, s2, T, k), sep='\n')
