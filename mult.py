import math
import operator
import sys
from collections import defaultdict, deque
from itertools import combinations, product

from Bio import SeqIO


def global_alignment(*seqs, score):
    """Global multiple sequence alignment using the Needleman-Wunsch dynamic
    programming algorithm."""
    M, T = defaultdict(lambda: -math.inf), defaultdict(lambda: (-1, -1))
    # All the different ways we can calculate a cell, we skip (0, 0, ..., 0)
    _, *calc_options = product((0, 1), repeat=len(seqs))

    def score_combinations(*chars):
        return sum(score(x, y) for x, y in combinations(chars, 2))

    for pos in product(*(range(-1, len(s)) for s in seqs)):
        is_origin = not any(map(lambda x: x + 1, pos))
        M[pos], T[pos] = (0, (-1, -1)) if is_origin else max((
            (M[tuple(x - y for x, y in zip(pos, o))] +
             score_combinations(*(s[x] if y == 1 else '-'
                                  for x, y, s in zip(pos, o, seqs))),
             tuple(x - y for x, y in zip(pos, o)))
            for o in calc_options), key=operator.itemgetter(0))
    return M, T


def traceback(T, *seqs):
    """Perform traceback for global multiple sequence alignment."""
    pos = tuple(len(s) - 1 for s in seqs)
    strings = [deque() for _ in seqs]
    while any(x + 1 for x in pos):
        npos = T[pos]
        for (i, s), (j, nj) in zip(enumerate(seqs), zip(pos, npos)):
            strings[i].appendleft(s[j] if j is not nj else '-')
        pos = npos
    return [''.join(s) for s in strings]


def score_function(x, y):
    return 0 if x == y else -1


if __name__ == '__main__':
    records = SeqIO.parse(sys.stdin, format='fasta')
    seqs = [s.seq for s in records]
    M, T = global_alignment(*seqs, score=score_function)
    pos = tuple(len(s) - 1 for s in seqs)
    print(M[pos], *traceback(T, *seqs), sep='\n')
