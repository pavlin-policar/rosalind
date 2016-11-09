import math
import sys
from collections import defaultdict
from itertools import combinations, product
from pprint import pprint

from Bio import SeqIO


# TODO This isn't finished


def global_alignment(*seqs, score):
    # Initialization phase
    M = defaultdict(lambda: -math.inf)
    M[tuple(len(seqs) * [-1])] = 0
    for i, s in enumerate(seqs):
        for j, cj in enumerate(s):
            idx = [-1] * len(seqs)
            idx[i] = j
            M[tuple(idx)] = (j + 1) * score(cj, '-')
    pprint(M)

    # All the different ways we can calculate a cell, we skip (0, 0, ..., 0)
    _, *calc_options = product((0, 1), repeat=len(seqs))

    for pos in product(*(range(len(s)) for s in seqs)):
        M[pos] = max(
            [M[tuple(x - y for x, y in zip(pos, o))] +
             score(*[s[x] if y == 1 else '-' for x, y, s in zip(pos, o, seqs)])
             for o in calc_options])

    return M


def score_function(*chars):
    return sum(0 if x == y else -1 for x, y in combinations(chars, 2))


if __name__ == '__main__':
    records = SeqIO.parse(sys.stdin, format='fasta')
    seqs = [s.seq for s in records]
    # multiple_alignment(*[s.seq for s in records], score=score_function)
    M = global_alignment(*seqs, score=score_function)
    # pprint(M)
    pos = tuple(len(s) - 1 for s in seqs)
    # pprint(M)
    print(M[pos])
