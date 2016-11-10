import math
import sys
from collections import defaultdict
from itertools import combinations, product

from Bio import SeqIO


def global_alignment(*seqs, score):
    M = defaultdict(lambda: -math.inf)
    # All the different ways we can calculate a cell, we skip (0, 0, ..., 0)
    _, *calc_options = product((0, 1), repeat=len(seqs))

    def score_combinations(*chars):
        return sum(score(x, y) for x, y in combinations(chars, 2))

    for pos in product(*(range(-1, len(s)) for s in seqs)):
        M[pos] = 0 if not any(map(lambda x: x + 1, pos)) else max(
            [M[tuple(x - y for x, y in zip(pos, o))] +
             score_combinations(*[s[x] if y == 1 else '-'
                                  for x, y, s in zip(pos, o, seqs)])
             for o in calc_options])
    return M


def score_function(x, y):
    return 0 if x == y else -1


if __name__ == '__main__':
    records = SeqIO.parse(sys.stdin, format='fasta')
    seqs = [s.seq for s in records]
    # multiple_alignment(*[s.seq for s in records], score=score_function)
    M = global_alignment(*seqs, score=score_function)
    # pprint(M)
    pos = tuple(len(s) - 1 for s in seqs)
    # pprint(M)
    print(M[pos])
