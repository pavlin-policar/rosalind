import math
import sys
from collections import defaultdict, namedtuple
from operator import itemgetter

from Bio import SeqIO


entry = namedtuple("dpt_entry", ["score", "n_paths"])


def count_optimal_alignments(seq1, seq2, score):
    # Prepend indel to strings to pad table
    seq1, seq2 = "-" + seq1, "-" + seq2

    M = defaultdict(lambda: entry(-math.inf, 1))
    M[-1, -1] = entry(0, 1)

    for i in range(len(seq1)):
        for j in range(len(seq2)):
            candidates = [
                (M[i - 1, j].score + score(seq1[i], "-"), M[i - 1, j].n_paths),
                (M[i, j - 1].score + score("-", seq2[j]), M[i, j - 1].n_paths),
                (M[i - 1, j - 1].score + score(seq1[i], seq2[j]), M[i - 1, j - 1].n_paths),
            ]
            max_val = max(candidates, key=itemgetter(0))[0]
            max_sources = [v[1] for v in candidates if v[0] == max_val]
            M[i, j] = entry(max_val, sum(max_sources))

    last_entry = M[len(seq1) - 1, len(seq2) - 1]
    return last_entry.n_paths, last_entry.score


def score_function(x, y):
    if x == "-" or y == "y":
        return -1
    return [-1, 0][x == y]


if __name__ == "__main__":
    records = list(SeqIO.parse(sys.stdin, format="fasta"))
    s1, s2 = records[0].seq, records[1].seq

    alignments, score = count_optimal_alignments(str(s1), str(s2), score_function)
    print(alignments % 134_217_727)
