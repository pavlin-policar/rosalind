import sys
from collections import defaultdict

from Bio.Seq import Seq


def de_bruijn(seqs, reverse_complement=False):
    paths = set()
    for seq in seqs:
        paths.add((str(seq[:-1]), str(seq[1:])))
        if reverse_complement:
            rev = seq.reverse_complement()
            paths.add((str(rev[:-1]), str(rev[1:])))
    return list(paths)


if __name__ == "__main__":
    seqs = [Seq(line.strip()) for line in sys.stdin]
    for x, y in sorted(de_bruijn(seqs, reverse_complement=True)):
        print("(%s, %s)" % (x, y))
