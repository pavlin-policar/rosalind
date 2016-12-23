import sys
from Bio.Seq import Seq


def de_bruijn(seqs):
    paths = set()
    for seq in seqs:
        paths.add((str(seq[:-1]), str(seq[1:])))
        rev = seq.reverse_complement()
        paths.add((str(rev[:-1]), str(rev[1:])))
    return paths


if __name__ == '__main__':
    seqs = [Seq(line.strip()) for line in sys.stdin]
    for x, y in sorted(de_bruijn(seqs)):
        print('(%s, %s)' % (x, y))
