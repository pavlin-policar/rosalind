import sys
from Bio.Seq import Seq


def de_bruijn(seqs, reverse_complement=False):
    paths = {}
    for seq in seqs:
        paths[str(seq[:-1])] = str(seq[1:])
        if reverse_complement:
            rev = seq.reverse_complement()
            paths[str(rev[:-1])] = str(rev[1:])
    return paths


if __name__ == '__main__':
    seqs = [Seq(line.strip()) for line in sys.stdin]
    for x, y in sorted(de_bruijn(seqs, reverse_complement=True).items()):
        print('(%s, %s)' % (x, y))
