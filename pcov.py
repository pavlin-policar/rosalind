import sys

from Bio.Seq import Seq

from dbru import de_bruijn


def cyclic(seqs):
    # Create a copy since we're going to mutate the input dict
    seqs = seqs.copy()
    start = list(seqs.keys())[0]
    final = ''
    while start in seqs:
        current = seqs[start]
        final += current[-1]
        del seqs[start]
        start = current
    return final


if __name__ == '__main__':
    seqs = [Seq(line.strip()) for line in sys.stdin]
    print(cyclic(de_bruijn(seqs)))
