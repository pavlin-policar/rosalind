import sys
from collections import defaultdict
from itertools import product

from Bio import SeqIO
from Bio.SubsMat.MatrixInfo import blosum62


def global_alignment(seq1, seq2, penalty=1):
    """Calculate global alignment with linear gap penalty using a given matrix.
    """
    dpt = defaultdict(int)
    dpt.update({(i, -1): ((i + 1) * -penalty) for i in range(len(seq1))})
    dpt.update({(-1, i): ((i + 1) * -penalty) for i in range(len(seq2))})
    for (i, ci), (j, cj) in product(enumerate(seq1), enumerate(seq2)):
        dpt[i, j] = max(
            dpt[i - 1, j] - penalty,
            dpt[i, j - 1] - penalty,
            dpt[i - 1, j - 1] + blosum62.get((ci, cj), blosum62.get((cj, ci))))
    return dpt


if __name__ == '__main__':
    records = list(SeqIO.parse(sys.stdin, format='fasta'))
    s1, s2 = records[0].seq, records[1].seq
    M = global_alignment(s1, s2, penalty=5)
    print(M[len(s1) - 1, len(s2) - 1])
