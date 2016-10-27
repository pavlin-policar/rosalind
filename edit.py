import sys
from itertools import product

from Bio import SeqIO


def print_table(table):
    max_x, max_y = 0, 0
    for i, j in table.keys():
        max_x, max_y = max(max_x, i), max(max_y, j)
    for row in range(max_x + 1):
        for cell in range(max_y + 1):
            if (row, cell) in table:
                value = table[row, cell]
            else:
                value = '-'
            print('%4s' % value, end='')
        print()


def edit_distance(sequences, score_function=lambda x, y: x != y, optimize=min):
    """Needleman Wunsch algorithm for sequence alignment."""
    x, y = sequences
    # Initialization phase
    dpt = {(0, 0): 0}
    for i in range(1, len(x) + 1):
        dpt[i, 0] = dpt[i - 1, 0] + score_function('-', x[i - 1])
    for j in range(1, len(y) + 1):
        dpt[0, j] = dpt[0, j - 1] + score_function('-', y[j - 1])
    # Computation phase
    for i, j in product(range(1, len(x) + 1), range(1, len(y) + 1)):
        dpt[i, j] = optimize(
            dpt[i - 1, j] + score_function(x[i - 1], '-'),
            dpt[i, j - 1] + score_function(y[j - 1], '-'),
            dpt[i - 1, j - 1] + score_function(x[i - 1], y[j - 1]))
    return dpt[len(x), len(y)]


if __name__ == '__main__':
    records = SeqIO.parse(sys.stdin, format='fasta')
    print(edit_distance(str(r.seq) for r in records))
