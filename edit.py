import sys
from collections import deque
from itertools import product

from Bio import SeqIO
from collections import namedtuple

Cell = namedtuple('Cell', ['distance', 'origin'])
Result = namedtuple('Result', ['distance', 'strings'])


def print_table(table):
    max_x, max_y = 0, 0
    for i, j in table.keys():
        max_x, max_y = max(max_x, i), max(max_y, j)
    for row in range(max_x + 1):
        for cell in range(max_y + 1):
            if (row, cell) in table:
                value = table[row, cell]
            else:
                value = Cell('-', '-')
            print('%8s' % str(value.origin), end='')
        print()


def edit_distance(sequences, score_function=lambda x, y: -int(x != y)):
    """Needleman Wunsch algorithm for sequence alignment."""
    x, y = sequences
    # Initialization phase
    dpt = {(0, 0): Cell(0, (-1, -1))}
    for i in range(1, len(x) + 1):
        dpt[i, 0] = Cell(dpt[i - 1, 0][0] + score_function('-', x[i - 1]), (i - 1, 0))
    for j in range(1, len(y) + 1):
        dpt[0, j] = Cell(dpt[0, j - 1][0] + score_function('-', y[j - 1]), (0, j - 1))
    # Computation phase
    for i, j in product(range(1, len(x) + 1), range(1, len(y) + 1)):
        dpt[i, j] = max(
            Cell(dpt[i - 1, j].distance + score_function(x[i - 1], '-'), (i - 1, j)),
            Cell(dpt[i, j - 1].distance + score_function(y[j - 1], '-'), (i, j - 1)),
            Cell(dpt[i - 1, j - 1].distance + score_function(x[i - 1], y[j - 1]), (i - 1, j - 1)))
    # Trace back phase
    max_x, max_y = len(x), len(y)
    i, j, str1, str2 = max_x, max_y, deque(), deque()
    while i > 0 and j > 0:
        ni, nj = dpt[i, j].origin
        str1.appendleft(x[i - 1] if i != ni else '-')
        str2.appendleft(y[j - 1] if j != nj else '-')
        i, j = ni, nj

    return Result(-dpt[len(x), len(y)].distance, (''.join(str1), ''.join(str2)))


if __name__ == '__main__':
    records = SeqIO.parse(sys.stdin, format='fasta')
    result = edit_distance(str(r.seq) for r in records)
    print(result.distance)
