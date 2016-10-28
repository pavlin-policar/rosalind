import sys
from collections import defaultdict
from itertools import product

import operator
from Bio import SeqIO


def longest_common_substring(*strings):
    table = defaultdict(int)
    for pos in product(*(range(len(s)) for s in strings)):
        same = len(set(s[i] for s, i in zip(strings, pos))) is 1
        table[pos] = table[tuple(map(lambda n: n - 1, pos))] + 1 if same else 0
    return max(table.items(), key=operator.itemgetter(1))


if __name__ == '__main__':
    records = list(SeqIO.parse(sys.stdin, format='fasta'))
    indices, length = longest_common_substring(*[str(r.seq) for r in records])
    start, *_ = indices
    print(records[0].seq[start - length + 1:start + 1])
