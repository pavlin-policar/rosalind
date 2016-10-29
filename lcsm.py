import sys
from collections import defaultdict
from itertools import product

import operator
from Bio import SeqIO


def longest_common_substring_dp(*strings):
    """Dynamic programming approach to finding longest common substrings."""
    table = defaultdict(int)
    for pos in product(*(range(len(s)) for s in strings)):
        same = len(set(s[i] for s, i in zip(strings, pos))) is 1
        table[pos] = table[tuple(map(lambda n: n - 1, pos))] + 1 if same else 0
    return max(table.items(), key=operator.itemgetter(1))


def longest_common_substring(*strings):
    shortest, *others = sorted(strings)
    length = len(shortest)
    m = ''
    for i in range(length):
        for j in range(length, i + len(m), -1):
            ss = shortest[i:j]
            matched = all(ss in s for s in others)
            if matched:
                m = ss
                break
    return m


if __name__ == '__main__':
    records = list(SeqIO.parse(sys.stdin, format='fasta'))
    print(longest_common_substring(*[str(r.seq) for r in records]))
