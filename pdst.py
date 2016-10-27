import sys
from itertools import product

from Bio import SeqIO

if __name__ == '__main__':
    records = list(SeqIO.parse(sys.stdin, 'fasta'))
    diffs = [[xi != yi for xi, yi in zip(str(x.seq), str(y.seq))]
             for x, y in product(records, repeat=2)]

    distances = [sum(diff) / len(diffs[0]) for diff in diffs]

    print(*[' '.join(map(str, distances[x:x + len(records)]))
            for x in range(0, len(records) ** 2, len(records))], sep='\n')
