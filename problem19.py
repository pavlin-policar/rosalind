import sys
from itertools import product

from Bio import SeqIO

if __name__ == '__main__':
    records = list(SeqIO.parse(sys.stdin, 'fasta'))
    diffs = [[int(str(x.seq)[i] != str(y.seq)[i]) for i in range(len(x.seq))]
             for x, y in product(records, records)]

    distances = [sum(diff) / len(diffs[0]) for diff in diffs]

    print(*[' '.join(map(str, distances[x:x + len(records)]))
            for x in range(0, len(records) ** 2, len(records))], sep='\n')
