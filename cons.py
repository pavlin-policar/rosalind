import sys
from collections import defaultdict

import operator
from Bio import SeqIO
from collections import Counter

if __name__ == '__main__':
    records = SeqIO.parse(sys.stdin, format='fasta')
    counts = defaultdict(Counter)
    for record in records:
        for i, c in enumerate(record.seq):
            counts[i].update(c)

    print(''.join([max(counts[i].items(), key=operator.itemgetter(1))[0] for i in counts]))
    for l in 'ACGT':
        print('%s:' % l, *[counts[i][l] for i in range(len(record.seq))], sep=' ')
