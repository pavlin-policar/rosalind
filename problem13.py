import sys
from collections import defaultdict
from itertools import product


if __name__ == '__main__':
    K = 4
    name = input()
    content = sys.stdin.read().replace('\n', '')
    possible_kmers = [''.join(x) for x in product('ATGC', repeat=K)]

    kmers = defaultdict(int)
    for i in range(len(content) - (K - 1)):
        kmer = content[i:i + K]
        kmers[kmer] += 1

    result = [str(kmers[key]) for key in sorted(kmers.keys())]
    print(' '.join(result))
