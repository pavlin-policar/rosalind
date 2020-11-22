import sys
from collections import defaultdict

from Bio.Seq import Seq

from dbru import de_bruijn


def cyclic(seqs):
    # Create a copy since we're going to mutate the input dict
    seqs = seqs.copy()
    start = list(seqs.keys())[0]
    final = ""
    while start in seqs:
        current = seqs[start][0]
        final += current[-1]
        del seqs[start]
        start = current
    return final


if __name__ == "__main__":
    seqs = [Seq(line.strip()) for line in sys.stdin]
    edges = de_bruijn(seqs)
    graph = defaultdict(list)
    for i, j in edges:
        graph[i].append(j)
    print(cyclic(dict(graph)))
