import math
import sys
from collections import Counter

from Bio import SeqIO

if __name__ == "__main__":
    seq = list(SeqIO.parse(sys.stdin, format="fasta"))[0]

    counts = Counter(seq.seq)
    max_counts_au = max(counts["A"], counts["U"])
    min_counts_au = min(counts["A"], counts["U"])

    max_counts_cg = max(counts["C"], counts["G"])
    min_counts_cg = min(counts["C"], counts["G"])

    result = (
        math.perm(max_counts_au, min_counts_au) *
        math.perm(max_counts_cg, min_counts_cg)
    )
    print(result)
