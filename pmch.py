import math
import sys
from collections import Counter

from Bio import SeqIO

if __name__ == "__main__":
    seq = list(SeqIO.parse(sys.stdin, format="fasta"))[0]

    counts = Counter(seq.seq)
    assert counts["A"] == counts["U"]
    assert counts["C"] == counts["G"]

    print(math.factorial(counts["A"]) * math.factorial(counts["C"]))
