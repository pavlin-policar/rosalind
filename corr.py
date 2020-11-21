import sys
from Bio import SeqIO
from collections import Counter


def hamming_distance(x, y):
    return sum(xi != yi for xi, yi in zip(x, y))


if __name__ == "__main__":
    records = list(SeqIO.parse(sys.stdin, format="fasta"))
    seqs = [str(r.seq) for r in records] + [str(r.seq.reverse_complement()) for r in records]

    seen = set(str(r.seq) for r in records)
    counter = Counter(seqs)
    correct = {k: v for k, v in counter.items() if v > 1}
    wrong = {k: v for k, v in counter.items() if v <= 1 and k in seen}

    corrections = set()
    for w in wrong:
        for c in correct:
            if hamming_distance(w, c) == 1:
                corrections.add((w, c))

    for w, c in corrections:
        print(f"{w}->{c}")
