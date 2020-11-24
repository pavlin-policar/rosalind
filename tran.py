import sys
from Bio import SeqIO

if __name__ == "__main__":
    r1, r2 = SeqIO.parse(sys.stdin, format="fasta")
    s1, s2 = r1.seq, r2.seq

    transitions = transversions = 0
    for x, y in zip(s1, s2):
        if x == y:
            pass
        elif (
            x == "A" and y == "G"
            or x == "G" and y == "A"
            or x == "C" and y == "T"
            or x == "T" and y == "C"
        ):
            transitions += 1
        else:
            transversions += 1

    print(transitions / transversions)
