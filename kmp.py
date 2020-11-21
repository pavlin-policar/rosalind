import sys

from Bio import SeqIO

if __name__ == "__main__":
    records = list(SeqIO.parse(sys.stdin, format="fasta"))
    content = str(records[0].seq)

    P, j = [0] * len(content), 0
    for i in range(1, len(content)):
        while j > 0 and content[i] != content[j]:
            j = P[j - 1]
        if content[i] == content[j]:
            j += 1
        P[i] = j

    print(*P)
