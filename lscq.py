import sys
from collections import defaultdict

from Bio import SeqIO


def print_dpt(dpt):
    max_i, max_j = -1, -1
    for i, j in dpt:
        max_i = max(max_i, i)
        max_j = max(max_j, j)

    for i in range(max_i + 1):
        for j in range(max_j + 1):
            print(f"{dpt[i, j]:2d}", end=" ")
        print()


def longest_common_subsequence(s1, s2):
    dpt = defaultdict(int)
    traceback = defaultdict(lambda: (-1, -1))

    for i in range(len(s1)):
        for j in range(len(s2)):
            dpt[i, j], traceback[i, j] = max(
                (dpt[i - 1, j], (i - 1, j)),
                (dpt[i, j - 1], (i, j - 1)),
                (dpt[i - 1, j - 1] + int(s1[i] == s2[j]), (i - 1, j - 1))
            )

    i, j = len(s1) - 1, len(s2) - 1
    lscq = []
    while i >= 0 or j >= 0:
        ni, nj = traceback[i, j]
        if ni < i and nj < j:
            lscq.insert(0, s1[i])
        i, j = ni, nj

    return "".join(lscq)


if __name__ == "__main__":
    s1, s2 = SeqIO.parse(sys.stdin, format="fasta")
    s1, s2 = str(s1.seq), str(s2.seq)

    lscq = longest_common_subsequence(s1, s2)
    print(len(s1), len(s2), len(lscq))
    print(lscq)
