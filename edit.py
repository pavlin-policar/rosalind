import sys

from Bio import SeqIO


def edit_distance(s1, s2):
    dpt = {(0, 0): 0}
    for i in range(1, len(s1)): dpt[i, 0] = dpt[i - 1, 0] + 1
    for j in range(1, len(s2)): dpt[0, j] = dpt[0, j - 1] + 1

    for i in range(1, len(s1)):
        for j in range(1, len(s2)):
            dpt[i, j] = min(
                dpt[i - 1, j] + 1,  # insertion
                dpt[i, j - 1] + 1,  # deletion
                dpt[i - 1, j - 1] + (s1[i] != s2[j])
            )

    return dpt[len(s1) - 1, len(s2) - 1]


if __name__ == "__main__":
    s1, s2 = SeqIO.parse(sys.stdin, format="fasta")
    result = edit_distance(str(s1.seq), str(s2.seq))
    print(result)
