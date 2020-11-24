import sys

from Bio import SeqIO


def edit_distance(s1, s2):
    dpt = {(0, 0): 0}
    traceback = {(0, 0): (-1, -1)}
    for i in range(1, len(s1)):
        dpt[i, 0] = dpt[i - 1, 0] + 1
        traceback[i, 0] = (i - 1, 0)
    for j in range(1, len(s2)):
        dpt[0, j] = dpt[0, j - 1] + 1
        traceback[0, j] = (0, j - 1)

    # Populate dynamic programming table
    for i in range(1, len(s1)):
        for j in range(1, len(s2)):
            dpt[i, j], traceback[i, j] = min(
                (dpt[i - 1, j] + 1, (i - 1, j)),  # insertion
                (dpt[i, j - 1] + 1, (i, j - 1)),  # deletion
                (dpt[i - 1, j - 1] + (s1[i] != s2[j]), (i - 1, j - 1)),
            )

    # Perform traceback
    i, j, alignment1, alignment2 = len(s1) - 1, len(s2) - 1, [], []
    score = dpt[i, j]
    while i >= 0 or j >= 0:
        ni, nj = traceback[i, j]
        alignment1.insert(0, s1[i] if ni != i else "-")
        alignment2.insert(0, s2[j] if nj != j else "-")
        i, j = ni, nj

    alignment1 = "".join(alignment1)
    alignment2 = "".join(alignment2)

    return score, alignment1, alignment2


if __name__ == "__main__":
    s1, s2 = SeqIO.parse(sys.stdin, format="fasta")
    result = edit_distance(str(s1.seq), str(s2.seq))
    print(*result, sep="\n")
