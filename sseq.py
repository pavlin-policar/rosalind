import sys
import unittest
from collections import defaultdict

from Bio import SeqIO


class TestGetTraceback(unittest.TestCase):
    def test_get_traceback_simple(self):
        traceback = {
            (-1, -1): [],
            (0, 0): [(-1, -1)],
            (1, 1): [(0, 0)],
            (2, 3): [(1, 1)],
        }
        paths = get_traceback(traceback, (2, 3))
        self.assertListEqual(paths, [[0, 1, 3]])

    def test_get_traceback_simple_multiple_branches(self):
        traceback = {
            (-1, -1): [],
            (0, 0): [(-1, -1)],
            (0, 1): [(-1, -1)],
            (1, 3): [(0, 0), (0, 1)],
            (2, 5): [(1, 3)],
            (2, 6): [(1, 3)],
        }
        paths = get_traceback(traceback, (2, 6))
        self.assertListEqual(paths, [[0, 3, 6], [1, 3, 6]])


def print_dpt(dpt):
    max_i, max_j = -1, -1
    for i, j in dpt:
        max_i = max(max_i, i)
        max_j = max(max_j, j)

    for i in range(max_i + 1):
        for j in range(max_j + 1):
            print(f"{dpt[i, j]:2d}", end=" ")
        print()


def get_traceback(traceback, seed, limit=None):
    if not len(traceback[seed]):
        return [[]]

    new_paths = []
    for new_seed in traceback[seed][:limit]:
        sub_paths = get_traceback(traceback, new_seed)
        for p in sub_paths[:limit]:
            new_paths.append(p + [seed[1]])

    return new_paths


def subsequences(string, subsequence, limit=1):
    dpt = defaultdict(int)
    dpt[-1, -1] = 1

    traceback = defaultdict(list)

    # Build up dynamic programming table
    for j in range(len(subsequence)):
        for i in range(len(string)):
            if subsequence[j] == string[i]:
                for k in range(-1, i):
                    dpt[j, i] += dpt[j - 1, k]
                    if dpt[j - 1, k] > 0:
                        traceback[j, i].append((j - 1, k))

    # Find all the end locations of the subsequences and trace back from each seed
    paths = []
    j = len(subsequence) - 1
    for i in range(len(string)):
        if dpt[j, i] > 0:
            paths.extend(get_traceback(traceback, (j, i), limit=limit))

    return paths


if __name__ == "__main__":
    seq, subseq = SeqIO.parse(sys.stdin, format="fasta")
    seq, subseq = str(seq.seq), str(subseq.seq)

    paths = subsequences(seq, subseq, limit=1)
    print(paths)
    for x in paths[0]:
        print(x + 1, end=" ")
    print()
