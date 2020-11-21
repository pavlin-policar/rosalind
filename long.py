import sys
import time
import unittest
from functools import reduce
from itertools import permutations
import numpy as np

from Bio import SeqIO


class Tests(unittest.TestCase):
    def test_join_overlapping(self):
        s1 = "ACTGGAT"
        s2 = "GATTTAATT"
        self.assertEqual(join_overlapping(s1, s2), ("ACTGGATTTAATT", 3))

        s1 = "ATGATGATG"
        s2 = "ATGATGTTT"
        self.assertEqual(join_overlapping(s1, s2), ("ATGATGATGTTT", 6))

        s1 = "AAAAAAA"
        s2 = "TTTT"
        self.assertEqual(join_overlapping(s1, s2), ("AAAAAAATTTT", 0))

        s1 = "TTAAAAAAAT"
        s2 = "TAAA"
        self.assertEqual(join_overlapping(s1, s2), ("TTAAAAAAAT", 4))

    def test_join_overlapping_multiple(self):
        strings = [
            "ACTGGAT",
            "GATTTAATT",
            "AATTA",
            "AATTAT",
        ]
        self.assertEqual(join_overlapping_multiple(strings), ("ACTGGATTTAATTAT", 12))

    def test_find_shortest_superstring(self):
        strings = [
            'this',
            'isaiah',
            'mourn',
            'ned',
            'dehou',
            'house',
            'used',
            'sed',
        ]
        self.assertEqual(shortest_superstring(strings), "thisaiahmournedehoused")


def join_overlapping_naive(s1, s2):
    """Naive implementation of shortest superstring generation."""
    if s1 in s2:
        return s2, len(s1)
    if s2 in s1:
        return s1, len(s2)

    best_str = s1 + s2
    for i in range(min(len(s1), len(s2)), -1, -1):
        if s1[-i:] == s2[:i]:
            best_str = s1 + s2[i:]
            break

    return best_str, i


def join_overlapping_dpt(seq1, seq2):
    """Dynamic programming approach to shortest superstring generation."""
    seq1, seq2 = "*" + seq1, "*" + seq2

    table = np.zeros(shape=(len(seq2), len(seq1)), dtype=np.int32)
    table[0, :] = 1

    for i in range(1, len(seq2)):
        for j in range(len(seq1) - 1, i - 1, -1):
            if table[i - 1, j - 1] > 0 and seq2[i] == seq1[j]:
                table[i, j] = table[i - 1, j - 1] + 1
            else:
                table[i, j] = 0

    # Try to find a postfix-prefix match
    col_idx_max = np.argmax(table[:, -1])
    prefix_score = table[col_idx_max, -1]

    # Check if seq2 is entirely contained within seq1
    row_idx_max = np.argmax(table[-1, :])
    containment_score = table[-1, row_idx_max]

    if prefix_score > 0 and prefix_score > containment_score:
        return seq1[1:] + seq2[col_idx_max + 1:]
    elif containment_score > 0 and containment_score >= prefix_score:
        return seq1[1:]

    # Otherwise, we have no overlap
    return seq1[1:] + seq2[1:]


join_overlapping = join_overlapping_naive


def join_overlapping_multiple(fragments):
    def _helper(acc, x):
        merged, overlap = join_overlapping(acc[0], x)
        return merged, acc[1] + overlap
    return reduce(_helper, fragments, ("", 0))


def shortest_superstring_bf(fragments):
    shortest = reduce(lambda x, r: r + x, fragments, "")
    for perm in permutations(seqs):
        new = join_overlapping_multiple(perm)
        if len(new) < len(shortest):
            shortest = new

    return shortest


def shortest_superstring_it(fragments):
    seqs = list(fragments)

    while len(seqs) > 1:
        curr_best_match = None
        curr_best_overlap = None
        best_i, best_j = -1, -1
        for i in range(len(seqs)):
            for j in range(len(seqs)):
                if i == j:
                    continue
                candidate, overlap = join_overlapping(seqs[i], seqs[j])
                if curr_best_overlap is None or curr_best_overlap < overlap:
                    curr_best_match = candidate
                    curr_best_overlap = overlap
                    best_i, best_j = i, j

        # Delete multiple indices from list at once
        for idx in sorted([best_i, best_j], reverse=True):
            del seqs[idx]
        seqs.append(curr_best_match)

    return seqs[0]


shortest_superstring = shortest_superstring_it


if __name__ == "__main__":
    records = list(SeqIO.parse(sys.stdin, format="fasta"))
    seqs = [str(r.seq) for r in records]

    # Brute force -- too slow
    start = time.time()
    shortest = shortest_superstring_bf(seqs)
    print(time.time() - start)

    print(len(shortest))
    print("=" * 80)

    # Iterative matching
    start = time.time()
    shortest = shortest_superstring(seqs)
    print(time.time() - start)

    print(len(shortest))
    print(shortest)
