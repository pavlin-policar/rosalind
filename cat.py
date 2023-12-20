import sys
import unittest
from collections import Counter

import numpy as np
from Bio import SeqIO


def catalan_numbers(n):
    """Recursive implementation of Catalan numbers."""
    if n <= 1:
        return 1

    result = 0
    for k in range(1, n + 1):
        result += catalan_numbers(k - 1) * catalan_numbers(n - k)

    return result


def catalan_numbers_dp(n):
    """Dynamic programming implementation of Catalan numbers."""
    catalan_nums = 0 * np.ones(n + 1, dtype=int)
    catalan_nums[0] = 1  # 0th Catalan number is 1
    table = 0 * np.ones((n + 2, n + 2), dtype=int)

    for i in range(1, n + 1):
        for k in range(1, i + 1):
            table[i, k] = catalan_nums[k - 1] * catalan_nums[i - k]
        catalan_nums[i] = table[i].sum()

    return catalan_nums[n]


def count_noncrossing_matches_full(seq):
    dpt = {}

    def _count_configurations(seq):
        nonlocal dpt
        if seq in dpt:
            return dpt[seq]

        if len(seq) == 0:
            dpt[seq] = 1
            return 1
        if len(seq) % 2 != 0:
            dpt[seq] = 0
            return 0

        counts = Counter(seq)
        if counts["A"] != counts["U"] or counts["C"] != counts["G"]:
            dpt[seq] = 0
            return 0
        if len(seq) == 2:
            dpt[seq] = 1
            return 1

        num_configurations = 0
        for i in range(len(seq)):
            for j in range(i + 1, len(seq), 2):
                a, b = seq[i], seq[j]
                if frozenset({a, b}) not in {frozenset({"A", "U"}), frozenset({"C", "G"})}:
                    continue
                left = seq[i + 1:j]
                right = seq[:i] + seq[j + 1:]
                result = _count_configurations(left) * _count_configurations(right)
                num_configurations += result
        num_configurations /= (len(seq) / 2)
        dpt[seq] = num_configurations

        return num_configurations

    return _count_configurations(seq)


def count_noncrossing_matches(seq):
    dpt = {}

    def _count_configurations(seq):
        if seq in dpt:
            return dpt[seq]

        if len(seq) == 0:
            dpt[seq] = 1
            return 1
        if len(seq) % 2 != 0:
            dpt[seq] = 0
            return 0

        counts = Counter(seq)
        if counts["A"] != counts["U"] or counts["C"] != counts["G"]:
            dpt[seq] = 0
            return 0
        if len(seq) == 2:
            dpt[seq] = 1
            return 1

        num_configurations = 0
        i = 0
        for j in range(i + 1, len(seq), 2):
            a, b = seq[i], seq[j]
            if frozenset({a, b}) not in {frozenset({"A", "U"}), frozenset({"C", "G"})}:
                continue
            left = seq[i + 1:j]
            right = seq[:i] + seq[j + 1:]
            result = _count_configurations(left) * _count_configurations(right)
            num_configurations += result
        # num_configurations /= (len(seq) / 2)
        dpt[seq] = num_configurations

        return num_configurations

    return _count_configurations(seq)


class TestCatalan(unittest.TestCase):
    def test_1(self):
        self.assertEqual(count_noncrossing_matches("AU"), 1)

    def test_2(self):
        self.assertEqual(count_noncrossing_matches("AUCG"), 1)

    def test_3(self):
        self.assertEqual(count_noncrossing_matches("AUAU"), 2)

    def test_4(self):
        self.assertEqual(count_noncrossing_matches("AUUA"), 1)

    def test_5(self):
        self.assertEqual(count_noncrossing_matches("AUUACG"), 1)

    def test_6(self):
        self.assertEqual(count_noncrossing_matches("AUAUCG"), 2)

    def test_7(self):
        self.assertEqual(count_noncrossing_matches("UAGCGUGAUCAC"), 2)


if __name__ == "__main__":
    seq = list(SeqIO.parse(sys.stdin, format="fasta"))[0]

    print(count_noncrossing_matches(str(seq.seq)) % 1_000_000)
