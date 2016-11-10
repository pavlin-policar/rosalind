import sys

import operator
from functools import reduce

if __name__ == '__main__':
    # Parse annoying format
    lines = [l for l in sys.stdin.read().splitlines() if not l.startswith('-')]
    x = lines[0]
    alphabet = lines[1].split()
    pi = lines[2]
    states = lines[3].split()
    r1, r2 = lines[5].split(), lines[6].split()
    # E will be the emission matrix
    E = {(alphabet[0], r1[0]): float(r1[1]),
         (alphabet[0], r2[0]): float(r2[1]),
         (alphabet[1], r1[0]): float(r1[2]),
         (alphabet[1], r2[0]): float(r2[2]),
         (alphabet[2], r1[0]): float(r1[3]),
         (alphabet[2], r2[0]): float(r2[3])}

    print(reduce(operator.mul, (E[i, j] for i, j in (zip(x, pi))), 1))
