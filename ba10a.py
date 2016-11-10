import operator
import sys
from functools import reduce

if __name__ == '__main__':
    # Parse annoying format
    lines = [l for l in sys.stdin.read().splitlines() if not l.startswith('-')]
    seq = lines[0]
    states = lines[1].split()
    r1, r2 = lines[3].split(), lines[4].split()
    # T will be the transistion matrix
    T = {(states[0], r1[0]): float(r1[1]),
         (states[0], r2[0]): float(r1[2]),
         (states[1], r1[0]): float(r2[1]),
         (states[1], r2[0]): float(r2[2])}

    print(reduce(operator.mul, (T[c] for c in zip(seq, seq[1:])), .5))
