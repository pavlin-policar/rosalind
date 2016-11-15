import sys
from collections import deque
from operator import itemgetter


def parse_matrix(lines):
    """Parse whitespace delimited matrices."""
    M = {}
    top, *body = [l.split() for l in lines]
    for line in body:
        r, *vals = line
        for i, c in enumerate(top):
            M[r, c] = float(vals[i])
    return M


def viterbi(x, T, E, B):
    states = {si for si, _ in T.keys()}
    M = {(0, s): E[s, x[0]] * B[s] for s in states}
    TR = {(0, s): None for s in states}
    for i, c in enumerate(x):
        for s in states:
            if i is 0:
                continue
            p, o = max([(M[i - 1, si] * T[si, s], si) for si in states],
                       key=itemgetter(0))
            M[i, s], TR[i, s] = E[s, c] * p, o
    return M, TR


def traceback(T, pos):
    i, s, string = *pos, deque()
    while i >= 0:
        string.appendleft(s)
        s, i = T[i, s], i - 1
    return ''.join(string)


if __name__ == '__main__':
    lines = [l for l in sys.stdin.read().splitlines() if not l.startswith('-')]
    x = lines[0]
    alphabet = lines[1].split()
    states = lines[2].split()
    T = parse_matrix(lines[3:3 + len(states) + 1])
    E = parse_matrix(lines[3 + len(states) + 1:3 + len(states) * 2 + 2])
    B = {s: 1 / len(states) for s in states}
    matrix, trace = viterbi(x, T, E, B)
    _, end_state = max([(matrix[len(x) - 1, s], s) for s in states])
    print(traceback(trace, (len(x) - 1, end_state)))
