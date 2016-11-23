import sys

from collections import Counter
from collections import defaultdict

from ba10c import viterbi, traceback


def parse_matrix(lines):
    """Parse whitespace delimited matrices."""
    M = {}
    top, *body = [l.split() for l in lines]
    for line in body:
        r, *vals = line
        for i, c in enumerate(top):
            M[r, c] = float(vals[i])
    return M


def print_matrix(M):
    cols, rows = [], []
    for (i, j) in M.keys():
        if j not in cols:
            cols.append(j)
        if i not in rows:
            rows.append(i)

    print('\t', end='')
    print(*sorted(cols), sep='\t')
    for i in sorted(rows):
        print(i, end='\t')
        for j in sorted(cols):
            print('%.3f' % M[i, j], end='\t')
        print()


if __name__ == '__main__':
    lines = [l for l in sys.stdin.read().splitlines() if not l.startswith('-')]
    iterations = int(lines[0])
    x = lines[1]
    states = lines[3].split()
    T = parse_matrix(lines[4:4 + len(states) + 1])
    E = parse_matrix(lines[4 + len(states) + 1:4 + len(states) * 2 + 2])
    B = {s: 1 / len(states) for s in states}

    len_seq = len(x)

    for i in range(iterations):
        M, TR = viterbi(x, T, E, B)
        _, end_state = max([(M[len(x) - 1, s], s) for s in states])
        trace = traceback(TR, (len(x) - 1, end_state))
        # Estimate transition probabilities
        combos = list(zip(trace, trace[1:]))
        t_counter, len_combos = Counter(combos), len(combos)

        occurence_count = defaultdict(int)
        for s in states:
            for (k, _), i in t_counter.items():
                if s == k:
                    occurence_count[s] += i

        for k, item in T.items():
            if k in t_counter:
                num_occ, _ = k
                T[k] = t_counter[k] / occurence_count[num_occ]
            else:
                T[k] = 0

        new_emissions = Counter(zip(trace, x))
        occurence_count = defaultdict(int)
        for s in states:
            for (k, _), i in new_emissions.items():
                if s == k:
                    occurence_count[s] += i

        for k, item in E.items():
            if k in new_emissions:
                num_occ, _ = k
                E[k] = new_emissions[k] / occurence_count[num_occ]
            else:
                E[k] = 0

    print_matrix(T)
    print('--------')
    print_matrix(E)

