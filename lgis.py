from collections import defaultdict

import numpy as np


def lgis(seq):
    table = defaultdict(list)
    for i in range(len(seq)):
        table[i].append(seq[i])
        for j in range(i):
            lst = table.get(j)
            if seq[i] > lst[-1] and len(lst) + 1 > len(table[i]):
                table[i] = [*lst, seq[i]]
    return max(table.values(), key=len)


if __name__ == "__main__":
    n_inputs = int(input())
    seq = np.array(list(map(int, input().strip().split(" "))))
    print(*lgis(seq))
    print(*[-i for i in lgis(-seq)])
