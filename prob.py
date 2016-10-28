import operator
from functools import reduce
from math import log10

if __name__ == '__main__':
    sequence = input()
    probabilites = map(float, input().split())

    prob = lambda l, p: p / 2 if l in ('C', 'G') else (1 - p) / 2

    logs = [log10(reduce(operator.mul, [prob(l, p) for l in sequence], 1))
            for p in probabilites]

    print(*logs)
