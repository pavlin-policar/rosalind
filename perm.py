from itertools import permutations
from math import factorial

if __name__ == '__main__':
    n = int(input())
    print(factorial(n))
    for p in permutations(range(1, n + 1)):
        print(*p)
