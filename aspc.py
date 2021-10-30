import math


if __name__ == "__main__":
    n, k = list(map(int, input().split()))
    print(sum(math.comb(n, i) for i in range(k, n + 1)) % 1_000_000)
