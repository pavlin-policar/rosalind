def partial_permutation(n, k):
    prod = 1
    for i in range(n, n - k, -1):
        prod *= i
    return prod


if __name__ == "__main__":
    n, k = list(map(int, input().split()))
    print(partial_permutation(n, k) % 1_000_000)
