from functools import lru_cache


@lru_cache()
def fibonacci(n, k=1):
    return n if n < 2 else fibonacci(n - 1, k) + k * fibonacci(n - 2, k)


if __name__ == '__main__':
    n, k = map(int, input().split(' '))
    print(fibonacci(n, k))
