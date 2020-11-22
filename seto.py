

if __name__ == "__main__":
    size = int(input())
    a = eval(input().strip())
    b = eval(input().strip())

    full_set = set(range(1, size + 1))
    print("=" * 30)
    print(a | b)
    print(a & b)
    print(a - b)
    print(b - a)
    print(full_set - a)
    print(full_set - b)
