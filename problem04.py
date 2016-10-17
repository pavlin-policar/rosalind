def iter(lower, upper):
    while lower < upper:
        if lower % 2 == 1:
            yield lower
        lower += 1

if __name__ == '__main__':
    args = input()
    a, b = [int(x) for x in args.split()]
    print(sum(iter(a, b)))
