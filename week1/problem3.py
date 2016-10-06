

if __name__ == '__main__':
    string = input()
    args = input()
    a, b, c, d = [int(x) for x in args.split()]

    print(string[a:b + 1], string[c:d + 1])
