def walk_n(seq, n):
    for i in range(len(seq) - 1 - n):
        yield seq[i:i + n]


if __name__ == '__main__':
    sequence = input()
    substring = input()
    occs = [i + 1 for i, x in enumerate(walk_n(sequence, len(substring)))
            if x == substring]
    print(' '.join(str(x) for x in occs))
