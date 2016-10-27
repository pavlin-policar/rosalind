def hamming_distance(x, y):
    return sum(xi != yi for xi, yi in zip(x, y))


if __name__ == '__main__':
    print(hamming_distance(input(), input()))
