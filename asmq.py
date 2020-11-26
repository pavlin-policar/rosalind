import sys

if __name__ == "__main__":
    seqs = []
    for line in sys.stdin:
        seqs.append(line.strip())

    lengths = sorted([len(s) for s in seqs])[::-1]
    total = sum(lengths)

    for length in lengths:  # from max to min
        length_sum = sum((l for l in lengths if l >= length))
        if length_sum > 0.5 * total:
            break
    n50 = length

    for length in lengths:  # from max to min
        length_sum = sum((l for l in lengths if l >= length))
        if length_sum > 0.75 * total:
            break
    n75 = length

    print(f"{n50} {n75}")
