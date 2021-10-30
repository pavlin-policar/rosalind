from scipy import stats


def sequence_probablity(seq, gc_content):
    prob = 1
    for l in seq:
        if l in ["G", "C"]:
            prob *= gc_content / 2
        else:
            prob *= (1 - gc_content) / 2
    return prob


if __name__ == "__main__":
    l = input().split()
    N, gc_content = int(l[0]), float(l[1])
    sequence = input()

    seq_prob = sequence_probablity(sequence, gc_content)
    print(stats.binom(N, seq_prob).sf(0))
