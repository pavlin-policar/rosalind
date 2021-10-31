from itertools import product


if __name__ == "__main__":
    alphabet = input().replace(" ", "")
    alphabet_ordering = {k: i for i, k in enumerate(alphabet)}

    def compare(x):
        return "".join(map(str, map(alphabet_ordering.get, x)))

    K = int(input())
    permutations = []
    for i in range(1, K + 1):
        permutations.extend(["".join(x) for x in product(alphabet, repeat=i)])
    print(*sorted(permutations, key=compare), sep="\n")
