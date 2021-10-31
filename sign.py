from itertools import product, permutations

if __name__ == "__main__":
    n = int(input())

    # Number of possible permutations
    num_possible_permutations = 1
    for i in range(2 * n, 0, -2):
        num_possible_permutations *= i
    print(num_possible_permutations)

    # Generate all possible permutations
    candidates = [[i, -i] for i in range(1, n + 1)]
    candidates = list(product(*candidates))
    for c in candidates:
        perms = list(permutations(c))
        for p in perms:
            print(" ".join(map(str, p)))
