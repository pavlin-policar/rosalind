from itertools import product
from collections import Counter, defaultdict

from scipy import stats

if __name__ == "__main__":
    num_generations, n = list(map(int, input().split()))
    # num_generations, n = 2, 1
    mate = ["Aa", "Bb"]

    population = {("Aa", "Bb"): 1}

    # Find the probabilities of each genotype in a given generation
    for generation in range(num_generations):
        new_population = defaultdict(float)
        for source in population:
            # Generate the different factor pairs e.g. AA, Aa, aA, aa
            combs = [
                list(map("".join, map(sorted, product(*factor_pair))))
                for factor_pair in zip(source, mate)
            ]

            organism_population = Counter(product(*combs))
            for k, v in organism_population.items():
                new_population[k] += v * population[source]

        N = sum(new_population.values())
        new_population = {k: n / N for k, n in new_population.items()}

        population = new_population

    p = population["Aa", "Bb"]
    num_children = 2 ** num_generations
    print(stats.binom(num_children, p).sf(n - 1))
