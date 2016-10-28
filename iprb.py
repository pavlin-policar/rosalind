def mendel(dominant, mixed, recessive):
    """Probability that a dominant allele will be expressed."""
    population = dominant + mixed + recessive

    p = (
        dominant / population,
        mixed / population * dominant / (population - 1),
        mixed / population * (mixed - 1) / (population - 1) * .75,
        mixed / population * recessive / (population - 1) * .5,
        recessive / population * dominant / (population - 1),
        recessive / population * mixed / (population - 1) * .5,
    )
    return sum(p)


if __name__ == '__main__':
    print(mendel(*map(int, input().split())))
