if __name__ == "__main__":
    generations, lifespan = tuple(map(int, input().split(" ")))
    population = {k: 0 for k in range(lifespan)}
    population[0] = 1

    for g in range(generations - 1):
        # age rabbits
        for p in range(len(population), 0, -1):
            population[p] = population[p - 1]
        # breed rabbits
        population[0] = 0
        for p in range(2, len(population)):
            population[0] += population[p]
        # kill off rabbits
        del population[lifespan]
    
    print(sum(population.values()))
