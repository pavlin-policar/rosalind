import numpy as np
from scipy import stats


if __name__ == "__main__":
    n = int(input())

    for i in range(2 * n):
        print(round(np.log10(stats.binom(2 * n, 0.5).sf(i)), 5), end=" ")
