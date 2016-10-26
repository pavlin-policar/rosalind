from itertools import product

if __name__ == '__main__':
    alphabet = input().replace(' ', '')
    K = int(input())
    permutations = [''.join(x) for x in product(alphabet, repeat=K)]
    print(*sorted(permutations), sep='\n')
