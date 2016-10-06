from collections import Counter


if __name__ == '__main__':
    string = input()
    occurences = Counter(list(string.split()[0]))
    print(occurences['A'], occurences['C'], occurences['G'], occurences['T'])
