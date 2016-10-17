from collections import Counter, OrderedDict


class OrderedCounter(Counter, OrderedDict):
    pass


if __name__ == '__main__':
    string = input()
    occurences = OrderedCounter(string.split())
    for key in occurences:
        print(key, occurences[key])
