import sys

if __name__ == '__main__':
    with open(sys.argv[1], 'r') as f:
        for num, line in enumerate(f, 1):
            if num % 2 == 0:
                print(line, end='')
