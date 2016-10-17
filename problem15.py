import sys

if __name__ == '__main__':
    name = input()
    content = sys.stdin.read().replace('\n', '')

    P, j = [0] * len(content), 0
    for i in range(1, len(content)):
        while j > 0 and content[i] != content[j]:
            j = P[j - 1]
        if content[i] == content[j]:
            j += 1
        P[i] = j

    print(*P)
