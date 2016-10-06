import sys
from collections import Counter, defaultdict


def gc(seq):
    counter = Counter(seq)
    return (counter['G'] + counter['C']) / len(seq)


if __name__ == '__main__':
    d, last_label = defaultdict(list), None
    for line in sys.stdin:
        if line.startswith('>'):
            last_label = line.strip('>\n')
        else:
            d[last_label].append(line)

    scores = {label: gc(''.join(e.strip() for e in d[label])) for label in d}
    max_gc = max(scores, key=scores.get)

    print(max_gc)
    print(scores[max_gc] * 100)

