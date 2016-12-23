import sys

from collections import defaultdict

from copy import deepcopy


def de_bruijn(seqs):
    nodes = defaultdict(list)
    for seq in seqs:
        nodes[seq[:-1]].append(seq[1:])
    return nodes


def path(nodes, start):
    sequences = []
    for next_node in nodes[start]:
        # Create a copy of the nodes and remove the current node
        next_nodes = deepcopy(nodes)
        for i in next_nodes[start]:
            if i == next_node:
                next_nodes[start].remove(i)
                break
        # Find the possible next sequences
        possible_seqs = path(next_nodes, next_node)
        if not len(possible_seqs):
            sequences.extend([start])
        else:
            sequences.extend([start[0] + x for x in possible_seqs])
    # Make sure all the possible paths are empty before adding to the seqs
    return sequences


if __name__ == '__main__':
    seqs = [seq.strip() for seq in sys.stdin]
    print(*path(de_bruijn(seqs), seqs[0][:-1]), sep='\n')
