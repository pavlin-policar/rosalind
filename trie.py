import sys
from collections import defaultdict


def build_trie(seqs):
    trie = defaultdict(set)
    edge_data = {}

    node_idx = 1
    for seq in seqs:
        node = 0  # root
        for l in seq:
            for neighbor in trie[node]:
                if edge_data[node, neighbor] == l:
                    node = neighbor
                    break
            # we haven't encoutered this path yet
            else:
                trie[node].add(node_idx)
                edge_data[node, node_idx] = l
                node = node_idx
                node_idx += 1

    return trie, edge_data


if __name__ == "__main__":
    seqs = []
    for line in sys.stdin:
        seqs.append(line.strip())

    trie, edge_data = build_trie(seqs)
    
    for i, j in edge_data:
        print(i + 1, j + 1, edge_data[i, j])

