import sys
from collections import defaultdict
from functools import reduce

from Bio.Seq import Seq


def kmers(s, length=3, prefix=2):
    for i in range(0, len(s) - (length - 1), length - prefix):
        yield s[i:i + length]


def generate_prefix_postfix_graph(seqs, prefix_length):
    num_nodes = len(seqs)

    node_data = defaultdict(dict)
    for idx, seq in enumerate(seqs):
        node_data[idx]["prefix"] = seq[:prefix_length]
        node_data[idx]["postfix"] = seq[-prefix_length:]

    # We can speed up graph construction by remembering which nodes have the
    # same prefix
    prefix_cache = defaultdict(list)
    for idx, seq in enumerate(seqs):
        prefix_cache[node_data[idx]["prefix"]].append(idx)

    graph = defaultdict(list)
    for i in range(num_nodes):
        graph[i] = [j for j in prefix_cache[node_data[i]["postfix"]] if j != i]

    return dict(graph)


def find_cycles(graph):
    cycles = []
    available = set(graph)
    while available:
        stack = [(next(iter(available)), [])]
        while stack:
            node, cycle = stack.pop()
            cycle.append(node)
            available -= {node}

            if len(cycle) > 1 and node == cycle[0]:
                cycles.append(cycle)
            else:
                for neighbor in graph[node]:
                    if neighbor in available or neighbor == cycle[0]:
                        stack.append((neighbor, list(cycle)))

    return cycles


def decode_cyclic_path_to_seq(path, seqs, prefix_length):
    ordered = [seqs[i] for i in path]
    return reduce(lambda acc, s: acc + s[prefix_length:], ordered, "")


if __name__ == "__main__":
    seqs = []
    for line in sys.stdin:
        seqs.append(line.strip())
    seqs += [str(Seq(seq).reverse_complement()) for seq in seqs]

    for k in range(len(seqs[0]) - 1, 1, -1):
        k_kmers = set()
        for s in seqs:
            k_kmers.update(list(kmers(s, k, k - 1)))
        k_kmers = list(k_kmers)

        graph = generate_prefix_postfix_graph(k_kmers, k - 1)
        cycles = find_cycles(graph)

        # If we've found exactly two cycles, reconstruct the string
        if len(cycles) == 2:
            print(decode_cyclic_path_to_seq(cycles[0][:-1], k_kmers, k - 1))
            break
