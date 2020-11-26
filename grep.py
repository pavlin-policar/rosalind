import sys

from collections import defaultdict

from copy import deepcopy


def generate_prefix_graph(seqs, prefix_length):
    graph = defaultdict(list)
    edge_data = defaultdict(dict)
    node_idx = 0
    for seq in seqs:
        graph[node_idx].append(node_idx + 1)
        edge_data[node_idx, node_idx + 1]["prefix"] = seq[:prefix_length]
        edge_data[node_idx, node_idx + 1]["postfix"] = seq[-prefix_length:]
        edge_data[node_idx, node_idx + 1]["seq"] = seq
        node_idx += 2

    # Find the nodes with the same prefix
    prefix_list = defaultdict(set)
    for i, j in edge_data:
        prefix_list[edge_data[i, j]["prefix"]].add(i)
        prefix_list[edge_data[i, j]["postfix"]].add(j)

    node_mapping = {}
    for idx, prefix in enumerate(prefix_list):
        for j in prefix_list[prefix]:
            node_mapping[j] = idx

    # Collapse graph nodes that share the same prefix
    collapsed_graph = defaultdict(list)
    for node, neighbors in graph.items():
        collapsed_graph[node_mapping[node]].extend([node_mapping[j] for j in neighbors])

    # If there is any node with no outgoing edges, it will only appear once
    all_neighbors = set()
    for edges in collapsed_graph.values():
        all_neighbors.update(edges)
    for neighbor in all_neighbors:
        collapsed_graph[neighbor]  # make sure the sink has a list

    # We also need to adapt node_data to the new node numbering
    collapsed_edge_data = {}
    for (i, j), data in edge_data.items():
        collapsed_edge_data[node_mapping[i], node_mapping[j]] = data

    return dict(collapsed_graph), dict(collapsed_edge_data)


def eulerian_path(graph, node):
    num_edges = sum(len(v) for v in graph.values())
    stack = []
    paths = []

    stack.append((node, graph, []))
    while len(stack):
        node, graph, curr_path = stack.pop()
        curr_path = list(curr_path) + [node]

        has_multiple_neighbors = len(set(graph[node])) > 1
        for neighbor in set(graph[node]):
            if has_multiple_neighbors:
                graph_ = deepcopy(graph)
            else:
                graph_ = graph
            graph_[node].remove(neighbor)
            stack.append((neighbor, graph_, curr_path))

        if not len(graph[node]):
            paths.append(curr_path)

    # Remove invalid paths
    paths = [p for p in paths if len(p) == num_edges + 1]

    return paths


def decode_path_to_seq(path, edge_data):
    s = ""
    for i in range(len(path) - 1):
        edge = path[i], path[i + 1]
        cutoff = len(edge_data[edge]["seq"]) - len(edge_data[edge]["postfix"])
        s += edge_data[edge]["seq"][:cutoff]
    return s


def decode_paths_to_seqs(paths, edge_data):
    results = set()
    for p in paths:
        results.add(decode_path_to_seq(p, edge_data))
    return results


if __name__ == "__main__":
    kmers = [seq.strip() for seq in sys.stdin]
    graph, edge_data = generate_prefix_graph(kmers, len(kmers[0]) - 1)
    paths = eulerian_path(graph, 0)
    assemblies = decode_paths_to_seqs(paths, edge_data)

    assemblies = {a for a in assemblies if a.startswith(kmers[0])}
    for a in assemblies:
        print(a)
