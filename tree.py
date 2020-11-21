import sys


def connected_components(num_nodes, edges):
    components = {i: {i} for i in range(1, num_nodes + 1)}
    nc_mapping = {i: i for i in range(1, num_nodes + 1)}
    component_id = num_nodes + 1

    for i, j in edges:
        components[component_id] = components[nc_mapping[i]] | components[nc_mapping[j]]
        del components[nc_mapping[i]], components[nc_mapping[j]]
        for k in components[component_id]:
            nc_mapping[k] = component_id
        component_id += 1

    return list(components.values())


if __name__ == "__main__":
    n_nodes = int(input())
    edges = set()
    for line in sys.stdin:
        edges.add(tuple(map(int, line.strip().split(" "))))

    ccs = connected_components(n_nodes, edges)
    print(ccs)
    print(len(ccs) - 1)

