import random


def generate_random_connected_graph(n, max_degree=3, extra_edge_probability=0.5, seed=None):
    """
    Build a random connected simple graph on n vertices with maximum
    degree max_degree.

    Construction:
      1. Build a random spanning tree respecting the degree bound: add
         vertices one at a time, attaching each new vertex to a randomly
         chosen already-connected vertex that still has spare degree
         capacity. This guarantees the result is connected. If no vertex
         has spare capacity at some point (can happen for very restrictive
         degree bounds), the whole construction restarts.
      2. Add random extra edges on top of the tree, also respecting the
         degree bound and simplicity (no repeated edges), so the result
         is not always a bare tree.

    Parameters:
        n: number of vertices
        max_degree: maximum allowed degree per vertex
        extra_edge_probability: probability of adding each valid extra
            edge in step 2 (0 = tree only, 1 = as dense as the degree
            bound allows)
        seed: optional seed for reproducibility

    Returns the adjacency matrix as a list of lists (0-indexed vertices).
    """
    if n < 1:
        raise ValueError("n must be at least 1")
    if n > 1 and max_degree < 1:
        raise ValueError("max_degree must be at least 1 to connect more than one vertex")

    rng = random.Random(seed)

    while True:
        adjacency = [[0] * n for _ in range(n)]
        degrees = [0] * n

        nodes = list(range(n))
        rng.shuffle(nodes)

        connected_nodes = [nodes[0]]
        tree_construction_failed = False

        # Step 1: random spanning tree respecting the degree bound.
        for new_node in nodes[1:]:
            available_parents = [v for v in connected_nodes if degrees[v] < max_degree]
            if not available_parents:
                tree_construction_failed = True
                break

            parent = rng.choice(available_parents)
            adjacency[parent][new_node] = adjacency[new_node][parent] = 1
            degrees[parent] += 1
            degrees[new_node] += 1
            connected_nodes.append(new_node)

        if tree_construction_failed:
            continue  # restart from scratch

        # Step 2: add random extra edges, respecting the degree bound.
        all_pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
        rng.shuffle(all_pairs)

        for i, j in all_pairs:
            if adjacency[i][j] == 1:
                continue
            if degrees[i] < max_degree and degrees[j] < max_degree:
                if rng.random() < extra_edge_probability:
                    adjacency[i][j] = adjacency[j][i] = 1
                    degrees[i] += 1
                    degrees[j] += 1

        return adjacency


# Example usage:
# random_matrix = generate_random_connected_graph(n=10, max_degree=3, seed=42)