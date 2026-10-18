from collections import defaultdict
from functools import lru_cache

import matplotlib.pyplot as plt
import numpy as np


def _build_neighbor_masks(adjacency_matrix):
    """Represent each node's neighborhood as a bitmask, for fast set operations."""
    n = len(adjacency_matrix)
    neighbor_masks = [0] * n
    for i in range(n):
        mask = 0
        for j in range(n):
            if adjacency_matrix[i][j] != 0:
                mask |= (1 << j)
        neighbor_masks[i] = mask
    return neighbor_masks


def compute_sequential_distribution(adjacency_matrix):
    """
    Compute the sequential model's occupied-count distribution using dynamic
    programming with memoization over the set of available nodes (as a
    bitmask), instead of enumerating all n! occupation sequences.

    The distribution of "how many more nodes end up occupied" from a given
    available set depends only on that set (its induced subgraph), not on
    how the process arrived there. This lets identical sub-problems, which
    are reached by many different sequences, be solved only once.

    Returns a dict {k: probability} of ending with exactly k occupied nodes.
    """
    n = len(adjacency_matrix)
    neighbor_masks = _build_neighbor_masks(adjacency_matrix)

    @lru_cache(maxsize=None)
    def solve(available_mask):
        if available_mask == 0:
            return ((0, 1.0),)  # tuple form so it is hashable for caching

        nodes = [i for i in range(n) if (available_mask >> i) & 1]
        num_nodes = len(nodes)

        combined = defaultdict(float)
        for node in nodes:
            blocked_mask = (1 << node) | neighbor_masks[node]
            next_mask = available_mask & ~blocked_mask
            for k, probability in solve(next_mask):
                combined[k + 1] += probability / num_nodes

        return tuple(combined.items())

    full_mask = (1 << n) - 1
    return dict(solve(full_mask))


def enumerate_maximal_independent_sets(adjacency_matrix):
    """
    Enumerate all maximal independent sets of the graph, each exactly once,
    using the Bron-Kerbosch algorithm with pivoting on the complement graph
    (a maximal independent set of G is a maximal clique of complement(G)).

    Returns a list of frozensets, one per maximal independent set.
    """
    n = len(adjacency_matrix)
    complement_neighbors = [
        {j for j in range(n) if j != i and adjacency_matrix[i][j] == 0}
        for i in range(n)
    ]

    results = []

    def bron_kerbosch(current_set, candidates, excluded):
        if not candidates and not excluded:
            results.append(frozenset(current_set))
            return

        pivot = max(candidates | excluded, key=lambda u: len(candidates & complement_neighbors[u]))
        remaining = candidates - complement_neighbors[pivot]

        for node in list(remaining):
            bron_kerbosch(
                current_set | {node},
                candidates & complement_neighbors[node],
                excluded & complement_neighbors[node],
            )
            candidates = candidates - {node}
            excluded = excluded | {node}

    bron_kerbosch(set(), set(range(n)), set())
    return results


def compute_static_distribution(adjacency_matrix):
    """
    Compute the static model's occupied-count distribution: the probability
    of size k is (number of maximal independent sets of size k) divided by
    the total number of maximal independent sets.

    Returns a dict {k: probability}.
    """
    maximal_sets = enumerate_maximal_independent_sets(adjacency_matrix)
    total = len(maximal_sets)

    counts_by_k = defaultdict(int)
    for independent_set in maximal_sets:
        counts_by_k[len(independent_set)] += 1

    return {k: count / total for k, count in counts_by_k.items()}


def verify_d1(adjacency_matrix, tolerance=1e-9, plot=True):
    """
    Check whether D_1 holds for a graph: the static model's cumulative
    distribution dominates the sequential model's, for every k.

    Uses compute_static_distribution and compute_sequential_distribution,
    both far more efficient than enumerating every occupation sequence.
    """
    stat_dist = compute_static_distribution(adjacency_matrix)
    seq_dist = compute_sequential_distribution(adjacency_matrix)

    ks = sorted(set(stat_dist) | set(seq_dist))
    stat_probs = [stat_dist.get(k, 0.0) for k in ks]
    seq_probs = [seq_dist.get(k, 0.0) for k in ks]

    cumulative_stat = np.cumsum(stat_probs)
    cumulative_seq = np.cumsum(seq_probs)

    d1_holds = all(cumulative_stat[i] >= cumulative_seq[i] - tolerance for i in range(len(ks)))

    if d1_holds:
        print("D_1 holds: the sequential model dominates the static model.")
    else:
        print("D_1 fails: the sequential model does not dominate the static model for every k.")
        for i, k in enumerate(ks):
            if cumulative_stat[i] < cumulative_seq[i] - tolerance:
                print(f"  Violated at k = {k}: cumulative_stat = {cumulative_stat[i]:.6f}, "
                      f"cumulative_seq = {cumulative_seq[i]:.6f}")

    if plot:
        plt.figure(figsize=(10, 6))
        plt.plot(ks, cumulative_stat, label="Static model", marker='o', color="blue")
        plt.plot(ks, cumulative_seq, label="Sequential model", marker='x', color="red")
        plt.xlabel("Maximum number of occupied nodes")
        plt.ylabel("Cumulative probability")
        plt.legend()
        plt.grid(True)
        plt.show()

    return {
        "d1_holds": d1_holds,
        "ks": ks,
        "cumulative_stat": cumulative_stat,
        "cumulative_seq": cumulative_seq,
    }


# Example usage:
# result = verify_d1(random_matrix)