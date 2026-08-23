from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np


def enumerate_blocking_sequences(adjacency_matrix):
    """
    Enumerate all sequential occupation orders on a graph via backtracking.

    At each step, a node is chosen uniformly at random among the still
    available (unblocked) nodes; the chosen node and its neighbors become
    blocked. The process ends when no nodes remain available.

    Returns a list of (sequence, probability, final_state) tuples, where
    final_state is a 0/1 tuple indicating which nodes ended up occupied.
    """
    n = len(adjacency_matrix)
    neighbors = {i: [j for j in range(n) if adjacency_matrix[i][j] != 0] for i in range(n)}

    results = []

    def backtrack(available, sequence, probability):
        if not available:
            occupied = set(sequence)
            final_state = tuple(1 if i in occupied else 0 for i in range(n))
            results.append((tuple(v + 1 for v in sequence), probability, final_state))
            return

        num_available = len(available)
        for node in sorted(available):
            next_probability = probability / num_available
            blocked_by_node = {node} | (set(neighbors[node]) & available)
            backtrack(available - blocked_by_node, sequence + [node], next_probability)

    backtrack(set(range(n)), [], 1.0)
    return results


def compute_occupied_count_distributions(sequences):
    """
    From the enumerated blocking sequences, compute two distributions over
    k = number of occupied nodes.

    - seq_dist[k]: probability of ending with exactly k occupied nodes,
      under the sequential model (uniformly random occupation order).
    - stat_dist[k]: probability of ending with exactly k occupied nodes,
      under the static model, which assumes every distinct final blocking
      state (maximal state reachable by the process) is equally likely.
    """
    state_probability = defaultdict(float)
    for _, probability, final_state in sequences:
        state_probability[final_state] += probability

    total_distinct_states = len(state_probability)

    states_by_k = defaultdict(set)
    seq_prob_by_k = defaultdict(float)

    for final_state, probability in state_probability.items():
        k = sum(final_state)
        states_by_k[k].add(final_state)
        seq_prob_by_k[k] += probability

    stat_dist = {k: len(states) / total_distinct_states for k, states in states_by_k.items()}
    seq_dist = dict(seq_prob_by_k)

    return stat_dist, seq_dist


def analyze_graph_dominance(adjacency_matrix, tolerance=1e-9):
    """
    Full analysis of the generalized unfriendly seating process on a graph.

    Given its adjacency matrix, computes the static and sequential
    occupied-count distributions, prints a table of their cumulative
    distributions, checks stochastic dominance (D_1, D_2, or neither),
    and plots both cumulative distributions (blue = static, red = sequential).
    """
    sequences = enumerate_blocking_sequences(adjacency_matrix)
    stat_dist, seq_dist = compute_occupied_count_distributions(sequences)

    ks = sorted(set(stat_dist) | set(seq_dist))
    stat_probs = [stat_dist.get(k, 0.0) for k in ks]
    seq_probs = [seq_dist.get(k, 0.0) for k in ks]

    cumulative_stat = np.cumsum(stat_probs)
    cumulative_seq = np.cumsum(seq_probs)

    # Print cumulative distribution table.
    print(f"{'k':>4} | {'Cum_stat':>12} | {'Cum_seq':>12}")
    print("-" * 34)
    for i, k in enumerate(ks):
        print(f"{k:>4} | {cumulative_stat[i]:>12.6f} | {cumulative_seq[i]:>12.6f}")

    # Check D_1 / D_2 stochastic dominance, with tolerance for floating-point noise.
    d1_holds = all(cumulative_stat[i] >= cumulative_seq[i] - tolerance for i in range(len(ks)))
    d2_holds = all(cumulative_seq[i] >= cumulative_stat[i] - tolerance for i in range(len(ks)))

    print()
    if d1_holds and d2_holds:
        print("Both cumulative distributions coincide (D_1 and D_2 both hold).")
    elif d1_holds:
        print("D_1 holds: the static model dominates the sequential model.")
    elif d2_holds:
        print("D_2 holds: the sequential model dominates the static model.")
    else:
        print("Neither D_1 nor D_2 holds: the cumulative distributions cross.")

    # Plot cumulative distributions.
    plt.figure(figsize=(10, 6))
    plt.plot(ks, cumulative_stat, label="Static model", marker='o', color="blue")
    plt.plot(ks, cumulative_seq, label="Sequential model", marker='x', color="red")
    plt.xlabel("Maximum number of occupied nodes")
    plt.ylabel("Cumulative probability")
    plt.legend()
    plt.grid(True)
    plt.show()

    return {
        "ks": ks,
        "stat_probs": stat_probs,
        "seq_probs": seq_probs,
        "cumulative_stat": cumulative_stat,
        "cumulative_seq": cumulative_seq,
    }


""""
M=[
  [0, 0, 0, 0, 1, 0, 0, 1],
  [0, 0, 0, 0, 0, 1, 0, 1],
  [0, 0, 0, 0, 0, 0, 1, 1],
  [0, 0, 0, 0, 0, 0, 1, 1],
  [1, 0, 0, 0, 0, 0, 0, 1],
  [0, 1, 0, 0, 0, 0, 0, 1],
  [0, 0, 1, 1, 0, 0, 0, 0],
  [1, 1, 1, 1, 1, 1, 0, 0]
  ]

analyze_graph_dominance(M)
"""