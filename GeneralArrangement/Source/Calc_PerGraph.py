from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np

def validate_adjacency_matrix(adjacency_matrix):
    """
    Verifica que adjacency_matrix corresponda a un grafo simple no dirigido:
    - Es cuadrada.
    - Todas las entradas son 0 o 1.
    - La diagonal es nula (sin lazos).
    - Es simétrica (grafo no dirigido).

    Lanza ValueError con un mensaje específico si alguna condición falla.
    """
    n = len(adjacency_matrix)

    # 1. Debe ser cuadrada.
    for i, row in enumerate(adjacency_matrix):
        if len(row) != n:
            raise ValueError(
                f"La matriz no es cuadrada: la fila {i} tiene longitud "
                f"{len(row)}, se esperaba {n}."
            )

    # 2. Entradas solo 0 o 1.
    for i in range(n):
        for j in range(n):
            val = adjacency_matrix[i][j]
            if val not in (0, 1):
                raise ValueError(
                    f"Entrada no binaria en ({i}, {j}): {val}. "
                    f"Se esperaban únicamente valores 0 o 1."
                )

    # 3. Diagonal nula (sin lazos).
    for i in range(n):
        if adjacency_matrix[i][i] != 0:
            raise ValueError(
                f"La diagonal no es nula: adjacency_matrix[{i}][{i}] = "
                f"{adjacency_matrix[i][i]}. El grafo no debe tener lazos."
            )

    # 4. Simetría (grafo no dirigido).
    for i in range(n):
        for j in range(i + 1, n):
            if adjacency_matrix[i][j] != adjacency_matrix[j][i]:
                raise ValueError(
                    f"La matriz no es simétrica en ({i}, {j}) vs ({j}, {i}): "
                    f"{adjacency_matrix[i][j]} != {adjacency_matrix[j][i]}. "
                    f"Se esperaba un grafo no dirigido."
                )

    return True
def enumerate_blocking_sequences(adjacency_matrix):
    validate_adjacency_matrix(adjacency_matrix)
    n = len(adjacency_matrix)
    """
    Enumerate all sequential occupation orders on a graph via backtracking.

    At each step, a node is chosen uniformly at random among the still
    available (unblocked) nodes; the chosen node and its neighbors become
    blocked. The process ends when no nodes remain available.

    Returns a list of (sequence, probability, final_state) tuples, where
    final_state is a 0/1 tuple indicating which nodes ended up occupied.
    """
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
    #print(f"{'k':>4} | {'Cum_stat':>12} | {'Cum_seq':>12}")
    #print("-" * 34)
    #for i, k in enumerate(ks):
    #    print(f"{k:>4} | {cumulative_stat[i]:>12.6f} | {cumulative_seq[i]:>12.6f}")

    # Check D_1 / D_2 stochastic dominance, with tolerance for floating-point noise.
    d1_holds = all(cumulative_stat[i] >= cumulative_seq[i] - tolerance for i in range(len(ks)))
    d2_holds = all(cumulative_seq[i] >= cumulative_stat[i] - tolerance for i in range(len(ks)))
     
    
    if d1_holds and d2_holds:
        Conclusion ="Both cumulative distributions coincide (D_1 and D_2 both hold)."
    elif d1_holds:
        Conclusion ="D_1 holds: the sequential model dominates the static model."
    elif d2_holds:
        Conclusion ="D_2 holds: the static model dominates the sequential model."
    else:
        Conclusion ="Neither D_1 nor D_2 holds: the cumulative distributions cross."

    # Plot cumulative distributions.
    #plt.figure(figsize=(10, 6))
    #plt.plot(ks, cumulative_stat, label="Static model", marker='o', color="blue")
    #plt.plot(ks, cumulative_seq, label="Sequential model", marker='x', color="red")
    #plt.xlabel("Maximum number of occupied nodes")
    #plt.ylabel("Cumulative probability")
    #plt.legend()
    #plt.grid(True)
    #plt.show()

    return {
        "ks": ks,
        "stat_probs": stat_probs,
        "seq_probs": seq_probs,
        "cumulative_stat": cumulative_stat,
        "cumulative_seq": cumulative_seq,
        "stochastic_dominance": Conclusion
    }


