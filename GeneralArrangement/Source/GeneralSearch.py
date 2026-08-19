import subprocess
from Calc_PerGraph import enumerate_blocking_sequences, compute_occupied_count_distributions


def _graph6_to_matrix(g6):
    """Decode a graph6-format string into an adjacency matrix."""
    data = g6.strip()
    bits = ''.join(f'{(ord(c) - 63):06b}' for c in data)

    if ord(data[0]) - 63 < 63:
        n = ord(data[0]) - 63
        edge_bits = bits[6:]
    else:
        n = int(bits[6:24], 2)
        edge_bits = bits[24:]

    matrix = [[0] * n for _ in range(n)]
    idx = 0
    for j in range(1, n):
        for i in range(j):
            if idx < len(edge_bits) and edge_bits[idx] == '1':
                matrix[i][j] = matrix[j][i] = 1
            idx += 1

    return matrix


def generate_symmetric_matrices(n, max_degree=None, min_degree=None, connected=True):
    """
    Generate all non-isomorphic simple graphs on n vertices (as adjacency
    matrices) using nauty-geng, optionally bounded by max/min degree.

    Parameters:
        n: number of vertices
        max_degree: if given, restrict to graphs with maximum degree <= max_degree
        min_degree: if given, restrict to graphs with minimum degree >= min_degree
        connected: if True (default), restrict to connected graphs

    This is a generator: matrices are decoded and yielded one at a time,
    so a caller can stop early (e.g. on the first match) without paying
    the cost of decoding the remaining graphs.
    """
    if n < 1:
        raise ValueError("n must be at least 1")

    command = ['nauty-geng', '-q']
    if connected:
        command.append('-c')
    if max_degree is not None:
        command.append(f'-D{max_degree}')
    if min_degree is not None:
        command.append(f'-d{min_degree}')
    command.append(str(n))

    result = subprocess.run(command, capture_output=True, text=True, check=True)

    for line in result.stdout.strip().split('\n'):
        if line:
            yield _graph6_to_matrix(line)


def find_first_dominance_failure(n, mode, max_degree=None, min_degree=None,
                                  connected=True, tolerance=1e-9):
    """
    Scan all graphs on n vertices (optionally degree-bounded) and return the
    first one that fails the requested stochastic dominance property.

    mode:
        0 - first graph where D_1 fails (static does NOT dominate sequential)
        1 - first graph where D_2 fails (sequential does NOT dominate static)
        2 - first graph where NEITHER D_1 nor D_2 holds (a genuine crossing)

    Returns a dict with the adjacency matrix, number of vertices, max degree
    and min degree of the first matching graph, or None if no graph in the
    generated set matches.
    """
    if mode not in (0, 1, 2):
        raise ValueError("mode must be 0, 1, or 2")

    for adjacency_matrix in generate_symmetric_matrices(n, max_degree, min_degree, connected):
        sequences = enumerate_blocking_sequences(adjacency_matrix)
        stat_dist, seq_dist = compute_occupied_count_distributions(sequences)

        ks = sorted(set(stat_dist) | set(seq_dist))
        if not ks:
            continue

        cumulative_stat = 0.0
        cumulative_seq = 0.0
        d1_holds = True
        d2_holds = True

        for k in ks:
            cumulative_stat += stat_dist.get(k, 0.0)
            cumulative_seq += seq_dist.get(k, 0.0)
            if cumulative_stat < cumulative_seq - tolerance:
                d1_holds = False
            if cumulative_seq < cumulative_stat - tolerance:
                d2_holds = False

        failure_found = (
            (mode == 0 and not d1_holds) or
            (mode == 1 and not d2_holds) or
            (mode == 2 and not d1_holds and not d2_holds)
        )

        if failure_found:
            degrees = [sum(row) for row in adjacency_matrix]
            return {
                "adjacency_matrix": adjacency_matrix,
                "num_vertices": len(adjacency_matrix),
                "max_degree": max(degrees),
                "min_degree": min(degrees),
            }

    return None


def report_first_dominance_failure(n, mode, max_degree=None, min_degree=None, connected=True):
    """Print the first graph on n vertices that fails the requested property."""
    mode_labels = {
        0: "D_1 (static dominates sequential)",
        1: "D_2 (sequential dominates static)",
        2: "neither D_1 nor D_2 (crossing)",
    }

    result = find_first_dominance_failure(n, mode, max_degree, min_degree, connected)

    if result is None:
        print(f"No graph on {n} vertices failed {mode_labels[mode]}.")
        return None

    print(f"First graph failing {mode_labels[mode]}:")
    print(f"  Number of vertices: {result['num_vertices']}")
    print(f"  Max degree: {result['max_degree']}")
    print(f"  Min degree: {result['min_degree']}")
    print("  Adjacency matrix:")
    for row in result["adjacency_matrix"]:
        print(f"    {row}")

    return result


# Example usage:
N = 15
MAX_DEGREE = 3     # set to None to leave the max degree unbounded
MIN_DEGREE = None  # set to None to leave the min degree unbounded
MODE = 0           # 0 = D_1 fails, 1 = D_2 fails, 2 = neither holds

# report_first_dominance_failure(N, MODE, max_degree=MAX_DEGREE, min_degree=MIN_DEGREE)