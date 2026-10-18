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
import shutil
import subprocess
import os

import os
import shutil
import subprocess


import os
import shutil
import subprocess


# Set this to the full path of a working geng/nauty-geng executable if it
# cannot be found automatically (e.g. when it's only available inside a
# conda environment not on this Python process's PATH). Leave as None to
# rely on automatic detection via PATH.
GENG_PATH_OVERRIDE = r"C:\Users\Lykos\miniconda3\Library\bin\geng.exe"


def _candidate_geng_paths():
    """
    Find every occurrence of 'nauty-geng' or 'geng' across all directories
    in PATH (not just the first one), since more than one installation
    may be present (e.g. a broken MSYS2 build alongside a working conda
    build). GENG_PATH_OVERRIDE, if set, is tried first.
    """
    candidates = []
    seen = set()

    if GENG_PATH_OVERRIDE:
        candidates.append(GENG_PATH_OVERRIDE)
        seen.add(GENG_PATH_OVERRIDE)

    path_dirs = os.environ.get("PATH", "").split(os.pathsep)
    for directory in path_dirs:
        for name in ("nauty-geng.exe", "geng.exe", "nauty-geng", "geng"):
            full_path = os.path.join(directory, name)
            if os.path.isfile(full_path) and full_path not in seen:
                seen.add(full_path)
                candidates.append(full_path)

    return candidates


def _find_geng_executable():
    """
    Locate a working nauty graph generator executable and its containing
    directory. Tries every candidate (GENG_PATH_OVERRIDE first, then PATH)
    with a trivial call (n=1) and returns the first one that actually runs
    successfully, skipping broken installations.

    Returns (executable_path, working_directory).
    """
    for candidate in _candidate_geng_paths():
        candidate_dir = os.path.dirname(candidate)
        try:
            subprocess.run(
                [candidate, "-q", "1"],
                capture_output=True, text=True, check=True, cwd=candidate_dir,
            )
            return candidate, candidate_dir
        except (subprocess.CalledProcessError, OSError):
            continue

    raise FileNotFoundError(
        "Could not find a working 'nauty-geng' or 'geng' executable. "
        "Set GENG_PATH_OVERRIDE in GeneralSearch.py to the full path of a "
        "working executable (e.g. from 'where geng' inside Anaconda Prompt), "
        "or make sure that directory is on this Python process's PATH."
    )


def generate_symmetric_matrices(n, max_degree=None, min_degree=None, connected=True):
    """
    Generate all non-isomorphic simple graphs on n vertices (as adjacency
    matrices) using nauty-geng, optionally bounded by max/min degree.
    """
    if n < 1:
        raise ValueError("n must be at least 1")

    # geng has a known edge case with n=1 combined with -c. Handle n=1
    # directly: the only simple graph on 1 vertex is an isolated vertex.
    if n == 1:
        degree = 0
        if max_degree is not None and degree > max_degree:
            return
        if min_degree is not None and degree < min_degree:
            return
        yield [[0]]
        return

    geng_path, geng_dir = _find_geng_executable()

    command = [geng_path, '-q']
    if connected:
        command.append('-c')
    if max_degree is not None:
        command.append(f'-D{max_degree}')
    if min_degree is not None:
        command.append(f'-d{min_degree}')
    command.append(str(n))

    try:
        result = subprocess.run(
            command, capture_output=True, text=True, check=True, cwd=geng_dir
        )
    except subprocess.CalledProcessError as error:
        raise RuntimeError(
            f"geng failed for command {command}.\n"
            f"stderr: {error.stderr.strip()}"
        ) from error

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
        0: "D_1 (sequential dominates static)",
        1: "D_2 (static dominates sequential)",
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