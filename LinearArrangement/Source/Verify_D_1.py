import math
import numpy as np

from Calc_Probabilities import (
    possible_indices,
    N_stat_k,
    probability_matrix_list,
)

# Numerical tolerance: N_stat_k comes from symbolic factorials (sympy)
# evaluated to float, while probability_matrix comes from a recursive
# sum in plain Python. These are two different numerical paths, so raw
# comparisons can show spurious differences of order 1e-9 to 1e-14 that
# are floating-point noise, not real violations of stochastic dominance.
TOLERANCE = 1e-9


def is_dominant(a, b, tol=TOLERANCE):
    """True if a >= b, allowing for floating-point noise of size tol."""
    return a >= b - tol


max_n = len(probability_matrix_list)

# D_1: for each fixed n, checks whether the static model's cumulative
# distribution dominates the sequential model's cumulative distribution.
# Stores full per-k data for later plotting; prints only the pass/fail
# summary per n.
accumulated_data = {}

for n in range(1, max_n + 1):
    ks = list(possible_indices(n))
    stat_probs = [N_stat_k(n, k) for k in ks]

    seq_probs = probability_matrix_list[n - 1]
    seq_probs = seq_probs[ks[0]:ks[0] + len(ks)]

    cumulative_stat = np.cumsum(stat_probs)
    cumulative_seq = np.cumsum(seq_probs)

    accumulated_data[n] = {
        "ks": ks,
        "cumulative_stat": cumulative_stat,
        "cumulative_seq": cumulative_seq,
    }

    d1_holds = all(
        is_dominant(cumulative_stat[i], cumulative_seq[i])
        for i in range(len(cumulative_seq))
    )

    if d1_holds:
        print(f"n = {n}: D_1 holds")
    else:
        print(f"n = {n}: D_1 fails")