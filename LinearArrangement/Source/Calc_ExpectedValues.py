from Calc_Probabilities import (
    padovan_sequence,
    probability_matrix,
    possible_indices,
    N_stat_k,
    N_seq_k,
)
import sympy as sp
import numpy as np
import math



expected_seats_sequential = []


for n_seats in range(0, len(probability_matrix), 5):
    print(f"Expected occupied seats for n_seats = {n_seats}")

    expected_value = 0
    for k in range(len(probability_matrix[n_seats])):
        if probability_matrix[n_seats][k] != 0:
            expected_value += k * probability_matrix[n_seats][k]

    if n_seats == 0:
        expected_value = 0
    else:
        expected_value = expected_value / n_seats

    print(f"expected_value = {expected_value}")
    expected_seats_sequential.append(expected_value)
    
    
    expected_seats_static = []

for n_seats in range(0, 200, 5):
    print(f"Expected occupied seats for n_seats = {n_seats}")

    if n_seats == 0:
        expected_value = 0
    else:
        ks = list(possible_indices(n_seats))
        stat_probs = [N_stat_k(n_seats, k) for k in ks]

        expected_value = sum(k * p for k, p in zip(ks, stat_probs))
        expected_value = expected_value / n_seats

    print(f"expected_value = {expected_value}")
    expected_seats_static.append(expected_value)