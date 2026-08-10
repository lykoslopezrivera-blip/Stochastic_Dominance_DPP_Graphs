"""
In this first code we will give functions to calculate probabilities of the random variables defined.
We will use the following notation:
X_stat(n,k) := Indicator random variable of the k-th position in the linear arrangement with n sites, in the static model.
X_sec(n,k) := Indicator random variable of the k-th position in the linear arrangement with n sites, in the sequential model.
N_sec(n,k) := Number of elements in the k-th position in the sequential model.
N_stat(n,k) := Number of elements in the k-th position in the static model.

"""
from math import factorial
import sympy as sp
import numpy as np
import math

def generate_padovan_sequence(n_terms: int) -> list[int]:
    """
    Generate the first n_terms of the Padovan sequence.
    Recurrence: P(n) = P(n-2) + P(n-3), with P(0)=P(1)=P(2)=1.
    """
    if n_terms <= 3:
        return [1] * n_terms

    padovan_sequence = [1, 1, 1]
    for i in range(3, n_terms):
        next_term = padovan_sequence[i - 2] + padovan_sequence[i - 3]
        padovan_sequence.append(next_term)

    return padovan_sequence


padovan_sequence = generate_padovan_sequence(300)



def X_stat(n: int, i: int) -> float:
    """
    Probability that seat i is occupied in a line of n seats,
    computed directly from the Padovan sequence.
    X_stat(n, i) = P(i-1) * P(n-i) / P(n+1)
    """
    return padovan_sequence[i - 1] * padovan_sequence[n - i] / padovan_sequence[n + 1]


def d(i: int, l: int) -> float:
    """
    Auxiliary term d_{i, i+l} used in the formula of X_seq.
    """
    total = 0
    for k in range(i):
        total += (((-1) ** (l + k)) / factorial(k)) * \
                  ((1 / factorial(i + l - k)) - (1 / factorial(l)))
    return total


def X_seq(n: int, i: int) -> float:
    """
    Probability that seat i is occupied in a line of n seats,
    computed as an alternating series expansion.
    """
    if (n - i) % 2 == 0:
        return 1 + sum(d(i, 2 * k) for k in range((n - i) // 2 + 1))
    else:
        return sum(d(i, 2 * k + 1) for k in range((n - i - 1) // 2 + 1))
    
    
    



x, y = sp.symbols('x y')

# Bivariate generating function T(x, y) for the occupied-seat probability sequence.
tanh_term = sp.tanh(y * sp.sqrt(x))
numerator = 1 + sp.sqrt(x) * tanh_term
denominator = 1 - y + (1 / sp.sqrt(x) - y * sp.sqrt(x)) * tanh_term
generating_function = numerator / denominator

# Expand generating_function as a power series around x=0, y=0.
# Order 16 in x and order 31 in y. This step is expensive (only run once);
# the resulting expansion is stored below in bivariate_series so it does not
# need to be recomputed.
#
# bivariate_series = sp.series(generating_function, x, 0, 16).removeO()
# bivariate_series = sp.series(bivariate_series, y, 0, 31).removeO()

# Precomputed result of the expansion above, used as the base case for the
# recursive computation of subsequent terms.
bivariate_series = x**2*y**4 + x*y**2 + x*y + y**30*(5628330861013796436068*x**15/136295711075535652265625 + 46349589728906493860692*x**14/136295711075535652265625 + 1299304445882170352*x**13/2762657567153859375 + 1740076870828*x**12/12342163640625 + 1135508*x**11/155003625 + x**10/59049) + y**29*(222148423805582000341*x**15/37643577344671751578125 + 344374976144793844256*x**14/1981240912877460609375 + 6199623310388210512*x**13/12362422773291215625 + 132442781175776*x**12/461596920159375 + 17927528*x**11/568346625 + 32*x**10/98415) + y**28*(18788936697793261621*x**14/336532619939594203125 + 133109442203412905999*x**13/336532619939594203125 + 11172459111055223*x**12/25011714599746875 + 4741197443*x**11/47888465625 + 3379*x**10/1148175) + y**27*(848070074996694008*x**14/99850337784275203125 + 21771302215301917807*x**13/99850337784275203125 + 151949719144261946*x**12/289421268939928125 + 6296869550203*x**11/27152760009375 + 24356*x**10/1476225 + x**9/19683) + y**26*(144554065878772801*x**13/1923043542511966875 + 870237226116499562*x**12/1923043542511966875 + 1394883397268*x**11/3419236445625 + 9351311*x**10/147349125 + 29*x**9/32805) + y**25*(886395722771204*x**13/72512954091703125 + 333744863127793772*x**12/1232720219558953125 + 212148290168599*x**11/397011342853125 + 3601239736*x**10/20523628125 + 394*x**9/54675) + y**24*(49832352543044*x**12/494737211863125 + 653278670528537*x**11/1286316750844125 + 71445168034*x**10/201131555625 + 13793*x**9/382725 + x**8/6561) + y**23*(78863079581996*x**12/4482618980214375 + 135383714754929*x**11/407510816383125 + 279382408*x**10/531349875 + 53934478*x**9/442047375 + 26*x**8/10935) + y**22*(13703146165504*x**11/102088631019375 + 11374671304873*x**10/20417726203875 + 2791221934*x**9/9577693125 + 88*x**8/5103) + y**21*(197394815194*x**11/7795859096025 + 15671100026704*x**10/38979295480125 + 2218701916*x**9/4469590125 + 87016*x**8/1148175 + x**7/2187) + y**20*(86814383846*x**10/488462349375 + 58077624047*x**9/97692469875 + 163163863*x**8/736745625 + 23*x**7/3645) + y**19*(7515775348*x**10/206239658625 + 98690223518*x**9/206239658625 + 8609126*x**8/19348875 + 341*x**7/8505) + y**18*(1342128374*x**9/5746615875 + 28192856*x**8/45972927 + 41506*x**7/273375 + x**6/729) + y**17*(569311642*x**9/10854718875 + 6059228432*x**8/10854718875 + 6104827*x**7/16372125 + 4*x**6/243) + y**16*(2589134*x**8/8513505 + 25799036*x**7/42567525 + 1273*x**6/14175) + y**15*(6885458*x**8/91216125 + 57933868*x**7/91216125 + 60656*x**6/212625 + x**5/243) + y**14*(1282348*x**7/3274425 + 1854632*x**6/3274425 + 17*x**5/405) + y**13*(220217*x**7/2027025 + 1418992*x**6/2027025 + 904*x**5/4725) + y**12*(21149*x**6/42525 + 20851*x**5/42525 + x**4/81) + y**11*(4876*x**6/31185 + 4615*x**5/6237 + 14*x**4/135) + y**10*(977*x**5/1575 + 598*x**4/1575) + y**9*(638*x**5/2835 + 2092*x**4/2835 + x**3/27) + y**8*(34*x**4/45 + 11*x**3/45) + y**7*(34*x**4/105 + 71*x**3/105) + y**6*(8*x**3/9 + x**2/9) + y**5*(7*x**3/15 + 8*x**2/15) + y**3*(2*x**2/3 + x/3) + 1





def possible_indices(L):
    """
    Range of valid k values (occupied seats) for a line of length L,
    based on the combinatorial bounds of the Padovan-based process.
    """
    return range(math.floor((L + 2) / 3), math.floor((L + 1) / 2) + 1)


def N_stat_k(n, k):
    """
    Probability of exactly k occupied seats in a line of n seats,
    computed directly with the closed formula (Padovan-based).
    """
    value = sp.factorial(k + 1) / (
        sp.factorial(3 * k - n) * sp.factorial(1 + n - 2 * k) * padovan_sequence[n + 1]
    )
    return value.evalf()


def N_seq_k(n, k):
    """
    Probability of exactly k occupied seats in a line of n seats,
    extracted as the coefficient of x^k * y^n in the bivariate
    generating function series expansion.
    """
    coeff = sp.Poly(bivariate_series, x, y).coeff_monomial(x**k * y**n)
    return coeff.evalf()



probability_matrix = []

# Base cases: n = 1, 2, 3, 4, computed directly from the series expansion.
for n in range(1, 5):
    possible_ks = list(possible_indices(n))
    seq_probs = [N_seq_k(n, k) for k in possible_ks]

    # Pad the start of the row with zeros so index k matches position k.
    leading_zeros = [0] * possible_ks[0]
    probability_matrix.append(leading_zeros + seq_probs)

max_seats = 200

# Pad every row on the right with zeros up to max_seats, so all rows
# have equal length and can later be stacked into a single array.
for row_index, row in enumerate(probability_matrix):
    remaining = max_seats - len(row)
    probability_matrix[row_index] = row + [0] * remaining

# Recursive computation for L = 5 up to max_seats.
# Each line of length L is split by an "unfriendly" seat at position k
# into a left sub-line (length k-1) and a right sub-line (length L-k).
# i and j are the occupied-seat counts in the left/right sub-lines,
# with i + j = n - 1 (n = total occupied seats in the line of length L).
for L in range(5, max_seats + 1):
    probability_vector = []
    possible_ns = list(possible_indices(L))

    for n in possible_ns:
        outer_sum = 0
        for k in range(1, L + 1):
            inner_sum = 0

            # j is fully determined by i (since i + j = n - 1).
            for i in range(0, n):
                j = n - 1 - i

                if k == 1 or k == 2:
                    # Left sub-line is empty or trivial: only the right
                    # sub-line contributes, and only when i = 0.
                    if i == 0 and math.floor((L - k) / 3) - 1 <= j <= math.floor((L - k - 1) / 2) + 1:
                        inner_sum += probability_matrix[L - k - 2][j]

                elif k == L or k == L - 1:
                    # Right sub-line is empty or trivial: only the left
                    # sub-line contributes, and only when j = 0.
                    if j == 0 and math.floor((k - 1) / 3) - 1 <= i <= math.floor((k - 2) / 2) + 1:
                        inner_sum += probability_matrix[k - 3][i]

                else:
                    # General case: both sub-lines contribute, independently.
                    left_valid = math.floor((k - 1) / 3) - 1 <= i <= math.floor((k - 2) / 2) + 1
                    right_valid = math.floor((L - k) / 3) - 1 <= j <= math.floor((L - k - 1) / 2) + 1
                    if left_valid and right_valid:
                        inner_sum += probability_matrix[k - 3][i] * probability_matrix[L - k - 2][j]

            outer_sum += inner_sum

        probability_vector.append(outer_sum / L)

    # Pad the new row (left with zeros, right up to max_seats) same as before.
    leading_zeros = [0] * possible_ns[0]
    probability_vector = leading_zeros + probability_vector
    remaining = max_seats - len(probability_vector)
    probability_vector += [0] * remaining

    probability_matrix.append(probability_vector)

probability_matrix_list = probability_matrix
probability_matrix = np.array(probability_matrix_list, dtype=float)