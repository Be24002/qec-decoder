"""The repetition code: the smallest quantum error-correcting code.

d physical qubits in a line protect one logical bit; a bit flip on any
qubit is corrected by majority vote. Exact failure probability (independent
flip rate p):

    P_L(d, p) = sum_{k > (d-1)/2} C(d, k) p^k (1-p)^(d-k)

The code has threshold p_th = 50%: below it, larger d suppresses the
logical error exponentially; above it, enlarging the code makes things
worse. Same dichotomy as the surface code, but analytically solvable --
the perfect warm-up and a validation anchor for the Monte Carlo machinery.
"""

from __future__ import annotations

from math import comb

import numpy as np


def repetition_logical_error_analytic(d: int, p: float) -> float:
    """Exact majority-vote failure probability for the d-repetition code."""
    return sum(
        comb(d, k) * p**k * (1 - p) ** (d - k)
        for k in range((d + 2) // 2, d + 1)
    )


def simulate_repetition(d: int, p: float, shots: int, seed: int = 11) -> float:
    """Monte Carlo logical error rate: decode by majority vote per shot."""
    rng = np.random.default_rng(seed)
    flips = rng.random((shots, d)) < p
    wrong_majority = flips.sum(axis=1) > d / 2
    return float(wrong_majority.mean())
