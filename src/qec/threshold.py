"""Monte Carlo threshold measurements for the surface code."""

from __future__ import annotations

import numpy as np

from .lattice import SurfaceCode
from .matching import LatticeContext, decode


def simulate_surface(d: int, p: float, shots: int, seed: int = 7) -> tuple[float, int]:
    """Logical X-failure rate of the rotated surface code at physical
    X-error rate p (perfect stabilizer measurement, X sector; the Z sector
    is symmetric for this noise model).

    Returns (logical_error_rate, failures).
    """
    code = SurfaceCode(d)
    ctx = LatticeContext(code.z_stabs)
    qubits = sorted(code.qubits)
    rng = np.random.default_rng(seed)

    failures = 0
    for _ in range(shots):
        errors = {q for q, r in zip(qubits, rng.random(len(qubits)) < p) if r}
        defects = code.syndrome(errors, "X")
        correction = decode([stab.center for stab in defects], ctx, d)
        if code.logical_x_failure(errors ^ correction):
            failures += 1
    return failures / shots, failures


def crossing_point(ps, rates_a, rates_b):
    """Linear-interpolation crossing of two error-rate curves."""
    ps = np.asarray(ps)
    ra, rb = np.asarray(rates_a), np.asarray(rates_b)
    diff = ra - rb
    crossings = []
    for i in range(len(ps) - 1):
        if diff[i] == 0:
            crossings.append(ps[i])
        elif diff[i] * diff[i + 1] < 0:
            t = diff[i] / (diff[i] - diff[i + 1])
            crossings.append(ps[i] + t * (ps[i + 1] - ps[i]))
    return float(np.mean(crossings)) if crossings else float("nan")
