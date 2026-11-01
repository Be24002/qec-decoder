"""Physics regression tests: lattice invariants, decoder validity,
analytic repetition anchors, and sub-threshold scaling."""

import numpy as np
import pytest

from qec.lattice import SurfaceCode
from qec.matching import LatticeContext, decode
from qec.repetition import repetition_logical_error_analytic, simulate_repetition
from qec.threshold import simulate_surface


def make(d=3):
    code = SurfaceCode(d)
    return code, LatticeContext(code.z_stabs)


# --------------------------------------------------------------- lattice


def test_stabilizer_counts():
    for d, total in [(3, 8), (5, 24), (7, 48)]:
        code = SurfaceCode(d)
        assert len(code.stabilizers) == total == d * d - 1
        assert len(code.z_stabs) == len(code.x_stabs) == (d * d - 1) // 2


def test_stabilizer_weights_are_2_or_4():
    code = SurfaceCode(5)
    assert all(len(s.qubits) in (2, 4) for s in code.stabilizers)


def test_every_qubit_is_covered():
    code = SurfaceCode(3)
    covered = set().union(*(s.qubits for s in code.stabilizers))
    assert covered == code.qubits


def test_logical_operators_are_weight_d():
    code = SurfaceCode(5)
    assert len(code.logical_z_row(0)) == 5
    assert len(code.logical_x_column(0)) == 5


def test_column_is_logical_row_is_detectable():
    # A full column of X commutes with every Z stabilizer (empty syndrome);
    # a full row of X does not (it is detectable, hence not the logical).
    code = SurfaceCode(3)
    column = set(code.logical_x_column(0))
    row = set(code.logical_z_row(0))
    assert code.syndrome(column, "X") == frozenset()
    assert len(code.syndrome(row, "X")) > 0


def test_single_error_syndrome_structure():
    code = SurfaceCode(3)
    # interior qubit -> two Z-type detectors fire
    assert len(code.syndrome({(1, 1)}, "X")) == 2
    # corner qubit -> one detector fires
    assert len(code.syndrome({(0, 0)}, "X")) == 1


# --------------------------------------------------------------- decoder


def test_zero_defects_give_empty_correction():
    _, ctx = make(3)
    assert decode([], ctx, 3) == set()


def test_correction_cancels_the_syndrome():
    # Core decoder-validity property: for any error pattern, the returned
    # correction must reproduce (and thus cancel) the observed syndrome.
    rng = np.random.default_rng(4)
    code, ctx = make(3)
    qubits = sorted(code.qubits)
    for _ in range(200):
        errors = {q for q, r in zip(qubits, rng.random(9) < 0.15) if r}
        defects = code.syndrome(errors, "X")
        correction = decode([s.center for s in defects], ctx, 3)
        residual = code.syndrome(errors ^ correction, "X")
        assert residual == frozenset()


def test_single_errors_are_always_corrected():
    code, ctx = make(3)
    for qubit in code.qubits:
        correction = decode([s.center for s in code.syndrome({qubit}, "X")], ctx, 3)
        assert not code.logical_x_failure({qubit} ^ correction)


def test_full_logical_column_causes_failure():
    code = SurfaceCode(3)
    errors = set(code.logical_x_column(1))
    assert code.syndrome(errors, "X") == frozenset()  # undetectable by design
    assert code.logical_x_failure(errors)  # and it flips the logical Z


def test_two_defects_pair_to_each_other():
    # Two nearby defects: matching should connect them, correcting both
    # with a short chain and no logical damage.
    code, ctx = make(5)
    errors = {(1, 1), (2, 1)}  # adjacent X errors inside one plaquette pair
    defects = code.syndrome(errors, "X")
    correction = decode([s.center for s in defects], ctx, 5)
    assert not code.logical_x_failure(errors ^ correction)


# ------------------------------------------------------------ repetition


@pytest.mark.parametrize("d,p", [(5, 0.10), (7, 0.30)])
def test_repetition_monte_carlo_matches_analytic(d, p):
    shots = 40_000
    measured = simulate_repetition(d, p, shots)
    exact = repetition_logical_error_analytic(d, p)
    sigma = np.sqrt(max(measured * (1 - measured), 1e-12) / shots)
    assert measured == pytest.approx(exact, abs=4 * sigma)


def test_repetition_improves_with_d_below_half():
    p = 0.30  # below the 50% threshold: bigger code = fewer logical errors
    assert repetition_logical_error_analytic(9, p) < repetition_logical_error_analytic(3, p)


def test_repetition_worsens_with_d_above_half():
    p = 0.60  # above threshold: bigger code = more logical errors
    assert repetition_logical_error_analytic(9, p) > repetition_logical_error_analytic(3, p)


# ------------------------------------------------- surface code scaling


def test_surface_code_subthreshold_scaling():
    # Well below threshold (~10%), distance 5 must beat distance 3.
    p = 0.03
    rate_3, _ = simulate_surface(3, p, shots=3000, seed=1)
    rate_5, _ = simulate_surface(5, p, shots=3000, seed=2)
    assert rate_5 < rate_3


def test_surface_code_above_threshold_scaling():
    # Well above threshold, growing the code makes things worse.
    p = 0.30
    rate_3, _ = simulate_surface(3, p, shots=1500, seed=3)
    rate_7, _ = simulate_surface(7, p, shots=1500, seed=4)
    assert rate_7 > rate_3
