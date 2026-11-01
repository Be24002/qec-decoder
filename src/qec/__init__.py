"""qec: quantum error correction simulated and measured, not asserted.

Modules
-------
lattice
    The rotated surface code: data qubits, X/Z stabilizers, syndromes,
    logical operators and logical-failure tests.
matching
    Minimum-weight perfect matching decoder (networkx blossom) with
    boundary handling, plus the lattice-path corrections.
repetition
    The 1D repetition code with majority-vote decoding and its exact
    analytic failure probability.
threshold
    Monte Carlo sweeps: logical error rate vs physical error rate and code
    distance, and the threshold (crossing point) estimation.
"""

__version__ = "0.1.0"
