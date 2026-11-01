"""The rotated surface code, built from scratch on a square lattice.

Layout (distance d): d^2 data qubits at integer coordinates (x, y) with
0 <= x, y < d. Stabilizers live at half-integer plaquette centers:

- interior (x+1/2, y+1/2), 0 <= x, y <= d-2: Z-type if (x+y) is even, X-type
  otherwise (a checkerboard of weight-4 plaquettes);
- top edge (x+1/2, -1/2) with x even, and bottom edge (x+1/2, d-1/2) with
  x odd: weight-2 X-type;
- left edge (-1/2, y+1/2) with y odd, and right edge (d-1/2, y+1/2) with
  y even: weight-2 Z-type.

The staggered edge pattern is what makes the boundaries "rough/smooth" so
that the code encodes exactly one logical qubit:

- logical Z = a full ROW of Z operators (joins the left/right edges),
- logical X = a full COLUMN of X operators (joins the top/bottom edges),

both of minimum weight d. Consequences used everywhere below:

- an X-error pattern fails correction iff it anticommutes with logical Z,
  i.e. iff the number of X errors in any single row is odd;
- a Z error on a data qubit anticommutes with (flips) the X-type
  stabilizers that contain it, and vice versa.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Stabilizer:
    """One stabilizer: its Pauli type and the data qubits it covers."""

    kind: str                 # "X" or "Z"
    center: tuple[float, float]
    qubits: frozenset[tuple[int, int]]


class SurfaceCode:
    """Rotated planar surface code of distance ``d``."""

    def __init__(self, d: int):
        if d < 3 or d % 2 == 0:
            raise ValueError("rotated surface code needs odd d >= 3")
        self.d = d
        self.qubits = frozenset((x, y) for x in range(d) for y in range(d))
        self.stabilizers = self._build_stabilizers()
        self.z_stabs = [s for s in self.stabilizers if s.kind == "Z"]
        self.x_stabs = [s for s in self.stabilizers if s.kind == "X"]

    # ----------------------------------------------------------- geometry

    def _plaquette_qubits(self, cx: float, cy: float) -> frozenset:
        candidates = []
        for qx in (cx - 0.5, cx + 0.5):
            for qy in (cy - 0.5, cy + 0.5):
                candidates.append((int(qx), int(qy)))
        return frozenset(q for q in candidates if q in self.qubits)

    def _build_stabilizers(self) -> list[Stabilizer]:
        d = self.d
        stabs: list[Stabilizer] = []
        # interior checkerboard
        for x in range(d - 1):
            for y in range(d - 1):
                kind = "Z" if (x + y) % 2 == 0 else "X"
                center = (x + 0.5, y + 0.5)
                stabs.append(Stabilizer(kind, center, self._plaquette_qubits(*center)))
        # staggered weight-2 boundaries
        for x in range(0, d - 1, 2):  # top: x even
            center = (x + 0.5, -0.5)
            stabs.append(Stabilizer("X", center, self._plaquette_qubits(*center)))
        for x in range(1, d - 1, 2):  # bottom: x odd
            center = (x + 0.5, d - 0.5)
            stabs.append(Stabilizer("X", center, self._plaquette_qubits(*center)))
        for y in range(1, d - 1, 2):  # left: y odd
            center = (-0.5, y + 0.5)
            stabs.append(Stabilizer("Z", center, self._plaquette_qubits(*center)))
        for y in range(0, d - 1, 2):  # right: y even
            center = (self.d - 0.5, y + 0.5)
            stabs.append(Stabilizer("Z", center, self._plaquette_qubits(*center)))
        return stabs

    # ---------------------------------------------------------- syndrome

    def syndrome(self, errors: set, kind: str) -> frozenset:
        """Stabilizers of the given type flipped by the error pattern.

        An X error flips the Z-type stabilizers containing it (odd overlap).
        """
        detectors = self.z_stabs if kind == "X" else self.x_stabs
        return frozenset(
            stab for stab in detectors if len(stab.qubits & errors) % 2 == 1
        )

    # ------------------------------------------------------ logical ops

    def logical_z_row(self, row: int = 0) -> frozenset:
        return frozenset((x, row) for x in range(self.d))

    def logical_x_column(self, column: int = 0) -> frozenset:
        return frozenset((column, y) for y in range(self.d))

    def logical_x_failure(self, x_errors: set) -> bool:
        """True iff the X-error pattern anticommutes with logical Z.

        Equivalently: odd number of X errors in any row (all rows agree
        because the pattern commutes with every Z stabilizer after
        correction; we test on row 0 directly on error+correction).
        """
        return len(x_errors & self.logical_z_row(0)) % 2 == 1
