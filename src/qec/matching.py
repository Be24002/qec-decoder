"""Minimum-weight perfect matching decoder for the surface code.

Syndrome decoding maps to matching: each defect (flipped stabilizer) is the
endpoint of an error chain; chains pair up (or run to a boundary), and the
most likely pairing minimizes total chain length (the Manhattan metric for
uniform error rates). MWPM is solved with networkx's Blossom algorithm.

Corrections are built **with syndrome verification**: a candidate chain is
only accepted if its syndrome is exactly the target defect set. Pure
geometry (staircase paths) can silently produce invalid Pauli operators on
this staggered-boundary lattice, so every chain is verified; long chains
are decomposed into verified hops between nearby Z-stabilizers, and all
chains are memoized (Monte Carlo hits the same pairs constantly).
"""

from __future__ import annotations

import itertools

import networkx as nx


def _manhattan(a, b) -> int:
    return int(abs(a[0] - b[0]) + abs(a[1] - b[1]))


class LatticeContext:
    """Decoder-side view of one code instance (X-error sector)."""

    def __init__(self, z_stabs):
        self.by_center = {stab.center: stab for stab in z_stabs}
        self.members = {}  # qubit -> centers of Z-stabs containing it
        for center, stab in self.by_center.items():
            for qubit in stab.qubits:
                self.members.setdefault(qubit, set()).add(center)
        centers = sorted(self.by_center)
        self.neighbors = {
            c: [c2 for c2 in centers if c != c2 and _manhattan(c, c2) <= 2]
            for c in centers
        }
        # BFS all-pairs distances on the Z-stabilizer sublattice: the true
        # chain cost between two defects is the number of lattice hops, not
        # the Manhattan distance (diagonal Z-neighbors connect with short
        # "kite" chains that Manhattan overcounts by 2x).
        self.dist = {}
        for source in centers:
            dist = {source: 0}
            queue = [source]
            while queue:
                nxt = []
                for node in queue:
                    for adj in self.neighbors[node]:
                        if adj not in dist:
                            dist[adj] = dist[node] + 1
                            nxt.append(adj)
                queue = nxt
            self.dist[source] = dist
        self._pair_cache: dict = {}
        self._boundary_cache: dict = {}

    # ------------------------------------------------------------ syndrome

    def syndrome(self, chain) -> frozenset:
        """Z-stabilizer centers flipped by an X-chain on these qubits."""
        counts: dict = {}
        for qubit in chain:
            for center in self.members.get(qubit, ()):  # noqa: B007
                counts[center] = counts.get(center, 0) + 1
        return frozenset(c for c, n in counts.items() if n % 2 == 1)

    # ------------------------------------------------------ verified chains

    def _hop_chain(self, a, b) -> frozenset:
        """Verified chain for a NEAR pair (Manhattan <= 2): search small
        qubit subsets of the two plaquettes until the syndrome is {a, b}."""
        target = frozenset({a, b})
        pool = sorted(set(self.by_center[a].qubits) | set(self.by_center[b].qubits))
        for size in (1, 2, 3, 4):
            for combo in itertools.combinations(pool, size):
                if self.syndrome(frozenset(combo)) == target:
                    return frozenset(combo)
        raise ValueError(f"no verified hop chain between {a} and {b}")

    def _hop_cost(self, a, b) -> int:
        return len(self._hop_chain(a, b))

    def _path_on_z_lattice(self, a, b) -> list:
        """Cheapest route a -> b on the Z-sublattice: Dijkstra with the
        REAL hop-chain lengths as edge costs (hop counts mislead: two-hop
        routes can differ by 2x in chain qubits)."""
        import heapq

        cost = {a: 0}
        parent = {a: None}
        heap = [(0, a)]
        while heap:
            total, node = heapq.heappop(heap)
            if node == b:
                break
            if total > cost[node]:
                continue
            for adj in self.neighbors[node]:
                new = total + self._hop_cost(node, adj)
                if adj not in cost or new < cost[adj]:
                    cost[adj] = new
                    parent[adj] = node
                    heapq.heappush(heap, (new, adj))
        if b not in parent:
            raise ValueError(f"no Z-lattice path {a} -> {b}")
        path = [b]
        while parent[path[-1]] is not None:
            path.append(parent[path[-1]])
        return list(reversed(path))

    def graph_distance(self, a, b) -> int:
        """BFS hop count between two Z-centers (~ chain qubits / 2)."""
        return self.dist[a][b]

    def pair_chain(self, a, b) -> frozenset:
        """Minimal verified chain with syndrome exactly {a, b}."""
        key = (a, b) if a <= b else (b, a)
        if key not in self._pair_cache:
            try:  # direct short chain first (often optimal)
                chain = self._hop_chain(a, b)
            except ValueError:  # noqa: TRY003 -- route through the sublattice
                chain = set()
                path = self._path_on_z_lattice(a, b)
                for u, v in zip(path, path[1:]):
                    chain ^= self._hop_chain(u, v)
                chain = frozenset(chain)
            self._pair_cache[key] = chain
        return self._pair_cache[key]

    def _search_boundary_chain(self, a, d: int, side: str) -> frozenset:
        """Straight vertical walk to the top/bottom boundary with a
        left/right choice at each crossing; keep the first verified one."""
        target = frozenset({a})
        cx, cy = a
        if side == "T":
            steps = int(cy + 0.5)
            rows = [int(cy - 0.5 - k) for k in range(steps)]
        else:
            steps = int(d - 0.5 - cy)
            rows = [int(cy + 0.5 + k) for k in range(steps)]
        crossings = []
        for y in rows:
            options = [
                q for q in ((int(cx - 0.5), y), (int(cx + 0.5), y))
                if 0 <= q[0] < d and 0 <= q[1] < d
            ]
            if not options:
                return frozenset()
            crossings.append(options)
        for choices in itertools.product(*crossings):
            chain = frozenset(choices)
            if self.syndrome(chain) == target:
                return chain
        return frozenset()

    def boundary_chain(self, a, d: int, side: str) -> frozenset:
        """Verified chain from defect a to the given boundary (T or B)."""
        if (a, side) not in self._boundary_cache:
            self._boundary_cache[(a, side)] = self._search_boundary_chain(a, d, side)
        chain = self._boundary_cache[(a, side)]
        if not chain:
            raise ValueError(f"no verified boundary chain from {a} to {side}")
        return chain

    def nearest_boundary_chain(self, a, d: int) -> frozenset:
        """Chain to the nearer boundary (used for single-defect syndromes)."""
        _, cy = a
        nearer = "T" if cy + 0.5 <= d - 0.5 - cy else "B"
        for side in (nearer, "T" if nearer == "B" else "B"):
            if (a, side) not in self._boundary_cache:
                self._boundary_cache[(a, side)] = self._search_boundary_chain(a, d, side)
            chain = self._boundary_cache[(a, side)]
            if chain:
                return chain
        raise ValueError(f"no verified boundary chain from {a}")


def _side_distances(center, d: int) -> tuple[int, int]:
    """Chain qubit counts from a defect to the top and bottom boundaries."""
    _, cy = center
    top = int(round(cy + 0.5))
    bottom = int(round(d - 0.5 - cy))
    return top, bottom


def decode(defects, ctx: LatticeContext, d: int) -> set:
    """Correction (set of data qubits) for the given Z-stabilizer defects.

    Matching graph (rigorous boundary handling):
    - defect_i -- defect_j with Manhattan weights;
    - defect_i -- T_i and defect_i -- B_i with the distance to that side
      (side-aware: a chain must exit the boundary the decoder chose);
    - T_i -- T_j and B_i -- B_j at weight 0 (any number of chains may exit
      the same side); T-B edges are forbidden (top and bottom must not mix,
      that would be a free logical flip);
    - one dummy node (weight-0 to all T/B) fixes the parity when needed.
    """
    defects = list(defects)
    if not defects:
        return set()
    if len(defects) == 1:
        return set(ctx.nearest_boundary_chain(defects[0], d))

    n = len(defects)
    graph = nx.Graph()
    for i in range(n):
        for j in range(i + 1, n):
            # exact chain length = -log-likelihood (not a geometric proxy)
            graph.add_edge(i, j, weight=len(ctx.pair_chain(defects[i], defects[j])))
        for side in ("T", "B"):
            try:
                weight = len(ctx.boundary_chain(defects[i], d, side))
            except ValueError:
                top, bottom = _side_distances(defects[i], d)
                weight = top if side == "T" else bottom
            graph.add_edge(i, (side, i), weight=weight)
    for i in range(n):
        for j in range(n):
            graph.add_edge(("T", i), ("T", j), weight=0)
            graph.add_edge(("B", i), ("B", j), weight=0)
            graph.add_edge(("T", i), ("B", j), weight=0)
    if n % 2 == 1:  # keep the total node count even
        for i in range(n):
            graph.add_edge("dummy", ("T", i), weight=0)
            graph.add_edge("dummy", ("B", i), weight=0)

    pairs = nx.algorithms.matching.min_weight_matching(graph, weight="weight")

    correction: set = set()
    for u, v in pairs:
        u_is_defect = isinstance(u, int)
        v_is_defect = isinstance(v, int)
        if u_is_defect and v_is_defect:
            correction ^= ctx.pair_chain(defects[u], defects[v])
        elif u_is_defect or v_is_defect:
            i, virtual = (u, v) if u_is_defect else (v, u)
            correction ^= ctx.boundary_chain(defects[i], d, virtual[0])
        # else: virtual-virtual pairs carry no chain
    return correction
