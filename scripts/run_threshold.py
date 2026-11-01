"""The surface code threshold: logical error rate vs physical error rate.

Sweeps the rotated surface code (perfect stabilizer measurements, X-error
sector) over code distances; the curves cross near the famous ~10.3%
matching threshold -- below it, larger codes suppress logical errors
exponentially; above it, they amplify them.

Produces figures/surface_threshold.png and results/surface_threshold.csv.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from qec.threshold import crossing_point, simulate_surface

ROOT = Path(__file__).resolve().parents[1]
FIGURES = ROOT / "figures"
RESULTS = ROOT / "results"
FIGURES.mkdir(exist_ok=True)
RESULTS.mkdir(exist_ok=True)

p_grid = [0.04, 0.06, 0.08, 0.095, 0.105, 0.12, 0.15, 0.20]
shots = {3: 4000, 5: 2500, 7: 1500}

curves = {}
for d in (3, 5, 7):
    rates = []
    for p in p_grid:
        rate, _ = simulate_surface(d, p, shots[d], seed=10 * d + int(1000 * p))
        rates.append(rate)
        print(f"d = {d}, p = {p:.3f}: logical error = {rate:.4f}")
    curves[d] = rates

threshold_est = crossing_point(p_grid, curves[3], curves[7])
print(f"threshold estimate (d=3 x d=7 crossing): {threshold_est:.4f}")
print("expected for MWPM on the code-capacity surface code: ~0.103")

fig, ax = plt.subplots(figsize=(7.5, 4.8))
for d, rates in curves.items():
    ax.plot(p_grid, rates, "o-", label=f"d = {d}")
ax.axvline(0.103, color="C3", ls=":", lw=1.2, label="known threshold ~10.3%")
ax.set_yscale("log")
ax.set_xlabel("physical error rate $p$")
ax.set_ylabel("logical error rate $p_L$ per shot")
ax.set_title(
    f"Rotated surface code threshold (MWPM): crossing at p ~ {threshold_est:.3f}"
)
ax.legend()
ax.grid(alpha=0.3, which="both")
fig.tight_layout()
fig.savefig(FIGURES / "surface_threshold.png", dpi=150)
plt.close(fig)

with (RESULTS / "surface_threshold.csv").open("w") as handle:
    handle.write("d,p,logical_error_rate,shots\n")
    for d, rates in curves.items():
        for p, rate in zip(p_grid, rates):
            handle.write(f"{d},{p:.4f},{rate:.6f},{shots[d]}\n")
print("Wrote figures/surface_threshold.png, results/surface_threshold.csv")
