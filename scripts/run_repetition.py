"""Repetition code: Monte Carlo vs the exact analytic failure formula.

Produces figures/repetition_code.png and results/repetition.csv.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from qec.repetition import repetition_logical_error_analytic, simulate_repetition

ROOT = Path(__file__).resolve().parents[1]
FIGURES = ROOT / "figures"
RESULTS = ROOT / "results"
FIGURES.mkdir(exist_ok=True)
RESULTS.mkdir(exist_ok=True)

distances = [1, 3, 5, 7, 9]
p_grid = np.linspace(0.02, 0.9, 45)
shots = 20_000

print("Repetition code: analytic (lines) vs Monte Carlo (points) ...")
fig, ax = plt.subplots(figsize=(7.5, 4.8))
rows = []
for d in distances:
    analytic = [repetition_logical_error_analytic(d, p) for p in p_grid]
    ax.plot(p_grid, analytic, "-", label=f"d = {d} (exact)")
    if d > 1:
        mc = [simulate_repetition(d, p, shots, seed=100 + d) for p in p_grid[::6]]
        ax.plot(p_grid[::6], mc, "o", ms=4, color="k",
                label="Monte Carlo" if d == 3 else None)
    rows.extend((d, float(p), repetition_logical_error_analytic(d, p))
                for p in p_grid)

ax.axvline(0.5, color="C3", ls=":", lw=1.2, label="threshold $p$ = 0.5")
ax.set_xlabel("physical flip probability $p$")
ax.set_ylabel("logical error rate $p_L$")
ax.set_title("Repetition code: below 50% bigger codes win, above they lose")
ax.legend(ncol=2, fontsize=9)
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(FIGURES / "repetition_code.png", dpi=150)
plt.close(fig)

with (RESULTS / "repetition.csv").open("w") as handle:
    handle.write("d,p,logical_error_analytic\n")
    for d, p, rate in rows:
        handle.write(f"{d},{p:.4f},{rate:.6f}\n")
print("Wrote figures/repetition_code.png, results/repetition.csv")
