import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.append(str(Path(__file__).resolve().parents[1]))
from python import helpers

# ============================================================================ #
# plot parameters

# plotting parameters
plt.rcParams.update(
    {
        # Figure
        "figure.figsize": (9, 6),
        "figure.dpi": 300,
        # Font
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Times", "Nimbus Roman"],
        "font.size": 14,
        # Titles and labels
        "axes.titlesize": 20,
        "axes.titleweight": "bold",
        "axes.labelsize": 15,
        # Tick labels
        "xtick.labelsize": 14,
        "ytick.labelsize": 14,
        # Lines and markers
        "lines.linewidth": 2,
        "lines.markersize": 6,
        # Legend
        "legend.fontsize": 15,
        "legend.frameon": True,
        # Grid
        "axes.grid": True,
        "grid.linestyle": "--",
        "grid.linewidth": 0.6,
        "grid.alpha": 0.5,
        # Axes
        "axes.linewidth": 1.0,
        # Saving figures
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        # Mathematical text
        "mathtext.fontset": "stix",
    }
)

# Colours I liked
COLORS = {
    "blue": "#5F8FAF",
    "orange": "#E6A15C",
    "red": "#D97A7A",
    "green": "#7FB685",
    "purple": "#A58CC2",
}


# ============================================================================ #
# Paths

PROJECT_DIR = Path(__file__).resolve().parents[2]

DATA_DIR = PROJECT_DIR / "data" / "CZM" / "mesh_convergence/"

Ns = [1000, 5000, 10000, 30000, 50000]

plt.figure()

time_cut = 0.25e-6

for N in Ns:
    time, num_fragments = helpers.fragment_evolution(DATA_DIR, N)
    mask = time <= time_cut
    plt.step(time[mask] * 1e6, num_fragments[mask], where="post", label=f"N = {N}")

plt.xlabel(r"Time [$\mu$s]")
plt.ylabel("Number of fragments")
plt.legend()
plt.tight_layout()
plt.show()
