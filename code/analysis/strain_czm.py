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

DATA_DIR_10000 = PROJECT_DIR / "data" / "CZM" / "strain_rate_sweep" / "N_10000"
DATA_DIR_40000 = PROJECT_DIR / "data" / "CZM" / "strain_rate_sweep" / "N_50000"


# ============================================================================ #
# Material parameters

first_run = next(DATA_DIR_10000.glob("rate_*"))
params = helpers.read_parameters_file(first_run / "parameters_used.dat")

E = params["elastic"]["E"]
rho = params["elastic"]["rho"]
sigma_c = params["cohesive"]["sigma_c"]
G_c = params["cohesive"]["G_c"]

eps_dot_0 = helpers.characteristic_strain_rate(E, rho, sigma_c, G_c)
s_0 = helpers.characteristic_fragment_size(E, sigma_c, G_c)


# ============================================================================ #
# N = 10000

strain_rates_10000 = []
fragment_sizes_10000 = []

for run_path in DATA_DIR_10000.glob("rate_*"):

    parameter_file = run_path / "parameters_used.dat"
    output_file = run_path / "energies.out"

    if not parameter_file.exists() or not output_file.exists():
        continue

    params = helpers.read_parameters_file(parameter_file)
    data = helpers.read_output_file(output_file)

    strain_rates_10000.append(params["user"]["strain_rate"])
    fragment_sizes_10000.append(data["mean_fragment_size"].iloc[-1])

strain_rates_10000 = np.asarray(strain_rates_10000)
fragment_sizes_10000 = np.asarray(fragment_sizes_10000)

order = np.argsort(strain_rates_10000)
strain_rates_10000 = strain_rates_10000[order]
fragment_sizes_10000 = fragment_sizes_10000[order]


# ============================================================================ #
# N = 40000

strain_rates_40000 = []
fragment_sizes_40000 = []

for run_path in DATA_DIR_40000.glob("rate_*"):

    parameter_file = run_path / "parameters_used.dat"
    output_file = run_path / "energies.out"

    if not parameter_file.exists() or not output_file.exists():
        continue

    params = helpers.read_parameters_file(parameter_file)
    data = helpers.read_output_file(output_file)

    strain_rates_40000.append(params["user"]["strain_rate"])
    fragment_sizes_40000.append(data["mean_fragment_size"].iloc[-1])

strain_rates_40000 = np.asarray(strain_rates_40000)
fragment_sizes_40000 = np.asarray(fragment_sizes_40000)

order = np.argsort(strain_rates_40000)
strain_rates_40000 = strain_rates_40000[order]
fragment_sizes_40000 = fragment_sizes_40000[order]


# ============================================================================ #
# Normalise CZM results

eps_norm_10000 = strain_rates_10000 / eps_dot_0
s_norm_10000 = fragment_sizes_10000 / s_0

eps_norm_40000 = strain_rates_40000 / eps_dot_0
s_norm_40000 = fragment_sizes_40000 / s_0


# ============================================================================ #
# Grady prediction

eps_dot_grady = np.logspace(
    np.log10(min(strain_rates_10000.min(), strain_rates_40000.min())),
    np.log10(max(strain_rates_10000.max(), strain_rates_40000.max())),
    300,
)

s_grady = np.array(
    [helpers.grady_fragment_size(rho, G_c, eps_dot) for eps_dot in eps_dot_grady]
)

eps_norm_grady = eps_dot_grady / eps_dot_0
s_norm_grady = s_grady / s_0


# ============================================================================ #
# Plot

plt.figure()
plt.loglog(
    eps_norm_grady,
    s_norm_grady,
    label="Grady",
    color=COLORS["green"],
    zorder=2,
    linewidth=2.5,
)
plt.scatter(
    eps_norm_10000,
    s_norm_10000,
    label="CZM, N = 10,000",
    color=COLORS["blue"],
    zorder=3,
)
plt.scatter(
    eps_norm_40000,
    s_norm_40000,
    label="CZM, N = 50,000",
    color=COLORS["orange"],
    zorder=3,
)

plt.xlabel(r"Normalised strain rate $\dot{\varepsilon}/\dot{\varepsilon}_0$")
plt.ylabel(r"Normalised mean fragment size $s/s_0$")
plt.legend()

plt.tight_layout()
plt.show()
