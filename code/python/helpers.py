### This file is fully ChatGPT generated for analytical purposes.
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt


# ============================================================================ #
# Readers
def read_output_file(filepath):
    filepath = Path(filepath)
    with open(filepath, "r") as file:
        header = file.readline().strip()

    if not header.startswith("#"):
        raise ValueError(f"No header found in {filepath}")

    columns = header[1:].split()
    return pd.read_csv(filepath, sep=r"\s+", comment="#", names=columns)


# ============================================================================ #
def get_run_path(path, N):
    return Path(path) / f"N_{N}"


# ============================================================================ #
from pathlib import Path
import re


# ============================================================================ #
def read_parameters_file(file_path):
    file_path = Path(file_path)

    params = {
        "seed": None,
        "elastic": {},
        "cohesive": {},
        "user": {},
    }

    section = None

    with open(file_path, "r") as file:
        for line in file:
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            # Seed
            if line.startswith("seed"):
                params["seed"] = int(line.split("=")[1].strip())
                continue

            # Sections
            if line.startswith("material elastic"):
                section = "elastic"
                continue

            if line.startswith("material cohesive_linear"):
                section = "cohesive"
                continue

            if line.startswith("user parameters"):
                section = "user"
                continue

            if line == "]":
                section = None
                continue

            if section is None or "=" not in line:
                continue

            key, value = line.split("=", 1)
            key = key.strip()
            value = value.strip()

            # Special case:
            # sigma_c = 262000000 uniform [0, 1000000]
            uniform_match = re.match(
                r"^([+-]?[0-9.eE+-]+)\s+uniform\s+\[([+-]?[0-9.eE+-]+),\s*([+-]?[0-9.eE+-]+)\]$",
                value,
            )

            if uniform_match:
                params[section][key] = float(uniform_match.group(1))
                params[section][f"{key}_distribution"] = "uniform"
                params[section][f"{key}_distribution_min"] = float(
                    uniform_match.group(2)
                )
                params[section][f"{key}_distribution_max"] = float(
                    uniform_match.group(3)
                )
                continue

            # Boolean
            if value.lower() == "true":
                params[section][key] = True
                continue

            if value.lower() == "false":
                params[section][key] = False
                continue

            # Number
            try:
                number = float(value)

                if number.is_integer():
                    number = int(number)

                params[section][key] = number
                continue

            except ValueError:
                pass

            # String
            params[section][key] = value

    return params


# ============================================================================ #
def characteristic_strain_rate(E, rho, sigma_c, G_c):
    c = np.sqrt(E / rho)
    eps_dot_0 = sigma_c**3 * c / (E**2 * G_c)

    return eps_dot_0


# ============================================================================ #
def characteristic_fragment_size(E, sigma_c, G_c):
    s_0 = E * G_c / sigma_c**2

    return s_0


# ============================================================================ #
def grady_fragment_size(rho, G_c, eps_dot):
    s = (24.0 * G_c / (rho * eps_dot**2)) ** (1.0 / 3.0)

    return s


# ============================================================================ #
def zmr_fragment_size(E, rho, sigma_c, G_c, eps_dot):
    eps_dot_0 = characteristic_strain_rate(E, rho, sigma_c, G_c)
    s_0 = characteristic_fragment_size(E, sigma_c, G_c)

    eps_dot_norm = eps_dot / eps_dot_0

    s_norm = 4.5 / (1.0 + 4.5 * eps_dot_norm ** (2.0 / 3.0))
    s = s_norm * s_0

    return s


# ============================================================================ #
def normalised_fragmentation_values(path, N, E, rho, sigma_c, G_c):
    run_path = get_run_path(path, N)

    if not (run_path / "energies.out").exists():
        data = read_output_file(run_path / "energies_pf.out")
    else:
        data = read_output_file(run_path / "energies.out")

    # Characteristic quantities
    eps_dot_0 = characteristic_strain_rate(E, rho, sigma_c, G_c)
    s_0 = characteristic_fragment_size(E, sigma_c, G_c)

    # Applied strain rate obtained from the simulation output
    eps_dot = (data["strain"].iloc[-1] - data["strain"].iloc[0]) / (
        data["time"].iloc[-1] - data["time"].iloc[0]
    )

    # Final average fragment size
    s = data["mean_fragment_size"].iloc[-1]

    # Normalised quantities
    eps_dot_norm = eps_dot / eps_dot_0
    s_norm = s / s_0

    return eps_dot_norm, s_norm


# ============================================================================ #
def locate_final_fragment_time(path, N):

    run_path = get_run_path(path, N)

    if (run_path / "energies.out").exists():
        data = read_output_file(run_path / "energies.out")
    elif (run_path / "energies_pf.out").exists():
        data = read_output_file(run_path / "energies_pf.out")
    else:
        raise FileNotFoundError(f"No energy output found in {run_path}")

    dt = data["dt"].iloc[0]

    log_file = run_path / f"run_{N}.log"

    if not log_file.exists():
        log_file = run_path / f"run_{N}.out"

    if not log_file.exists():
        raise FileNotFoundError(f"No run log found in {run_path}")

    steps = []
    fragments = []

    with open(log_file, "r") as file:
        for line in file:
            if line.startswith("Step:"):
                parts = line.split("--")

                step = int(parts[0].split(":")[1].strip())
                num_fragments = int(parts[1].split("=")[1].strip())

                steps.append(step)
                fragments.append(num_fragments)

    if not fragments:
        raise ValueError(f"No fragment-count checks found in {log_file}")

    final_num_fragments = fragments[-1]

    i = len(fragments) - 1

    while i > 0 and fragments[i - 1] == final_num_fragments:
        i -= 1

    t_frag = steps[i] * dt

    if len(steps) > 1:
        sampling_interval = (steps[1] - steps[0]) * dt
    else:
        sampling_interval = np.nan

    return t_frag, sampling_interval


# ============================================================================ #
def peak_average_stress(path, N):
    run_path = get_run_path(path, N)

    if (run_path / "energies.out").exists():
        data = read_output_file(run_path / "energies.out")
    elif (run_path / "energies_pf.out").exists():
        data = read_output_file(run_path / "energies_pf.out")
    else:
        raise FileNotFoundError(f"No energy output found in {run_path}")

    return data["average_stress"].max()


# ============================================================================ #
def fracture_energy_at_time(path, N, time):
    run_path = get_run_path(path, N)

    if (run_path / "energies.out").exists():
        data = read_output_file(run_path / "energies.out")
    elif (run_path / "energies_pf.out").exists():
        data = read_output_file(run_path / "energies_pf.out")
    else:
        raise FileNotFoundError(f"No energy output found in {run_path}")

    fracture_energy = np.interp(time, data["time"], data["fracture_energy"])

    return fracture_energy


# ============================================================================ #
def max_energy_error_until_time(path, N, time):
    run_path = get_run_path(path, N)

    if (run_path / "energies.out").exists():
        data = read_output_file(run_path / "energies.out")
    elif (run_path / "energies_pf.out").exists():
        data = read_output_file(run_path / "energies_pf.out")
    else:
        raise FileNotFoundError(f"No energy output found in {run_path}")

    data_until_time = data[data["time"] <= time]

    if data_until_time.empty:
        raise ValueError(f"No diagnostic data available before t = {time}")

    return np.abs(data_until_time["energy_error"]).max()


# ============================================================================ #
def bar_height(path, N):
    run_path = get_run_path(path, N)

    params = read_parameters_file(run_path / "parameters_used.dat")

    L = params["user"]["bar_length"]
    N_elements = params["user"]["number_elements"]

    return 2.0 * L / N_elements


# ============================================================================ #
def mesh_convergence_values(path, N):
    t_frag, sampling_interval = locate_final_fragment_time(path, N)

    sigma_peak = peak_average_stress(path, N)

    fracture_energy = fracture_energy_at_time(path, N, t_frag)
    energy_error_max = max_energy_error_until_time(path, N, t_frag)

    H = bar_height(path, N)

    fracture_energy_per_height = fracture_energy / H
    energy_error_per_height = energy_error_max / H

    return (
        sigma_peak,
        t_frag,
        sampling_interval,
        fracture_energy_per_height,
        energy_error_per_height,
    )


# ============================================================================ #
def fragment_evolution(path, N):
    run_path = get_run_path(path, N)
    data = read_output_file(run_path / "energies.out")

    time = data["time"].to_numpy()
    num_fragments = data["num_fragments"].to_numpy()

    return time, num_fragments


# ============================================================================ #
def energy_evolution(path, N):
    run_path = get_run_path(path, N)

    data = read_output_file(run_path / "energies.out")
    params = read_parameters_file(run_path / "parameters_used.dat")

    L = params["user"]["bar_length"]
    N_elements = params["user"]["number_elements"]

    H = 2.0 * L / N_elements

    time = data["time"].to_numpy()

    energies = {
        "kinetic": data["kinetic_energy"].to_numpy() / H,
        "elastic": data["elastic_energy"].to_numpy() / H,
        "reversible": data["reversible_cohesive_energy"].to_numpy() / H,
        "fracture": data["fracture_energy"].to_numpy() / H,
        "contact": data["contact_energy"].to_numpy() / H,
        "external_work": data["external_work"].to_numpy() / H,
        "total": data["total_energy"].to_numpy() / H,
        "energy_error": data["energy_error"].to_numpy() / H,
    }

    return time, energies


# ============================================================================ #
def plot_energy_evolution(path, N, time_cut=None):
    time, energies = energy_evolution(path, N)

    if time_cut is not None:
        mask = time <= time_cut
        time = time[mask]

        for key in energies:
            energies[key] = energies[key][mask]

    plt.figure()

    plt.plot(time * 1e6, energies["kinetic"], label="Kinetic")
    plt.plot(time * 1e6, energies["elastic"], label="Elastic")
    plt.plot(time * 1e6, energies["reversible"], label="Reversible cohesive")
    plt.plot(time * 1e6, energies["fracture"], label="Fracture")
    plt.plot(time * 1e6, energies["contact"], label="Contact")
    plt.plot(time * 1e6, energies["external_work"], label="External work")
    plt.plot(time * 1e6, energies["total"], label="Total", linestyle="--")

    plt.xlabel(r"Time [$\mu$s]")
    plt.ylabel("Energy / bar height")
    plt.legend()
    plt.tight_layout()
    plt.show()


# ============================================================================ #
def plot_energy_for_meshes(path, Ns, energy_type, time_cut=None):
    plt.figure()

    for N in Ns:
        time, energies = energy_evolution(path, N)

        if energy_type not in energies:
            raise ValueError(
                f"Unknown energy type '{energy_type}'. "
                f"Choose from: {list(energies.keys())}"
            )

        energy = energies[energy_type]

        if time_cut is not None:
            mask = time <= time_cut
            time = time[mask]
            energy = energy[mask]

        plt.plot(time * 1e6, energy, label=f"N = {N}")

    plt.xlabel(r"Time [$\mu$s]")
    plt.ylabel(f"{energy_type.replace('_', ' ').title()} / bar height")
    plt.legend()
    plt.tight_layout()
    plt.show()


# ============================================================================
# Phase-field analysis
# ============================================================================


def get_pf_case_path(path, N, sweep=None, value=None):
    """
    Return the directory of a phase-field simulation.

    Examples
    --------
    Mesh convergence:
        get_pf_case_path(path, 10000)

    Time-step sweep:
        get_pf_case_path(path, 10000, "time_step", "0.25")

    l0 sweep:
        get_pf_case_path(path, 10000, "l0", "1e-4")

    Damage-limit sweep:
        get_pf_case_path(path, 10000, "damage_limit", "0.95")

    Strain-rate sweep:
        get_pf_case_path(path, 10000, "strain_rate", "1e4")
    """

    run_path = Path(path) / f"N_{N}"

    if sweep is None:
        return run_path

    prefixes = {
        "time_step": "dt",
        "l0": "l0",
        "damage_limit": "damage",
        "strain_rate": "rate",
    }

    if sweep not in prefixes:
        raise ValueError(
            f"Unknown PF sweep '{sweep}'. " f"Choose from {list(prefixes.keys())}"
        )

    return run_path / f"{prefixes[sweep]}_{value}"


# ============================================================================


def read_scalar_parameter(file_path, key):
    """
    Read a scalar parameter from material.dat / parameters_used.dat.

    Useful for material parameters such as l0 that are not necessarily
    inside the user-parameters block.
    """

    import re

    file_path = Path(file_path)

    pattern = re.compile(rf"^\s*{re.escape(key)}\s*=\s*([^\s#;]+)")

    with open(file_path, "r") as file:
        for line in file:
            match = pattern.search(line)

            if match:
                value = match.group(1)

                try:
                    return float(value)
                except ValueError:
                    return value

    raise KeyError(f"Parameter '{key}' not found in {file_path}")


# ============================================================================ #

def pf_run_values(run_path):
    """
    Extract the main quantities required to analyse one phase-field run.
    """

    run_path = Path(run_path)

    params = read_parameters_file(run_path / "parameters_used.dat")

    L = params["user"]["bar_length"]
    N = int(params["user"]["number_elements"])

    H = 2.0 * L / N

    diagnostics_file = params["user"].get("diagnostics_file", "energies.out")

    data = read_output_file(run_path / diagnostics_file)

    # ------------------------------------------------------------------------
    # Fragmentation

    final_num_fragments = int(data["num_fragments"].iloc[-1])

    final_mean_fragment_size = data["mean_fragment_size"].iloc[-1]
    final_min_fragment_size = data["min_fragment_size"].iloc[-1]
    final_max_fragment_size = data["max_fragment_size"].iloc[-1]

    fragmentation = data[data["num_fragments"] > 1]

    if len(fragmentation) > 0:
        fragmentation_time = fragmentation["time"].iloc[0]
    else:
        fragmentation_time = np.nan

    # ------------------------------------------------------------------------
    # Stress

    peak_average_stress = data["average_stress"].max()
    peak_max_stress = data["max_stress"].max()

    # ------------------------------------------------------------------------
    # Energy

    final_elastic_energy = data["elastic_energy"].iloc[-1] / H
    final_kinetic_energy = data["kinetic_energy"].iloc[-1] / H
    final_phase_energy = data["fracture_energy"].iloc[-1] / H
    final_total_energy = data["total_energy"].iloc[-1] / H
    final_external_work = data["external_work"].iloc[-1] / H

    max_energy_error = np.abs(data["energy_error"]).max() / H

    # ------------------------------------------------------------------------

    return {
        "N": N,
        "mesh_size": H,
        "bar_height": H,
        "final_time": data["time"].iloc[-1],
        "num_fragments": final_num_fragments,
        "mean_fragment_size": final_mean_fragment_size,
        "min_fragment_size": final_min_fragment_size,
        "max_fragment_size": final_max_fragment_size,
        "fragmentation_time": fragmentation_time,
        "peak_average_stress": peak_average_stress,
        "peak_max_stress": peak_max_stress,
        "elastic_energy_per_height": final_elastic_energy,
        "kinetic_energy_per_height": final_kinetic_energy,
        "phase_energy_per_height": final_phase_energy,
        "total_energy_per_height": final_total_energy,
        "external_work_per_height": final_external_work,
        "max_energy_error_per_height": max_energy_error,
    }


# ============================================================================ #
# PF mesh convergence


def pf_mesh_convergence_values(path, Ns):
    rows = []

    for N in Ns:
        run_path = get_pf_case_path(path, N)

        values = pf_run_values(run_path)

        rows.append(values)

    return pd.DataFrame(rows)


# ============================================================================ #
# PF time-step convergence


def pf_time_step_convergence_values(path, N, time_step_factors):
    rows = []

    for factor in time_step_factors:
        run_path = get_pf_case_path(path, N, "time_step", factor)

        values = pf_run_values(run_path)

        values["time_step_factor"] = float(factor)

        rows.append(values)

    return pd.DataFrame(rows)


# ============================================================================ #
# PF l0 sweep

def pf_l0_sweep_values(path, N, l0_values):
    rows = []

    for l0 in l0_values:
        run_path = get_pf_case_path(path, N, "l0", l0)

        values = pf_run_values(run_path)

        l0_float = float(l0)

        values["l0"] = l0_float
        values["h_over_l0"] = values["mesh_size"] / l0_float

        rows.append(values)

    return pd.DataFrame(rows)


# ============================================================================ #
# PF damage-limit sweep


def pf_damage_limit_sweep_values(path, N, damage_limits):
    rows = []

    for damage_limit in damage_limits:
        run_path = get_pf_case_path(path, N, "damage_limit", damage_limit)

        values = pf_run_values(run_path)

        values["damage_limit"] = float(damage_limit)

        rows.append(values)

    return pd.DataFrame(rows)


# ============================================================================ #
# PF strain-rate sweep


def pf_strain_rate_sweep_values(path, N, strain_rates):
    rows = []

    for strain_rate in strain_rates:
        run_path = get_pf_case_path(path, N, "strain_rate", strain_rate)

        values = pf_run_values(run_path)

        values["strain_rate"] = float(strain_rate)

        rows.append(values)

    return pd.DataFrame(rows)


# ============================================================================ #
# PF fragment evolution


def pf_fragment_evolution(run_path):
    run_path = Path(run_path)

    params = read_parameters_file(run_path / "parameters_used.dat")

    diagnostics_file = params["user"].get("diagnostics_file", "energies.out")

    data = read_output_file(run_path / diagnostics_file)

    time = data["time"].to_numpy()
    num_fragments = data["num_fragments"].to_numpy()
    mean_fragment_size = data["mean_fragment_size"].to_numpy()

    return time, num_fragments, mean_fragment_size


# ============================================================================ #
# PF energy evolution


def pf_energy_evolution(run_path):
    run_path = Path(run_path)

    params = read_parameters_file(run_path / "parameters_used.dat")

    L = params["user"]["bar_length"]
    N = int(params["user"]["number_elements"])

    H = 2.0 * L / N

    diagnostics_file = params["user"].get("diagnostics_file", "energies.out")

    data = read_output_file(run_path / diagnostics_file)

    time = data["time"].to_numpy()

    energies = {
        "kinetic": data["kinetic_energy"].to_numpy() / H,
        "elastic": data["elastic_energy"].to_numpy() / H,
        "phase": data["fracture_energy"].to_numpy() / H,
        "external_work": data["external_work"].to_numpy() / H,
        "total": data["total_energy"].to_numpy() / H,
        "energy_error": data["energy_error"].to_numpy() / H,
    }

    return time, energies


# ============================================================================ #
# PF stress evolution


def pf_stress_evolution(run_path):
    run_path = Path(run_path)

    params = read_parameters_file(run_path / "parameters_used.dat")

    diagnostics_file = params["user"].get("diagnostics_file", "energies.out")

    data = read_output_file(run_path / diagnostics_file)

    time = data["time"].to_numpy()
    average_stress = data["average_stress"].to_numpy()
    max_stress = data["max_stress"].to_numpy()

    return time, average_stress, max_stress
