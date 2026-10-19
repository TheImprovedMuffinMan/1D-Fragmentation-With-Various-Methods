import sys
from pathlib import Path
import matplotlib.pyplot as plt

sys.path.append(str(Path(__file__).resolve().parents[1]))

from python import helpers

PROJECT_DIR = Path(__file__).resolve().parents[2]
DATA_DIR = PROJECT_DIR / "data" / "CZM" / "mesh_convergence"


# ============================================================================ #
if __name__ == "__main__":

    Ns = [1000, 1500, 2000, 3000, 5000, 7500, 10000, 15000, 20000, 30000, 50000]
    eps_dot_norms = []
    s_norms = []

    for N in Ns:

        run_path = helpers.get_run_path(DATA_DIR, N)

        params = helpers.read_parameters_file(run_path / "parameters_used.dat")

        rho = params["elastic"]["rho"]
        E = params["elastic"]["E"]

        sigma_c = params["cohesive"]["sigma_c"]
        G_c = params["cohesive"]["G_c"]

        L = params["user"]["bar_length"]
        eps_dot_input = params["user"]["strain_rate"]
        N_file = params["user"]["number_elements"]

        eps_dot_0 = helpers.characteristic_strain_rate(E, rho, sigma_c, G_c)

        s_0 = helpers.characteristic_fragment_size(E, sigma_c, G_c)

        eps_dot_norm, s_norm = helpers.normalised_fragmentation_values(
            DATA_DIR, N, E, rho, sigma_c, G_c
        )

        eps_dot_norms.append(eps_dot_norm)
        s_norms.append(s_norm)

        s_grady = helpers.grady_fragment_size(rho, G_c, eps_dot_input)
        s_zmr = helpers.zmr_fragment_size(E, rho, sigma_c, G_c, eps_dot_input)

        s_grady_norm = s_grady / s_0
        s_zmr_norm = s_zmr / s_0

        sigma_peak, t_frag, t_frag_error, E_fracture_per_height, E_error_per_height = (
            helpers.mesh_convergence_values(DATA_DIR, N)
        )

        print(f"N = {N_file:6d}")
        print(f"  eps_dot/eps_dot_0         = {eps_dot_norm:.6f}")
        print(f"  s/s_0                     = {s_norm:.6f}")
        print(
            f"  Grady: s/s_0              = {s_grady_norm:.6f}, s = {s_grady * 1e3:.6f} mm"
        )
        print(
            f"  ZMR:   s/s_0              = {s_zmr_norm:.6f}, s = {s_zmr * 1e3:.6f} mm"
        )
        print(f"  peak stress               = {sigma_peak / 1e6:.6f} MPa")
        print(f"  t_frag                    = {t_frag * 1e6:.6f} us")
        print(f"  fracture energy / height  = {E_fracture_per_height:.6e}")
        print(f"  max energy error / height = {E_error_per_height:.6e}")
        print()
