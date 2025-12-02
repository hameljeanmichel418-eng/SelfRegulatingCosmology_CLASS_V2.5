#!/usr/bin/env python3
"""
Compute effective stability parameter kappa_eff(z) for SRC gamma0 scan.

We use the already computed ΔH/H[%] from:
  output/src/deltaH_gamma_scan_z0_1.dat

Definitions:
  - ΔH/H is stored in percent (%)
  - We convert to a fraction: delta_frac = (ΔH/H[%]) / 100
  - kappa_eff(z) = d(delta_frac)/d ln a
  - With a = 1/(1+z), ln a = -ln(1+z) -> d ln a = - dz/(1+z)
    => d/d ln a = -(1+z) d/dz
  - Hence:
      kappa_eff(z) = -(1+z) * d(delta_frac)/dz
                   = -(1+z)/100 * d(ΔH/H[%])/dz

We compute kappa_eff(z) for runs:
  A2 (gamma0=0.020)
  A3 (gamma0=0.026)
  A4 (gamma0=0.030)
"""

import numpy as np
import os


def load_deltaH_table(filename="output/src/deltaH_gamma_scan_z0_1.dat"):
    """
    Load zoomed ΔH/H[%] table for the SRC gamma scan.

    Columns:
      0 : z
      1 : ΔH/H_A1[%]
      2 : ΔH/H_A2[%]
      3 : ΔH/H_A3[%]
      4 : ΔH/H_A4[%]
      5 : ΔH/H_A5[%]

    Returns
    -------
    z : ndarray
        Redshift array.
    delta_dict : dict
        Mapping tag -> ΔH/H[%] array for that run.
    """
    if not os.path.exists(filename):
        raise FileNotFoundError(f"Could not find file: {filename}")

    data = np.loadtxt(filename, comments="#")

    z = data[:, 0]
    dA1 = data[:, 1]
    dA2 = data[:, 2]
    dA3 = data[:, 3]
    dA4 = data[:, 4]
    dA5 = data[:, 5]

    # Just in case: sort by increasing z
    sort_indices = np.argsort(z)
    z_sorted = z[sort_indices]

    delta_dict = {
        "A1": dA1[sort_indices],
        "A2": dA2[sort_indices],
        "A3": dA3[sort_indices],
        "A4": dA4[sort_indices],
        "A5": dA5[sort_indices],
    }

    return z_sorted, delta_dict


def compute_kappa(z, delta_percent):
    """
    Compute kappa_eff(z) from ΔH/H[%] using finite differences.

    Parameters
    ----------
    z : ndarray
        Redshift array (assumed increasing).
    delta_percent : ndarray
        ΔH/H in percent (%).

    Returns
    -------
    z_sorted : ndarray
        Redshift array (increasing, same as input).
    kappa : ndarray
        kappa_eff(z), dimensionless.
    """
    # Ensure increasing z
    sort_indices = np.argsort(z)
    z_sorted = z[sort_indices]
    d_sorted = delta_percent[sort_indices]

    n = len(z_sorted)
    if n < 3:
        raise ValueError("Need at least 3 points to compute finite differences.")

    # Finite difference for d(ΔH/H[%])/dz
    d_delta_dz = np.zeros(n)

    # Central differences for interior points
    for i in range(1, n - 1):
        dz = z_sorted[i + 1] - z_sorted[i - 1]
        d_delta = d_sorted[i + 1] - d_sorted[i - 1]
        d_delta_dz[i] = d_delta / dz

    # Boundary: copy neighbors
    d_delta_dz[0] = d_delta_dz[1]
    d_delta_dz[-1] = d_delta_dz[-2]

    # kappa_eff(z) = -(1+z)/100 * d(ΔH/H[%])/dz
    kappa = -(1.0 + z_sorted) * d_delta_dz / 100.0

    return z_sorted, kappa


def main():
    """Main driver for kappa_eff computation."""

    print("Computing kappa_eff(z) for SRC gamma0 scan (A2, A3, A4)")
    print("=" * 70)

    # Load ΔH/H table
    z, delta_dict = load_deltaH_table()
    print(f"Loaded ΔH/H table with {len(z)} points in z ∈ [{z.min():.3f}, {z.max():.3f}]")

    # Runs of interest
    runs = [
        ("A2", 0.020),
        ("A3", 0.026),
        ("A4", 0.030),
    ]

    # Storage for summary
    results = {}

    # Compute kappa for each run
    for tag, gamma0 in runs:
        print(f"\nProcessing {tag} (gamma0={gamma0:.3f})...")
        delta_percent = delta_dict[tag]
        z_k, kappa_k = compute_kappa(z, delta_percent)

        results[tag] = {"z": z_k, "kappa": kappa_k, "gamma0": gamma0}

        # Save full kappa(z) curve
        outfile = f"output/src/kappa_{tag}.dat"
        header = (
            "# z   kappa_eff(z)\n"
            f"# tag = {tag}, gamma0 = {gamma0:.3f}\n"
            "# kappa_eff is dimensionless (uses ΔH/H as a fraction)."
        )
        out_array = np.column_stack([z_k, kappa_k])
        np.savetxt(outfile, out_array, header=header, fmt="%.10e")
        print(f"  Saved full kappa(z) curve to: {outfile}")

    # Build summary table at key redshifts
    print("\n" + "=" * 70)
    print("Building summary table at key redshifts...")
    print("=" * 70)

    z_targets = np.array([0.0, 0.25, 0.5, 0.75, 1.0])

    kA2_t = np.interp(z_targets, results["A2"]["z"], results["A2"]["kappa"])
    kA3_t = np.interp(z_targets, results["A3"]["z"], results["A3"]["kappa"])
    kA4_t = np.interp(z_targets, results["A4"]["z"], results["A4"]["kappa"])

    summary_array = np.column_stack([z_targets, kA2_t, kA3_t, kA4_t])

    summary_file = "output/src/kappa_gamma_points.txt"
    header_summary = (
        "# z   kappa_A2   kappa_A3   kappa_A4\n"
        "# gamma0 values: 0.020 0.026 0.030\n"
        "# kappa_eff is dimensionless (ΔH/H treated as a fraction)."
    )
    np.savetxt(summary_file, summary_array, header=header_summary, fmt="%.10e")
    print(f"Summary saved to: {summary_file}")

    # Print summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Redshift targets: {z_targets}")
    print()

    gamma_values = [0.020, 0.026, 0.030]
    run_tags = ["A2", "A3", "A4"]
    k_arrays = [kA2_t, kA3_t, kA4_t]

    for tag, gamma0, k_t in zip(run_tags, gamma_values, k_arrays):
        k_str = ", ".join([f"{val:+.6e}" for val in k_t])
        print(f"{tag} (gamma0={gamma0:.3f}): kappa_eff(z) = [{k_str}]")

    print("=" * 70)


if __name__ == "__main__":
    main()

