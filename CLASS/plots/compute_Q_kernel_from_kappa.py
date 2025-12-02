#!/usr/bin/env python3
"""
Build properly normalized SRC energy–transfer kernel from kappa_eff(z)
for the fiducial SRC run A3 (gamma0 = 0.026).

Input
-----
    output/src/src_kernel_A3_normalized.dat

    Columns (from previous step A.6):
        z            : redshift
        a            : scale factor = 1/(1+z)
        kappa_eff    : effective stability parameter kappa_eff(z)
        kernel_norm  : kappa_eff(z) / max(kappa_eff) over 0 <= z <= 1

Output
------
    output/src/Q_kernel_A3.dat
        Columns:
            z
            a
            kappa_eff
            kernel_norm          (raw peak-normalized kernel)
            K_ln_norm            (Q-kernel normalized such that ∫ K_ln_norm d ln a = 1)

    output/src/Q_kernel_A3_points.txt
        Values of K_ln_norm at key redshifts:
            z = [0.0, 0.25, 0.5, 0.75, 1.0]
"""

import numpy as np


def load_src_kernel(filename):
    """
    Load the fiducial SRC kernel for run A3 from the normalized kappa file.

    Parameters
    ----------
    filename : str
        Path to output/src/src_kernel_A3_normalized.dat

    Returns
    -------
    z : ndarray
        Redshift array
    a : ndarray
        Scale factor array
    kappa_eff : ndarray
        Effective kappa_eff(z)
    kernel_norm : ndarray
        Peak-normalized kernel: kappa_eff / max(kappa_eff)
    """
    data = np.loadtxt(filename, comments="#")

    z = data[:, 0]
    a = data[:, 1]
    kappa_eff = data[:, 2]
    kernel_norm = data[:, 3]

    return z, a, kappa_eff, kernel_norm


def normalize_kernel_in_ln_a(a, kernel_norm):
    """
    Normalize the kernel in ln(a) such that:

        ∫ K_ln_norm(a) d ln a = 1

    We work in a-grid order (increasing a), compute a finite-difference
    approximation to d ln a, then rescale the kernel.

    Parameters
    ----------
    a : ndarray
        Scale factor array (not necessarily sorted)
    kernel_norm : ndarray
        Peak-normalized kernel values on the same grid as a

    Returns
    -------
    a_sorted : ndarray
        Scale factor sorted in increasing order
    K_ln_norm_sorted : ndarray
        Kernel normalized in ln(a) such that ∫ K_ln_norm d ln a = 1
    """
    # Sort by increasing a for a well-defined integration in ln(a)
    sort_idx = np.argsort(a)
    a_sorted = a[sort_idx]
    K_sorted = kernel_norm[sort_idx]

    # Compute d ln a with finite differences
    ln_a = np.log(a_sorted)
    n = len(a_sorted)
    dln_a = np.zeros(n)

    # Interior points: central differences
    for i in range(1, n - 1):
        dln_a[i] = 0.5 * (ln_a[i + 1] - ln_a[i - 1])

    # Boundaries: one-sided differences
    dln_a[0] = ln_a[1] - ln_a[0]
    dln_a[-1] = ln_a[-1] - ln_a[-2]

    # Total area in ln(a)
    area = np.sum(K_sorted * dln_a)

    if area <= 0.0:
        raise ValueError(f"Non-positive area in ln(a) normalization: area = {area:.6e}")

    # Normalized kernel in ln(a)
    K_ln_norm_sorted = K_sorted / area

    # Sanity check: recompute integral
    check_area = np.sum(K_ln_norm_sorted * dln_a)
    print(f"  Normalization check: ∫ K_ln_norm d ln a = {check_area:.6e}")

    return a_sorted, K_ln_norm_sorted


def main():
    """
    Main function to build the Q-kernel for the fiducial SRC run A3.
    """
    print("Building Q-kernel (normalized in ln a) from kappa_eff(z) for A3")
    print("=" * 70)

    # ------------------------------------------------------------------
    # 1) Load the fiducial SRC kernel from previous step A.6
    # ------------------------------------------------------------------
    src_kernel_file = "output/src/src_kernel_A3_normalized.dat"
    print(f"Loading SRC kernel from: {src_kernel_file}")
    z, a, kappa_eff, kernel_norm = load_src_kernel(src_kernel_file)

    print(f"Loaded {len(z)} points")
    print(f"Redshift range: z ∈ [{z.min():.6f}, {z.max():.6f}]")
    print(f"Scale-factor range: a ∈ [{a.min():.6f}, {a.max():.6f}]")

    # ------------------------------------------------------------------
    # 2) Normalize the kernel in ln(a)
    # ------------------------------------------------------------------
    print("\nNormalizing kernel in ln(a) such that ∫ K_ln_norm d ln a = 1 ...")
    a_sorted, K_ln_norm_sorted = normalize_kernel_in_ln_a(a, kernel_norm)

    # We also want the corresponding z and other columns reordered consistently.
    # We sorted by 'a', so we reuse that ordering to sort z, kappa_eff, kernel_norm.
    sort_idx = np.argsort(a)
    z_sorted = z[sort_idx]
    kappa_sorted = kappa_eff[sort_idx]
    kernel_raw_sorted = kernel_norm[sort_idx]

    print("  Kernel normalization in ln(a) completed.")

    # ------------------------------------------------------------------
    # 3) Save full Q-kernel table (sorted by increasing z for readability)
    # ------------------------------------------------------------------
    print("\nSaving full Q-kernel table ...")

    # For output, sort by increasing z (more natural for cosmology plots)
    sort_z_idx = np.argsort(z_sorted)
    z_out = z_sorted[sort_z_idx]
    a_out = a_sorted[sort_z_idx]
    kappa_out = kappa_sorted[sort_z_idx]
    kernel_raw_out = kernel_raw_sorted[sort_z_idx]
    K_ln_norm_out = K_ln_norm_sorted[sort_z_idx]

    Q_kernel_file = "output/src/Q_kernel_A3.dat"
    header = (
        "# z   a   kappa_eff   kernel_norm   K_ln_norm\n"
        "#\n"
        "# Q-kernel for SRC fiducial run A3 (gamma0 = 0.026).\n"
        "# kernel_norm  : peak-normalized (kappa_eff / max(kappa_eff)) in 0 <= z <= 1.\n"
        "# K_ln_norm    : normalized such that ∫ K_ln_norm(a) d ln a = 1.\n"
        "# This kernel encodes the shape of the late-time PBH→phi energy-transfer window.\n"
    )

    Q_array = np.column_stack([z_out, a_out, kappa_out, kernel_raw_out, K_ln_norm_out])
    np.savetxt(Q_kernel_file, Q_array, header=header, fmt="%.10e")
    print(f"  Saved Q-kernel table to: {Q_kernel_file}")

    # ------------------------------------------------------------------
    # 4) Save key points at specific redshifts (for Volume 2.5 tables)
    # ------------------------------------------------------------------
    print("\nExtracting key Q-kernel values at z = 0, 0.25, 0.5, 0.75, 1.0 ...")
    z_targets = np.array([0.0, 0.25, 0.5, 0.75, 1.0])

    # We interpolate K_ln_norm_out as a function of z_out
    K_targets = np.interp(z_targets, z_out, K_ln_norm_out)

    Q_points_file = "output/src/Q_kernel_A3_points.txt"
    header_points = (
        "# z   K_ln_norm(z)\n"
        "# Q-kernel for A3 (gamma0 = 0.026), normalized in ln a.\n"
        "# Defined such that ∫ K_ln_norm(a) d ln a = 1 over the range 0 <= z <= 1.\n"
    )

    points_array = np.column_stack([z_targets, K_targets])
    np.savetxt(Q_points_file, points_array, header=header_points, fmt="%.10e")
    print(f"  Saved key points to: {Q_points_file}")

    # ------------------------------------------------------------------
    # 5) Print short summary
    # ------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Redshift targets: {z_targets}")
    print("K_ln_norm(z) at targets:")
    for z_t, K_t in zip(z_targets, K_targets):
        print(f"  z = {z_t:4.2f}  ->  K_ln_norm = {K_t:.6e}")
    print("=" * 70)


if __name__ == "__main__":
    main()

