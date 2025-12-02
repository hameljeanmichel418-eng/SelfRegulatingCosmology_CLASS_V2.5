#!/usr/bin/env python3
"""
Reconstruct cumulative SRC injection fraction F(z) from normalized Q-kernel.

Computes F(z) = ∫[z to z_max] K_ln(z') / (1+z') dz' and compares with ΔH/H
for the fiducial A3 run.
"""

import numpy as np


def load_Q_kernel(filename="output/src/Q_kernel_A3.dat"):
    """
    Load the normalized Q-kernel for A3.
    
    Parameters
    ----------
    filename : str
        Path to Q_kernel_A3.dat
    
    Returns
    -------
    z : ndarray
        Redshift array
    a : ndarray
        Scale factor array
    K_ln : ndarray
        K_ln_norm, normalized such that ∫ K_ln d ln a = 1
    """
    data = np.loadtxt(filename, comments="#")
    
    z = data[:, 0]
    a = data[:, 1]
    K_ln = data[:, 4]  # K_ln_norm column
    
    # Sort by increasing z
    sort_idx = np.argsort(z)
    z_sorted = z[sort_idx]
    a_sorted = a[sort_idx]
    K_ln_sorted = K_ln[sort_idx]
    
    return z_sorted, a_sorted, K_ln_sorted


def compute_cumulative_F(z, K_ln):
    """
    Compute cumulative injection fraction F(z) from Q-kernel.
    
    F(z) = ∫[z to z_max] K_ln(z') / (1+z') dz'
    
    Normalized such that F(z_min) = 1.0, F(z_max) ≈ 0.0
    
    Parameters
    ----------
    z : ndarray
        Redshift array (sorted ascending)
    K_ln : ndarray
        K_ln_norm array on same grid
    
    Returns
    -------
    F_z : ndarray
        Cumulative fraction F(z), same length as z
    """
    n = len(z)
    if n < 2:
        raise ValueError("Need at least 2 points to compute cumulative integral")
    
    # Weight: K_ln / (1+z)
    weight = K_ln / (1.0 + z)
    
    # Compute cumulative integral from z to z_max (backwards from each point)
    # F(z_i) = ∫[z_i to z_max] weight(z') dz'
    F_z = np.zeros(n)
    
    # Use trapezoidal rule: integrate backwards from each z_i to z_max
    for i in range(n):
        # Integrate from z[i] to z[-1] (z_max)
        if i == n - 1:
            # At z_max, F = 0 (no injection beyond this point)
            F_z[i] = 0.0
        else:
            # Trapezoidal integration from z[i] to z[-1]
            z_segment = z[i:]
            weight_segment = weight[i:]
            
            # Trapezoidal rule
            dz_segments = np.diff(z_segment)
            avg_weights = 0.5 * (weight_segment[:-1] + weight_segment[1:])
            F_z[i] = np.sum(avg_weights * dz_segments)
    
    # Normalize so F(z_min) = 1.0
    if F_z[0] > 0:
        F_z = F_z / F_z[0]
    else:
        raise ValueError("F(z_min) is zero or negative, cannot normalize")
    
    return F_z


def main():
    """Main function to reconstruct cumulative Q profile."""
    
    print("Reconstructing cumulative SRC injection fraction F(z) from Q-kernel")
    print("="*70)
    
    # ------------------------------------------------------------------
    # 1) Load Q-kernel
    # ------------------------------------------------------------------
    Q_kernel_file = "output/src/Q_kernel_A3.dat"
    print(f"\nLoading Q-kernel from: {Q_kernel_file}")
    z, a, K_ln = load_Q_kernel(Q_kernel_file)
    
    print(f"Loaded {len(z)} points")
    print(f"Redshift range: z ∈ [{z.min():.6f}, {z.max():.6f}]")
    print(f"Scale-factor range: a ∈ [{a.min():.6f}, {a.max():.6f}]")
    
    # ------------------------------------------------------------------
    # 2) Compute cumulative F(z)
    # ------------------------------------------------------------------
    print("\nComputing cumulative injection fraction F(z)...")
    F_z = compute_cumulative_F(z, K_ln)
    
    print(f"F(z) range: [{F_z.min():.6f}, {F_z.max():.6f}]")
    print(f"  F(z_min) = {F_z[0]:.6f} (should be 1.0)")
    print(f"  F(z_max) = {F_z[-1]:.6f} (should be ≈ 0.0)")
    
    # ------------------------------------------------------------------
    # 3) Save full cumulative curve
    # ------------------------------------------------------------------
    print("\nSaving full cumulative curve...")
    cumulative_file = "output/src/Q_cumulative_A3.dat"
    
    output_array = np.column_stack([z, a, F_z])
    header = (
        "# z   a   F_z\n"
        "# Cumulative SRC injection fraction for A3 (gamma0 = 0.026)\n"
        "# F_z(z) ≈ fraction of total Q injected between z and 0 (today)"
    )
    
    np.savetxt(cumulative_file, output_array, header=header, fmt="%.10e")
    print(f"  Saved to: {cumulative_file}")
    
    # ------------------------------------------------------------------
    # 4) Load ΔH/H points and compare
    # ------------------------------------------------------------------
    print("\nLoading ΔH/H points for comparison...")
    deltaH_file = "output/src/deltaH_gamma_points.txt"
    deltaH_data = np.loadtxt(deltaH_file, comments="#")
    
    z_pts = deltaH_data[:, 0]
    dH_A3 = deltaH_data[:, 3]  # A3 column in percent
    
    print(f"Loaded {len(z_pts)} ΔH/H points")
    print(f"Redshift targets: {z_pts}")
    
    # Interpolate F_z at key redshifts
    F_interp = np.interp(z_pts, z, F_z)
    
    # ------------------------------------------------------------------
    # 5) Save comparison table
    # ------------------------------------------------------------------
    print("\nSaving comparison table...")
    comparison_file = "output/src/Q_cumulative_A3_points.txt"
    
    comparison_array = np.column_stack([z_pts, F_interp, dH_A3])
    header_comp = (
        "# z   F_z   dH_over_H_A3[%]\n"
        "# F_z = cumulative SRC injection fraction from z to today\n"
        "# dH_over_H_A3[%] from deltaH_gamma_points.txt"
    )
    
    np.savetxt(comparison_file, comparison_array, header=header_comp, fmt="%.10e")
    print(f"  Saved to: {comparison_file}")
    
    # ------------------------------------------------------------------
    # 6) Print summary
    # ------------------------------------------------------------------
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    print(f"Redshift range for cumulative table: z ∈ [{z.min():.6f}, {z.max():.6f}]")
    print()
    print("Comparison at key redshifts:")
    print("-" * 70)
    for z_t, F_val, dH_val in zip(z_pts, F_interp, dH_A3):
        print(f"z = {z_t:.2f}: F_z = {F_val:.4f},  ΔH/H_A3 = {dH_val:.4f}%")
    print("="*70)


if __name__ == "__main__":
    main()

