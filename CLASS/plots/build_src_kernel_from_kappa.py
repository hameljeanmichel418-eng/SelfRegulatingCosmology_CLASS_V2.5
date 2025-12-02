#!/usr/bin/env python3
"""
Build normalized SRC kernel from kappa_eff(z) for fiducial run A3.

Extracts the shape of the energy-transfer window from the fiducial
SRC run (A3, gamma0 = 0.026) for use in Volume 2.5 analysis.
"""

import numpy as np
import os


def main():
    """Main function to build SRC kernel from kappa data."""
    
    print("Building SRC kernel from kappa_eff(z) for fiducial run A3")
    print("="*70)
    
    # Load kappa_eff(z) for A3
    kappa_file = "output/src/kappa_A3.dat"
    
    if not os.path.exists(kappa_file):
        raise FileNotFoundError(f"Could not find file: {kappa_file}")
    
    print(f"\nLoading kappa data from: {kappa_file}")
    data = np.loadtxt(kappa_file, comments="#")
    
    z = data[:, 0]
    k = data[:, 1]  # kappa_eff(z)
    
    print(f"Loaded {len(z)} data points")
    print(f"Redshift range: z ∈ [{z.min():.6f}, {z.max():.6f}]")
    
    # Sort by increasing z
    sort_indices = np.argsort(z)
    z_sorted = z[sort_indices]
    k_sorted = k[sort_indices]
    
    # Optionally restrict to physical range 0 <= z <= 1
    mask_phys = (z_sorted >= 0.0) & (z_sorted <= 1.0)
    z_phys = z_sorted[mask_phys]
    k_phys = k_sorted[mask_phys]
    
    print(f"Physical range (0 <= z <= 1): {len(z_phys)} points")
    print(f"  z ∈ [{z_phys.min():.6f}, {z_phys.max():.6f}]")
    
    # Compute scale factor
    a_phys = 1.0 / (1.0 + z_phys)
    
    # Find peak location
    k_max = np.max(k_phys)
    idx_max = np.argmax(k_phys)
    z_peak = z_phys[idx_max]
    a_peak = a_phys[idx_max]
    
    print(f"\nPeak location:")
    print(f"  kappa_peak = {k_max:.6e}")
    print(f"  z_peak = {z_peak:.6f}")
    print(f"  a_peak = {a_peak:.6f}")
    
    # Build normalized kernel
    kernel_norm = k_phys / k_max
    
    print(f"\nNormalized kernel: min = {kernel_norm.min():.6f}, max = {kernel_norm.max():.6f}")
    
    # Characterize window width
    thr_half = 0.5  # FWHM threshold
    thr_ten = 0.1   # 10% support threshold
    
    mask_half = (kernel_norm >= thr_half)
    mask_ten = (kernel_norm >= thr_ten)
    
    # FWHM window
    if np.any(mask_half):
        z_min_half = z_phys[mask_half].min()
        z_max_half = z_phys[mask_half].max()
        a_min_half = 1.0 / (1.0 + z_max_half)
        a_max_half = 1.0 / (1.0 + z_min_half)
    else:
        z_min_half = z_max_half = np.nan
        a_min_half = a_max_half = np.nan
    
    # 10% support window
    if np.any(mask_ten):
        z_min_ten = z_phys[mask_ten].min()
        z_max_ten = z_phys[mask_ten].max()
        a_min_ten = 1.0 / (1.0 + z_max_ten)
        a_max_ten = 1.0 / (1.0 + z_min_ten)
    else:
        z_min_ten = z_max_ten = np.nan
        a_min_ten = a_max_ten = np.nan
    
    # Save normalized kernel table
    output_file = "output/src/src_kernel_A3_normalized.dat"
    
    output_array = np.column_stack([
        z_phys,
        a_phys,
        k_phys,
        kernel_norm,
    ])
    
    header = (
        "# z   a   kappa_eff   kernel_norm\n"
        "# Fiducial SRC kernel from run A3 (gamma0 = 0.026)\n"
        "# kappa_eff is dimensionless; kernel_norm = kappa_eff / max(kappa_eff) over 0<=z<=1"
    )
    
    np.savetxt(output_file, output_array, header=header, fmt="%.10e")
    print(f"\nSaved normalized kernel to: {output_file}")
    
    # Print detailed summary
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    
    print(f"\nPeak:")
    print(f"  z_peak = {z_peak:.6f}")
    print(f"  a_peak = {a_peak:.6f}")
    print(f"  kappa_peak = {k_max:.6e}")
    
    print(f"\nFWHM window (kernel_norm >= 0.5):")
    if not np.isnan(z_min_half):
        print(f"  z in [{z_min_half:.6f}, {z_max_half:.6f}]")
        print(f"  a in [{a_min_half:.6f}, {a_max_half:.6f}]")
    else:
        print("  No points above threshold")
    
    print(f"\n10% support window (kernel_norm >= 0.1):")
    if not np.isnan(z_min_ten):
        print(f"  z in [{z_min_ten:.6f}, {z_max_ten:.6f}]")
        print(f"  a in [{a_min_ten:.6f}, {a_max_ten:.6f}]")
    else:
        print("  No points above threshold")
    
    print("="*70)


if __name__ == "__main__":
    main()

