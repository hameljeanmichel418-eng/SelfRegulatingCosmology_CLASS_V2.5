#!/usr/bin/env python3
"""
Extract zoomed ΔH/H data and key points from gamma0 scan results.

1) Extracts the physical range 0 <= z <= 1
2) Interpolates ΔH/H at key redshifts: z = [0.0, 0.25, 0.5, 0.75, 1.0]
3) Saves both zoomed table and key points to disk
"""

import numpy as np


def main():
    """Main extraction and interpolation function."""
    
    # Load the full deltaH scan data
    print("Loading deltaH_gamma_scan.dat...")
    data = np.loadtxt("output/src/deltaH_gamma_scan.dat", comments="#")
    
    # Extract columns
    z = data[:, 0]
    dA1 = data[:, 1]
    dA2 = data[:, 2]
    dA3 = data[:, 3]
    dA4 = data[:, 4]
    dA5 = data[:, 5]
    
    print(f"Loaded {len(z)} data points")
    print(f"Redshift range in file: z_min = {z.min():.6e}, z_max = {z.max():.6e}")
    
    # Sort all arrays by z in ascending order
    # (currently sorted high-to-low)
    sort_indices = np.argsort(z)
    z_sorted = z[sort_indices]
    dA1_sorted = dA1[sort_indices]
    dA2_sorted = dA2[sort_indices]
    dA3_sorted = dA3[sort_indices]
    dA4_sorted = dA4[sort_indices]
    dA5_sorted = dA5[sort_indices]
    
    # Build mask for physical range 0 <= z <= 1
    mask = (z_sorted >= 0.0) & (z_sorted <= 1.0)
    
    # Extract physical range data
    z_phys = z_sorted[mask]
    dA1_phys = dA1_sorted[mask]
    dA2_phys = dA2_sorted[mask]
    dA3_phys = dA3_sorted[mask]
    dA4_phys = dA4_sorted[mask]
    dA5_phys = dA5_sorted[mask]
    
    print(f"\nExtracted {len(z_phys)} points in range 0 <= z <= 1")
    print(f"Physical range: z_min = {z_phys.min():.6f}, z_max = {z_phys.max():.6f}")
    
    # Save zoomed data to file
    zoomed_output = "output/src/deltaH_gamma_scan_z0_1.dat"
    
    # Build output array for zoomed data
    zoomed_array = np.column_stack([
        z_phys,
        dA1_phys,
        dA2_phys,
        dA3_phys,
        dA4_phys,
        dA5_phys,
    ])
    
    # Create header
    header_zoomed = (
        "# z   dH_over_H_A1[%]   dH_over_H_A2[%]   dH_over_H_A3[%]   "
        "dH_over_H_A4[%]   dH_over_H_A5[%]\n"
        "# gamma0 values: 0.015 0.020 0.026 0.030 0.035"
    )
    
    np.savetxt(zoomed_output, zoomed_array, header=header_zoomed, fmt="%.10e")
    print(f"Zoomed data saved to: {zoomed_output}")
    
    # Define target redshifts for interpolation
    z_targets = np.array([0.0, 0.25, 0.5, 0.75, 1.0])
    
    # Interpolate ΔH/H at target redshifts for each run
    print("\nInterpolating at target redshifts...")
    dA1_t = np.interp(z_targets, z_phys, dA1_phys)
    dA2_t = np.interp(z_targets, z_phys, dA2_phys)
    dA3_t = np.interp(z_targets, z_phys, dA3_phys)
    dA4_t = np.interp(z_targets, z_phys, dA4_phys)
    dA5_t = np.interp(z_targets, z_phys, dA5_phys)
    
    # Build output array for key points
    points_array = np.column_stack([
        z_targets,
        dA1_t,
        dA2_t,
        dA3_t,
        dA4_t,
        dA5_t,
    ])
    
    # Save key points to file
    points_output = "output/src/deltaH_gamma_points.txt"
    
    # Create header
    header_points = (
        "# z   dH_over_H_A1[%]   dH_over_H_A2[%]   dH_over_H_A3[%]   "
        "dH_over_H_A4[%]   dH_over_H_A5[%]\n"
        "# gamma0 values: 0.015 0.020 0.026 0.030 0.035"
    )
    
    np.savetxt(points_output, points_array, header=header_points, fmt="%.10e")
    print(f"Key points saved to: {points_output}")
    
    # Print summary
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    
    # Define gamma0 values for each run
    gamma_values = [0.015, 0.020, 0.026, 0.030, 0.035]
    run_tags = ["A1", "A2", "A3", "A4", "A5"]
    deltaH_arrays = [dA1_t, dA2_t, dA3_t, dA4_t, dA5_t]
    
    print(f"Redshift targets: {z_targets}")
    print()
    
    for tag, gamma0, deltaH in zip(run_tags, gamma_values, deltaH_arrays):
        # Format values as strings for clean output
        deltaH_str = ", ".join([f"{val:+.6f}" for val in deltaH])
        print(f"{tag} (gamma0={gamma0:.3f}): ΔH/H[%] = [{deltaH_str}]")
    
    print("="*70)


if __name__ == "__main__":
    main()

