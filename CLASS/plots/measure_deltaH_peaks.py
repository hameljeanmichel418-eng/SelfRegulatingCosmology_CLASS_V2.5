#!/usr/bin/env python3
"""
Measure peak ΔH/H values for SRC gamma0 scan.

For each run A1..A5, finds the maximum ΔH/H[%] in the range 0 <= z <= 1
and records the corresponding redshift z_peak.
"""

import numpy as np


def main():
    """Main peak measurement function."""
    
    # Load the zoomed deltaH data
    print("Loading deltaH_gamma_scan_z0_1.dat...")
    data = np.loadtxt("output/src/deltaH_gamma_scan_z0_1.dat", comments="#")
    
    # Extract columns
    z = data[:, 0]
    dA1 = data[:, 1]
    dA2 = data[:, 2]
    dA3 = data[:, 3]
    dA4 = data[:, 4]
    dA5 = data[:, 5]
    
    print(f"Loaded {len(z)} data points in range z = [{z.min():.6f}, {z.max():.6f}]")
    
    # Define runs and gamma0 values
    runs = [
        ("A1", 0.015),
        ("A2", 0.020),
        ("A3", 0.026),
        ("A4", 0.030),
        ("A5", 0.035),
    ]
    
    # Store deltaH arrays
    deltaH_arrays = [dA1, dA2, dA3, dA4, dA5]
    
    # Find peaks for each run
    results = []
    
    print("\nFinding peaks...")
    for (tag, gamma0), deltaH in zip(runs, deltaH_arrays):
        # Find index of maximum ΔH/H
        peak_idx = np.argmax(deltaH)
        
        # Record z_peak and deltaH_peak
        z_peak = z[peak_idx]
        deltaH_peak = deltaH[peak_idx]
        
        results.append((tag, gamma0, z_peak, deltaH_peak))
        
        print(f"  {tag} (gamma0={gamma0:.3f}): z_peak = {z_peak:.6f}, "
              f"ΔH/H_max = {deltaH_peak:.6f}%")
    
    # Build output array
    # Columns: tag, gamma0, z_peak, deltaH_peak
    # We'll save as text with formatted strings for tag
    output_data = []
    for tag, gamma0, z_peak, deltaH_peak in results:
        output_data.append([gamma0, z_peak, deltaH_peak])
    
    output_array = np.array(output_data)
    
    # Save to file
    output_file = "output/src/deltaH_gamma_peaks.txt"
    
    # Create header
    header = "# tag gamma0 z_peak deltaH_over_H_max[%]"
    
    # Save with formatted output including tags
    with open(output_file, 'w') as f:
        f.write(header + "\n")
        for (tag, gamma0, z_peak, deltaH_peak), _ in zip(results, output_array):
            f.write(f"{tag:3s} {gamma0:8.3f} {z_peak:12.6e} {deltaH_peak:15.6f}\n")
    
    print(f"\nResults saved to: {output_file}")
    
    # Print summary
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    
    for tag, gamma0, z_peak, deltaH_peak in results:
        print(f"{tag} (gamma0={gamma0:.3f}): z_peak = {z_peak:.6f}, "
              f"ΔH/H_max = {deltaH_peak:.6f}%")
    
    print("="*70)


if __name__ == "__main__":
    main()

