#!/usr/bin/env python3
"""
Analyze ΔH/H(z) for SRC gamma0 scan.

Computes the fractional change in Hubble parameter H(z) relative to
the reference run A0 (gamma0=0.000) for runs A1 through A5.
"""

import numpy as np
import os


def load_background(tag):
    """
    Load background data from CLASS output file.
    
    Parameters:
    -----------
    tag : str
        Run tag (e.g., "A0", "A1", etc.)
    
    Returns:
    --------
    z_array : ndarray
        1D array of redshift values
    H_array : ndarray
        1D array of Hubble parameter H(z) [1/Mpc]
    """
    filename = f"output/src/{tag}__background.dat"
    
    if not os.path.exists(filename):
        raise FileNotFoundError(f"Background file not found: {filename}")
    
    # Load data, skipping comment lines
    # Column indices: 0=z, 3=H [1/Mpc]
    data = np.loadtxt(filename, comments="#")
    
    z_array = data[:, 0]
    H_array = data[:, 3]
    
    return z_array, H_array


def main():
    """Main analysis function."""
    
    # Define runs and gamma0 values
    runs = [
        ("A0", 0.000),
        ("A1", 0.015),
        ("A2", 0.020),
        ("A3", 0.026),
        ("A4", 0.030),
        ("A5", 0.035),
    ]
    
    # Load reference (A0)
    print("Loading reference run A0 (gamma0=0.000)...")
    z_ref, H_ref = load_background("A0")
    
    # Dictionary to store results
    results = {}
    
    # Process each run (A1 through A5)
    for tag, gamma0 in runs[1:]:  # Skip A0, already loaded
        print(f"Processing {tag} (gamma0={gamma0:.3f})...")
        
        # Load data
        z_k, H_k = load_background(tag)
        
        # Check array shapes match
        if z_k.shape != z_ref.shape:
            raise ValueError(
                f"Shape mismatch: {tag} has {z_k.shape[0]} points, "
                f"reference has {z_ref.shape[0]} points"
            )
        
        # Check redshift arrays match (within numerical precision)
        if not np.allclose(z_k, z_ref, rtol=1e-10):
            raise ValueError(
                f"Redshift arrays differ for {tag} compared to reference"
            )
        
        # Compute fractional change in percent
        dH_over_H = (H_k - H_ref) / H_ref * 100.0
        
        # Store results
        results[tag] = dH_over_H
    
    # Build output array
    # Columns: z_ref, dH_over_H_A1, dH_over_H_A2, ..., dH_over_H_A5
    output_array = np.column_stack([
        z_ref,
        results["A1"],
        results["A2"],
        results["A3"],
        results["A4"],
        results["A5"],
    ])
    
    # Save to file
    output_file = "output/src/deltaH_gamma_scan.dat"
    
    # Create header
    header = (
        "# z   dH_over_H_A1[%]   dH_over_H_A2[%]   dH_over_H_A3[%]   "
        "dH_over_H_A4[%]   dH_over_H_A5[%]\n"
        "# gamma0 values: 0.015 0.020 0.026 0.030 0.035"
    )
    
    np.savetxt(output_file, output_array, header=header, fmt="%.10e")
    
    print(f"\nResults saved to: {output_file}")
    
    # Print summary
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    print(f"Redshift range: z_min = {z_ref.min():.6e}, z_max = {z_ref.max():.6e}")
    print(f"Number of points: {len(z_ref)}")
    print("\nΔH/H [%] statistics:")
    print("-" * 60)
    
    for tag, gamma0 in runs[1:]:
        dH = results[tag]
        print(f"{tag} (gamma0={gamma0:.3f}): "
              f"min = {dH.min():.6f}%, max = {dH.max():.6f}%")
    
    print("="*60)


if __name__ == "__main__":
    main()

