#!/usr/bin/env python3
"""
Compute dH/dz for SRC gamma0 scan.

Computes the derivative dH/dz from H(z) using finite differences.
"""

import numpy as np


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
    
    # Load data, skipping comment lines
    # Column indices: 0=z, 3=H [1/Mpc]
    data = np.loadtxt(filename, comments="#")
    
    z_array = data[:, 0]
    H_array = data[:, 3]
    
    return z_array, H_array


def compute_dH_dz(z, H):
    """
    Compute dH/dz from H(z) using finite differences.
    
    Parameters:
    -----------
    z : ndarray
        1D array of redshift values
    H : ndarray
        1D array of Hubble parameter H(z) [1/Mpc]
    
    Returns:
    --------
    z_sorted : ndarray
        Sorted redshift array (increasing)
    dH_dz : ndarray
        Derivative dH/dz [1/Mpc per unit z]
    """
    # Sort by increasing z
    sort_indices = np.argsort(z)
    z_sorted = z[sort_indices]
    H_sorted = H[sort_indices]
    
    # Compute dH/dz using finite differences
    n = len(z_sorted)
    dH_dz = np.zeros(n)
    
    # Central differences for interior points
    for i in range(1, n - 1):
        dH_dz[i] = (H_sorted[i+1] - H_sorted[i-1]) / (z_sorted[i+1] - z_sorted[i-1])
    
    # Boundary conditions: copy neighbor value
    dH_dz[0] = dH_dz[1]
    dH_dz[-1] = dH_dz[-2]
    
    return z_sorted, dH_dz


def main():
    """Main computation function."""
    
    print("Computing dH/dz for A0..A5 from CLASS background outputs")
    print("="*70)
    
    # Define runs and gamma0 values
    runs = [
        ("A0", 0.000),
        ("A1", 0.015),
        ("A2", 0.020),
        ("A3", 0.026),
        ("A4", 0.030),
        ("A5", 0.035),
    ]
    
    # Load reference (A0) to get shape
    print("\nLoading reference run A0 (gamma0=0.000)...")
    z_ref, H_ref = load_background("A0")
    n_ref = len(z_ref)
    print(f"Reference has {n_ref} data points")
    
    # Dictionary to store results
    results = {}
    
    # Process each run
    for tag, gamma0 in runs:
        print(f"\nProcessing {tag} (gamma0={gamma0:.3f})...")
        
        # Load data
        z, H = load_background(tag)
        
        # Check shape consistency
        if len(z) != n_ref:
            raise ValueError(
                f"Shape mismatch: {tag} has {len(z)} points, "
                f"reference has {n_ref} points"
            )
        
        # Compute dH/dz
        z_sorted, dH_dz = compute_dH_dz(z, H)
        
        # Store results
        results[tag] = {
            'z': z_sorted,
            'dH_dz': dH_dz,
            'gamma0': gamma0
        }
        
        # Save individual curve
        output_file = f"output/src/dHdz_{tag}.dat"
        header = (
            f"# z   dH_dz(z) [1/Mpc per unit z]\n"
            f"# tag = {tag}, gamma0 = {gamma0:.3f}"
        )
        
        output_array = np.column_stack([z_sorted, dH_dz])
        np.savetxt(output_file, output_array, header=header, fmt="%.10e")
        print(f"  Saved to: {output_file}")
    
    # Build summary table at key redshifts
    print("\n" + "="*70)
    print("Building summary table at key redshifts...")
    print("="*70)
    
    z_targets = np.array([0.0, 0.25, 0.5, 0.75, 1.0])
    
    # Interpolate dH/dz for A0, A2, A3, A4 at target redshifts
    dH_A0_t = np.interp(z_targets, results["A0"]["z"], results["A0"]["dH_dz"])
    dH_A2_t = np.interp(z_targets, results["A2"]["z"], results["A2"]["dH_dz"])
    dH_A3_t = np.interp(z_targets, results["A3"]["z"], results["A3"]["dH_dz"])
    dH_A4_t = np.interp(z_targets, results["A4"]["z"], results["A4"]["dH_dz"])
    
    # Build output array
    summary_array = np.column_stack([
        z_targets,
        dH_A0_t,
        dH_A2_t,
        dH_A3_t,
        dH_A4_t,
    ])
    
    # Save summary
    summary_file = "output/src/dHdz_gamma_points.txt"
    header_summary = (
        "# z   dH_dz_A0   dH_dz_A2   dH_dz_A3   dH_dz_A4\n"
        "# gamma0 values: A0=0.000 (LCDM), A2=0.020, A3=0.026, A4=0.030"
    )
    
    np.savetxt(summary_file, summary_array, header=header_summary, fmt="%.10e")
    print(f"Summary saved to: {summary_file}")
    
    # Print summary
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    
    print(f"Redshift targets: {z_targets}")
    print()
    
    # Compute relative differences for A2, A3, A4
    gamma_values = [0.020, 0.026, 0.030]
    run_tags = ["A2", "A3", "A4"]
    dH_arrays = [dH_A2_t, dH_A3_t, dH_A4_t]
    
    for tag, gamma0, dH_k in zip(run_tags, gamma_values, dH_arrays):
        # Compute relative difference in percent
        delta = (dH_k - dH_A0_t) / np.abs(dH_A0_t) * 100.0
        delta_str = ", ".join([f"{val:+.6f}" for val in delta])
        print(f"{tag} (gamma0={gamma0:.3f}): Δ(dH/dz)/|dH/dz| [%] = [{delta_str}]")
    
    print("="*70)


if __name__ == "__main__":
    main()

