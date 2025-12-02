#!/usr/bin/env python3
"""
Compute total equation of state w_tot(z) for SRC gamma0 scan.

Reconstructs w_tot(z) from H(z) using finite differences on ln(H) vs ln(a).
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


def compute_wtot(z, H):
    """
    Compute total equation of state w_tot(z) from H(z).
    
    Uses finite differences on ln(H) vs ln(a) to compute d(ln H)/d(ln a),
    then w_tot = -1 - (2/3) * d(ln H)/d(ln a).
    
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
    w_tot : ndarray
        Total equation of state w_tot(z)
    """
    # Sort by increasing z
    sort_indices = np.argsort(z)
    z_sorted = z[sort_indices]
    H_sorted = H[sort_indices]
    
    # Compute scale factor and logarithms
    a = 1.0 / (1.0 + z_sorted)
    ln_a = np.log(a)
    ln_H = np.log(H_sorted)
    
    # Compute d(ln H)/d(ln a) using finite differences
    n = len(z_sorted)
    dlnH_dlna = np.zeros(n)
    
    # Central differences for interior points
    for i in range(1, n - 1):
        dlnH_dlna[i] = (ln_H[i+1] - ln_H[i-1]) / (ln_a[i+1] - ln_a[i-1])
    
    # Boundary conditions: copy neighbor value
    dlnH_dlna[0] = dlnH_dlna[1]
    dlnH_dlna[-1] = dlnH_dlna[-2]
    
    # Compute w_tot = -1 - (2/3) * d(ln H)/d(ln a)
    w_tot = -1.0 - (2.0 / 3.0) * dlnH_dlna
    
    return z_sorted, w_tot


def main():
    """Main computation function."""
    
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
    print("Loading reference run A0 (gamma0=0.000)...")
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
        
        # Compute w_tot
        z_sorted, w_tot = compute_wtot(z, H)
        
        # Store results
        results[tag] = {
            'z': z_sorted,
            'w_tot': w_tot,
            'gamma0': gamma0
        }
        
        # Save individual curve
        output_file = f"output/src/wtot_{tag}.dat"
        header = (
            f"# z   w_tot(z)\n"
            f"# tag = {tag}, gamma0 = {gamma0:.3f}"
        )
        
        output_array = np.column_stack([z_sorted, w_tot])
        np.savetxt(output_file, output_array, header=header, fmt="%.10e")
        print(f"  Saved to: {output_file}")
    
    # Build summary table at key redshifts
    print("\n" + "="*70)
    print("Building summary table at key redshifts...")
    print("="*70)
    
    z_targets = np.array([0.0, 0.25, 0.5, 0.75, 1.0])
    
    # Interpolate w_tot for A1..A5 at target redshifts
    w_A1_t = np.interp(z_targets, results["A1"]["z"], results["A1"]["w_tot"])
    w_A2_t = np.interp(z_targets, results["A2"]["z"], results["A2"]["w_tot"])
    w_A3_t = np.interp(z_targets, results["A3"]["z"], results["A3"]["w_tot"])
    w_A4_t = np.interp(z_targets, results["A4"]["z"], results["A4"]["w_tot"])
    w_A5_t = np.interp(z_targets, results["A5"]["z"], results["A5"]["w_tot"])
    
    # Build output array
    summary_array = np.column_stack([
        z_targets,
        w_A1_t,
        w_A2_t,
        w_A3_t,
        w_A4_t,
        w_A5_t,
    ])
    
    # Save summary
    summary_file = "output/src/wtot_gamma_points.txt"
    header_summary = (
        "# z   w_tot_A1   w_tot_A2   w_tot_A3   w_tot_A4   w_tot_A5\n"
        "# gamma0 values: 0.015 0.020 0.026 0.030 0.035"
    )
    
    np.savetxt(summary_file, summary_array, header=header_summary, fmt="%.10e")
    print(f"Summary saved to: {summary_file}")
    
    # Print summary
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    
    print(f"Redshift targets: {z_targets}")
    print()
    
    gamma_values = [0.015, 0.020, 0.026, 0.030, 0.035]
    run_tags = ["A1", "A2", "A3", "A4", "A5"]
    w_arrays = [w_A1_t, w_A2_t, w_A3_t, w_A4_t, w_A5_t]
    
    for tag, gamma0, w_t in zip(run_tags, gamma_values, w_arrays):
        w_str = ", ".join([f"{val:+.6f}" for val in w_t])
        print(f"{tag} (gamma0={gamma0:.3f}): w_tot(z) = [{w_str}]")
    
    print("="*70)


if __name__ == "__main__":
    main()

