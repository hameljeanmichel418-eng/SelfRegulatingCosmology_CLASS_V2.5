#!/usr/bin/env python3
"""
Compute ISW proxy S(z) = d(1/H)/dz for SRC gamma0 scan.

S(z) = - dH/dz / H^2
"""

import numpy as np


def load_background_H(tag):
    """
    Load z, H(z) from output/src/{tag}__background.dat.
    
    Parameters:
    -----------
    tag : str
        Run tag (e.g., "A0", "A2", etc.)
    
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


def load_dH_dz(tag):
    """
    Load z, dH/dz from output/src/dHdz_{tag}.dat.
    
    Parameters:
    -----------
    tag : str
        Run tag (e.g., "A0", "A2", etc.)
    
    Returns:
    --------
    z_array : ndarray
        1D array of redshift values
    dH_dz_array : ndarray
        1D array of dH/dz [1/Mpc per unit z]
    """
    filename = f"output/src/dHdz_{tag}.dat"
    
    # Load data, skipping comment lines
    # Column indices: 0=z, 1=dH_dz
    data = np.loadtxt(filename, comments="#")
    
    z_array = data[:, 0]
    dH_dz_array = data[:, 1]
    
    return z_array, dH_dz_array


def build_isw_proxy(z_H, H, z_dH, dH_dz):
    """
    Build S(z) = - dH/dz / H^2 on a common z grid.
    
    Parameters:
    -----------
    z_H : ndarray
        Redshift array from H data
    H : ndarray
        Hubble parameter H(z) [1/Mpc]
    z_dH : ndarray
        Redshift array from dH/dz data
    dH_dz : ndarray
        dH/dz [1/Mpc per unit z]
    
    Returns:
    --------
    z_sorted : ndarray
        Sorted redshift array (increasing)
    S : ndarray
        ISW proxy S(z) = - dH/dz / H^2 [Mpc]
    """
    # Sort both arrays by increasing z
    sort_indices_H = np.argsort(z_H)
    z_H_sorted = z_H[sort_indices_H]
    H_sorted = H[sort_indices_H]
    
    sort_indices_dH = np.argsort(z_dH)
    z_dH_sorted = z_dH[sort_indices_dH]
    dH_dz_sorted = dH_dz[sort_indices_dH]
    
    # Check that z grids match
    if len(z_H_sorted) != len(z_dH_sorted):
        raise ValueError(
            f"Length mismatch: H has {len(z_H_sorted)} points, "
            f"dH/dz has {len(z_dH_sorted)} points"
        )
    
    if not np.allclose(z_H_sorted, z_dH_sorted, rtol=1e-10):
        raise ValueError("Redshift grids from H and dH/dz do not match")
    
    # Compute S(z) = - dH_dz / H^2
    S = -dH_dz_sorted / (H_sorted ** 2)
    
    return z_H_sorted, S


def main():
    """Main computation function."""
    
    print("Computing ISW proxy S(z) = d(1/H)/dz for A0, A2, A3, A4")
    print("="*70)
    
    # Define runs of interest
    runs = [
        ("A0", 0.000),  # LCDM
        ("A2", 0.020),
        ("A3", 0.026),
        ("A4", 0.030),
    ]
    
    # Dictionary to store results
    results = {}
    
    # Process each run
    for tag, gamma0 in runs:
        print(f"\nProcessing {tag} (gamma0={gamma0:.3f})...")
        
        # Load H and dH/dz
        z_H, H = load_background_H(tag)
        z_dH, dH_dz = load_dH_dz(tag)
        
        # Build ISW proxy (handles sorting and consistency checks)
        z_S, S = build_isw_proxy(z_H, H, z_dH, dH_dz)
        
        # Store results
        results[tag] = {
            'z': z_S,
            'S': S,
            'gamma0': gamma0
        }
        
        # Save individual curve
        output_file = f"output/src/iswproxy_{tag}.dat"
        header = (
            f"# z   S(z) = d(1/H)/dz [Mpc]\n"
            f"# tag = {tag}, gamma0 = {gamma0:.3f}"
        )
        
        output_array = np.column_stack([z_S, S])
        np.savetxt(output_file, output_array, header=header, fmt="%.10e")
        print(f"  Saved to: {output_file}")
    
    # Build summary table at key redshifts
    print("\n" + "="*70)
    print("Building summary table at key redshifts...")
    print("="*70)
    
    z_targets = np.array([0.0, 0.25, 0.5, 0.75, 1.0])
    
    # Interpolate S(z) for A0, A2, A3, A4 at target redshifts
    S_A0_t = np.interp(z_targets, results["A0"]["z"], results["A0"]["S"])
    S_A2_t = np.interp(z_targets, results["A2"]["z"], results["A2"]["S"])
    S_A3_t = np.interp(z_targets, results["A3"]["z"], results["A3"]["S"])
    S_A4_t = np.interp(z_targets, results["A4"]["z"], results["A4"]["S"])
    
    # Build output array
    summary_array = np.column_stack([
        z_targets,
        S_A0_t,
        S_A2_t,
        S_A3_t,
        S_A4_t,
    ])
    
    # Save summary
    summary_file = "output/src/iswproxy_gamma_points.txt"
    header_summary = (
        "# z   S_A0   S_A2   S_A3   S_A4\n"
        "# S(z) = d(1/H)/dz [Mpc]\n"
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
    S_arrays = [S_A2_t, S_A3_t, S_A4_t]
    
    for tag, gamma0, S_k in zip(run_tags, gamma_values, S_arrays):
        # Compute relative difference in percent
        delta = (S_k - S_A0_t) / np.abs(S_A0_t) * 100.0
        delta_str = ", ".join([f"{val:+.6f}" for val in delta])
        print(f"{tag} (gamma0={gamma0:.3f}): ΔS/|S|_LCDM [%] = [{delta_str}]")
    
    print("="*70)


if __name__ == "__main__":
    main()

