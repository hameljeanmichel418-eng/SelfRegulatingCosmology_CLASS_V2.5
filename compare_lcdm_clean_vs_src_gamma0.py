#!/usr/bin/env python3
"""
Compare CLASS_clean (LCDM) vs CLASS (SRC with gamma0=0) outputs.

R1 Regression Test: Verify that SRC-modified CLASS reduces to standard LCDM
when src_gamma0 = 0, and that it agrees with CLASS_clean at the 10^-3 to 10^-4 level.

Run from: ~/dev/cosmo-class/
"""

import numpy as np
import sys
import os
from pathlib import Path


# Numerical tolerances for R1 validation
TOL_H = 1e-4      # max ΔH/H
TOL_PK = 1e-3     # max ΔP/P
TOL_CL = 1e-3     # max ΔC_ℓ/C_ℓ


def read_class_background(filename):
    """
    Read CLASS background file and extract H(z).
    
    Returns:
    --------
    z : ndarray
        Redshift array
    H : ndarray
        Hubble parameter H(z) [1/Mpc]
    """
    if not os.path.exists(filename):
        return None, None
    
    data = np.loadtxt(filename, comments="#")
    
    # Column indices: 0=z, 3=H [1/Mpc]
    z = data[:, 0]
    H = data[:, 3]
    
    # Sort by increasing z
    sort_idx = np.argsort(z)
    z_sorted = z[sort_idx]
    H_sorted = H[sort_idx]
    
    return z_sorted, H_sorted


def read_class_pk(filename):
    """
    Read CLASS power spectrum file (single redshift per file).
    
    Parameters:
    -----------
    filename : str
        Path to CLASS pk file (e.g., lcdm_z1_pk.dat)
    
    Returns:
    --------
    k : ndarray
        Wavenumber array [h/Mpc]
    Pk : ndarray
        Power spectrum P(k) [(Mpc/h)^3]
    """
    if not os.path.exists(filename):
        return None, None
    
    # CLASS pk files have format:
    # Header lines starting with '#'
    # Data: k (h/Mpc) in column 0, P (Mpc/h)^3 in column 1
    data = np.loadtxt(filename, comments="#")
    
    k = data[:, 0]
    Pk = data[:, 1]
    
    return k, Pk


def read_class_cl(filename, spectrum='TT'):
    """
    Read CLASS CMB power spectrum file.
    
    Parameters:
    -----------
    filename : str
        Path to CLASS cl file
    spectrum : str
        Spectrum type: 'TT', 'EE', 'TE', etc.
    
    Returns:
    --------
    ell : ndarray
        Multipole array
    Cl : ndarray
        C_ℓ values
    """
    if not os.path.exists(filename):
        return None, None
    
    # CLASS cl files have format:
    # Header lines starting with '#'
    # Data line with column names: "1:l  2:TT  3:EE  4:TE ..."
    # Data: ell in column 0, Cl_TT in column 1, Cl_EE in column 2, etc.
    
    with open(filename, 'r') as f:
        # Find the header line with column names
        spectrum_idx = None
        for line in f:
            if line.startswith('#') and ':' in line:
                # Parse column names
                parts = line.strip().split()
                for i, part in enumerate(parts):
                    if spectrum.upper() in part.upper() and ':' in part:
                        # Extract column number (format: "2:TT")
                        try:
                            col_num = int(part.split(':')[0])
                            spectrum_idx = col_num - 1  # Convert to 0-based index
                            break
                        except (ValueError, IndexError):
                            continue
                if spectrum_idx is not None:
                    break
        
        if spectrum_idx is None:
            return None, None
        
        # Read data (skip header lines)
        f.seek(0)  # Reset to beginning
        data = np.loadtxt(f, comments="#")
        ell = data[:, 0]
        Cl = data[:, spectrum_idx]
    
    return ell, Cl


def compute_relative_diff(x_ref, y_ref, x_test, y_test):
    """
    Compute relative difference between two functions on potentially different grids.
    
    Parameters:
    -----------
    x_ref, y_ref : ndarray
        Reference function (x_ref, y_ref)
    x_test, y_test : ndarray
        Test function (x_test, y_test)
    
    Returns:
    --------
    max_rel_diff : float
        Maximum relative difference |y_test - y_ref| / |y_ref|
    """
    # Find common range
    x_min = max(x_ref.min(), x_test.min())
    x_max = min(x_ref.max(), x_test.max())
    
    if x_min >= x_max:
        return np.inf
    
    # Create common grid
    x_common = np.linspace(x_min, x_max, min(len(x_ref), len(x_test), 1000))
    
    # Interpolate both functions
    y_ref_interp = np.interp(x_common, x_ref, y_ref)
    y_test_interp = np.interp(x_common, x_test, y_test)
    
    # Compute relative difference
    mask = np.abs(y_ref_interp) > 0
    if not np.any(mask):
        return np.inf
    
    rel_diff = np.abs((y_test_interp[mask] - y_ref_interp[mask]) / y_ref_interp[mask])
    max_rel_diff = np.max(rel_diff)
    
    return max_rel_diff


def main():
    """Main comparison function."""
    
    # Get script directory
    script_dir = Path(__file__).parent.absolute()
    os.chdir(script_dir)
    
    print("="*70)
    print("R1 Regression Test: CLASS_clean vs CLASS (SRC gamma0=0)")
    print("="*70)
    print(f"Working directory: {script_dir}")
    print()
    
    # Define file paths
    base_clean = Path("CLASS_clean/output/lcdm")
    base_src = Path("CLASS/output/src_gamma0")
    
    # Check if directories exist
    if not base_clean.exists():
        print(f"Error: {base_clean} not found")
        print("Please run CLASS_clean/class src_science_LCDM.ini first")
        sys.exit(1)
    
    if not base_src.exists():
        print(f"Error: {base_src} not found")
        print("Please run CLASS/class src_science_SRC_gamma0.ini first")
        sys.exit(1)
    
    # File paths
    # Note: CLASS names files as z1_pk.dat for first z_pk value, z2_pk.dat for second, etc.
    # With z_pk = 0, 0.5, 1.0, 2.0: z1=0, z2=0.5, z3=1.0, z4=2.0
    bg_clean = base_clean / "lcdm_background.dat"
    bg_src = base_src / "src_gamma0_background.dat"
    
    pk_clean_z0 = base_clean / "lcdm_z1_pk.dat"  # z1 corresponds to first z_pk value (0.0)
    pk_src_z0 = base_src / "src_gamma0_z1_pk.dat"
    
    pk_clean_z1 = base_clean / "lcdm_z3_pk.dat"  # z3 corresponds to third z_pk value (1.0)
    pk_src_z1 = base_src / "src_gamma0_z3_pk.dat"
    
    cl_clean = base_clean / "lcdm_cl.dat"
    cl_src = base_src / "src_gamma0_cl.dat"
    
    results = {}
    
    # 1. Compare H(z)
    print("Comparing H(z) from background files...")
    z_clean, H_clean = read_class_background(str(bg_clean))
    z_src, H_src = read_class_background(str(bg_src))
    
    if z_clean is None or z_src is None:
        print("  ✗ Error: Could not read background files")
        results['H'] = (False, np.inf, "File read error")
    else:
        max_delta_H = compute_relative_diff(z_clean, H_clean, z_src, H_src)
        pass_H = max_delta_H < TOL_H
        results['H'] = (pass_H, max_delta_H, f"max ΔH/H = {max_delta_H:.2e}")
        status = "✓" if pass_H else "✗"
        print(f"  {status} {results['H'][2]} (tolerance: {TOL_H:.0e})")
    
    # 2. Compare P(k) at z=0 (from z1_pk.dat files)
    print("\nComparing P(k) at z=0...")
    k_clean_z0, Pk_clean_z0 = read_class_pk(str(pk_clean_z0))
    k_src_z0, Pk_src_z0 = read_class_pk(str(pk_src_z0))
    
    if k_clean_z0 is None or k_src_z0 is None:
        print("  ✗ Error: Could not read P(k) files at z=0")
        results['Pk_z0'] = (False, np.inf, "File read error")
    else:
        max_delta_Pk_z0 = compute_relative_diff(k_clean_z0, Pk_clean_z0, k_src_z0, Pk_src_z0)
        pass_Pk_z0 = max_delta_Pk_z0 < TOL_PK
        results['Pk_z0'] = (pass_Pk_z0, max_delta_Pk_z0, f"max ΔP/P = {max_delta_Pk_z0:.2e}")
        status = "✓" if pass_Pk_z0 else "✗"
        print(f"  {status} {results['Pk_z0'][2]} (tolerance: {TOL_PK:.0e})")
    
    # 3. Compare P(k) at z=1 (from z3_pk.dat files)
    print("\nComparing P(k) at z=1...")
    k_clean_z1, Pk_clean_z1 = read_class_pk(str(pk_clean_z1))
    k_src_z1, Pk_src_z1 = read_class_pk(str(pk_src_z1))
    
    if k_clean_z1 is None or k_src_z1 is None:
        print("  ✗ Error: Could not read P(k) files at z=1")
        results['Pk_z1'] = (False, np.inf, "File read error")
    else:
        max_delta_Pk_z1 = compute_relative_diff(k_clean_z1, Pk_clean_z1, k_src_z1, Pk_src_z1)
        pass_Pk_z1 = max_delta_Pk_z1 < TOL_PK
        results['Pk_z1'] = (pass_Pk_z1, max_delta_Pk_z1, f"max ΔP/P = {max_delta_Pk_z1:.2e}")
        status = "✓" if pass_Pk_z1 else "✗"
        print(f"  {status} {results['Pk_z1'][2]} (tolerance: {TOL_PK:.0e})")
    
    # 4. Compare C_ℓ^TT
    print("\nComparing C_ℓ^TT from CMB files...")
    ell_clean, Cl_clean = read_class_cl(str(cl_clean), 'TT')
    ell_src, Cl_src = read_class_cl(str(cl_src), 'TT')
    
    if ell_clean is None or ell_src is None:
        print("  ✗ Error: Could not read CMB files")
        results['Cl'] = (False, np.inf, "File read error")
    else:
        max_delta_Cl = compute_relative_diff(ell_clean, Cl_clean, ell_src, Cl_src)
        pass_Cl = max_delta_Cl < TOL_CL
        results['Cl'] = (pass_Cl, max_delta_Cl, f"max ΔC_ℓ/C_ℓ = {max_delta_Cl:.2e}")
        status = "✓" if pass_Cl else "✗"
        print(f"  {status} {results['Cl'][2]} (tolerance: {TOL_CL:.0e})")
    
    # Summary
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    
    all_pass = all(r[0] for r in results.values())
    
    print(f"max ΔH/H        = {results['H'][1]:.2e} (tolerance: {TOL_H:.0e}) {'✓' if results['H'][0] else '✗'}")
    print(f"max ΔP/P (z=0)  = {results['Pk_z0'][1]:.2e} (tolerance: {TOL_PK:.0e}) {'✓' if results['Pk_z0'][0] else '✗'}")
    print(f"max ΔP/P (z=1)  = {results['Pk_z1'][1]:.2e} (tolerance: {TOL_PK:.0e}) {'✓' if results['Pk_z1'][0] else '✗'}")
    print(f"max ΔC_ℓ/C_ℓ    = {results['Cl'][1]:.2e} (tolerance: {TOL_CL:.0e}) {'✓' if results['Cl'][0] else '✗'}")
    
    print("="*70)
    
    if all_pass:
        print("RESULT: ✓ PASS - SRC implementation is numerically indistinguishable")
        print("        from standard CLASS when src_gamma0 = 0")
        return 0
    else:
        print("RESULT: ✗ FAIL - Differences exceed tolerances")
        return 1


if __name__ == "__main__":
    sys.exit(main())

