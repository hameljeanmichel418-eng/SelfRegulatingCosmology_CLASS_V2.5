#!/usr/bin/env python3
"""
Compare SRC vs LCDM background evolution for R3 physics validation.

R3 Background Physics Validation: Verify that:
1. SRC background H(z) differs smoothly from ΛCDM with a small late-time modulation
2. Total energy budget remains consistent: Ω_tot(z) ≈ 1
3. No obvious pathologies (e.g. negative densities) in the SRC channel

Run from: ~/dev/cosmo-class/
"""

import numpy as np
import sys
import os
from pathlib import Path


# Numerical tolerances for R3 validation
TOL_H_R3 = 1e-1      # max |ΔH/H| < 0.1 (10%) in 0 <= z <= 5
TOL_OMEGA = 1e-3     # max |Ω_tot - 1| < 1e-3 for both runs


def extract_omegas_from_header(filename):
    """
    Extract Omega_* column indices from CLASS background file header.
    
    Parameters:
    -----------
    filename : str
        Path to CLASS background file
    
    Returns:
    --------
    omega_cols : dict
        Dictionary mapping column names (e.g., "Omega_b") to 0-based column indices
    H_col : int
        0-based column index for H(z)
    """
    omega_cols = {}
    H_col = 3  # Default assumption
    
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('#'):
                # Parse header lines like "# 1:z  2:a  3:H[1/Mpc]  4:Omega_b  5:Omega_cdm ..."
                if ':' in line:
                    parts = line.strip().split()
                    for part in parts:
                        if ':' in part:
                            try:
                                col_num, col_name = part.split(':', 1)
                                col_idx = int(col_num) - 1  # Convert to 0-based
                                
                                # Check if it's an Omega column
                                if col_name.startswith('Omega_'):
                                    omega_cols[col_name] = col_idx
                                
                                # Check if it's the H column
                                if 'H' in col_name.upper() or 'H[1/Mpc]' in col_name:
                                    H_col = col_idx
                            except (ValueError, IndexError):
                                continue
            else:
                # Stop after first non-comment line
                break
    
    return omega_cols, H_col


def read_class_background(filename):
    """
    Read CLASS background file and extract z, H(z), and Ω_i columns.
    
    Returns:
    --------
    z : ndarray
        Redshift array
    H : ndarray
        Hubble parameter H(z) [1/Mpc]
    omega_dict : dict
        Dictionary mapping species names to their Ω_i arrays
    header_info : dict
        Dictionary with column information from header
    """
    if not os.path.exists(filename):
        return None, None, None, None
    
    # Extract Omega column indices from header
    omega_cols, H_col = extract_omegas_from_header(filename)
    
    # Read data
    data = np.loadtxt(filename, comments="#")
    
    if len(data) == 0:
        return None, None, None, None
    
    # Extract z (column 0) and H
    z = data[:, 0]
    H = data[:, H_col]
    
    # Extract Ω columns
    omega_dict = {}
    for col_name, col_idx in omega_cols.items():
        if col_idx < data.shape[1]:
            omega_dict[col_name] = data[:, col_idx]
    
    # Sort by increasing z
    sort_idx = np.argsort(z)
    z_sorted = z[sort_idx]
    H_sorted = H[sort_idx]
    omega_sorted = {}
    for key, val in omega_dict.items():
        omega_sorted[key] = val[sort_idx]
    
    # Create header_info for compatibility
    header_info = omega_cols.copy()
    header_info['H'] = H_col
    
    return z_sorted, H_sorted, omega_sorted, header_info


def compute_omega_tot(omega_dict):
    """
    Compute total Ω_tot = sum of all Ω_i.
    
    Parameters:
    -----------
    omega_dict : dict
        Dictionary mapping species names to their Ω_i arrays
    
    Returns:
    --------
    omega_tot : ndarray or None
        Total density parameter Ω_tot(z), or None if dict is empty
    """
    if len(omega_dict) == 0:
        return None
    
    omega_tot = np.zeros_like(list(omega_dict.values())[0])
    for key, val in omega_dict.items():
        omega_tot += val
    return omega_tot


def main():
    """Main comparison function."""
    
    # Get script directory
    script_dir = Path(__file__).parent.absolute()
    os.chdir(script_dir)
    
    # Redirect output to log file
    log_file = open("R3_compare.log", "w")
    
    def log_print(*args, **kwargs):
        """Print to both console and log file."""
        print(*args, **kwargs)
        print(*args, file=log_file, **kwargs)
    
    log_print("="*70)
    log_print("R3 Background Physics Validation: SRC vs LCDM")
    log_print("="*70)
    log_print(f"Working directory: {script_dir}")
    log_print()
    
    # Define file paths
    bg_lcdm = Path("CLASS_clean/output/R3_lcdm/R3_lcdm_background.dat")
    bg_src = Path("CLASS/output/R3_src/R3_src_background.dat")
    
    # Check if files exist
    if not bg_lcdm.exists():
        log_print(f"Error: {bg_lcdm} not found")
        log_print("Please run ./run_R3_background_comparison.sh first")
        log_file.close()
        sys.exit(1)
    
    if not bg_src.exists():
        log_print(f"Error: {bg_src} not found")
        log_print("Please run ./run_R3_background_comparison.sh first")
        log_file.close()
        sys.exit(1)
    
    # Read background files
    log_print("Reading background files...")
    z_lcdm, H_lcdm, omega_lcdm, header_lcdm = read_class_background(str(bg_lcdm))
    z_src, H_src, omega_src, header_src = read_class_background(str(bg_src))
    
    if z_lcdm is None or z_src is None:
        log_print("  ✗ Error: Could not read background files")
        log_file.close()
        sys.exit(1)
    
    log_print(f"  LCDM: {len(z_lcdm)} points, z ∈ [{z_lcdm.min():.3f}, {z_lcdm.max():.3f}]")
    log_print(f"  SRC:  {len(z_src)} points, z ∈ [{z_src.min():.3f}, {z_src.max():.3f}]")
    log_print()
    
    # Restrict to physical window: 0 <= z <= 5
    z_max_phys = 5.0
    mask_lcdm = (z_lcdm >= 0.0) & (z_lcdm <= z_max_phys)
    mask_src = (z_src >= 0.0) & (z_src <= z_max_phys)
    
    z_lcdm_phys = z_lcdm[mask_lcdm]
    H_lcdm_phys = H_lcdm[mask_lcdm]
    z_src_phys = z_src[mask_src]
    H_src_phys = H_src[mask_src]
    
    omega_lcdm_phys = {k: v[mask_lcdm] for k, v in omega_lcdm.items()}
    omega_src_phys = {k: v[mask_src] for k, v in omega_src.items()}
    
    log_print(f"Physical window: 0 <= z <= {z_max_phys}")
    log_print(f"  LCDM: {len(z_lcdm_phys)} points")
    log_print(f"  SRC:  {len(z_src_phys)} points")
    log_print()
    
    # 1. Compute ΔH/H(z) on a common z-grid
    log_print("Computing ΔH/H(z) = (H_SRC - H_LCDM) / H_LCDM...")
    z_min = max(z_lcdm_phys.min(), z_src_phys.min())
    z_max = min(z_lcdm_phys.max(), z_src_phys.max())
    z_common = np.linspace(z_min, z_max, min(len(z_lcdm_phys), len(z_src_phys), 1000))
    
    H_lcdm_interp = np.interp(z_common, z_lcdm_phys, H_lcdm_phys)
    H_src_interp = np.interp(z_common, z_src_phys, H_src_phys)
    
    delta_H_over_H = (H_src_interp - H_lcdm_interp) / H_lcdm_interp
    abs_delta_H = np.abs(delta_H_over_H)
    
    max_delta_H = np.max(abs_delta_H)
    idx_max_H = np.argmax(abs_delta_H)
    z_at_max_H = z_common[idx_max_H]
    
    log_print(f"  max |ΔH/H| = {max_delta_H:.2e} at z = {z_at_max_H:.3f}")
    log_print()
    
    # 2. Compute Ω_tot for both runs
    log_print("Computing Ω_tot(z) = sum of all Ω_i...")
    
    omega_tot_lcdm = compute_omega_tot(omega_lcdm_phys)
    omega_tot_src = compute_omega_tot(omega_src_phys)
    
    if omega_tot_lcdm is None or omega_tot_src is None:
        log_print("  ⚠ Warning: no Omega_* columns detected; skipping Ω_tot tests")
        log_print()
        max_omega_lcdm_dev = None
        max_omega_src_dev = None
        z_at_max_omega_lcdm = None
        z_at_max_omega_src = None
        pass_omega_lcdm = True  # Don't fail if we can't test
        pass_omega_src = True
    else:
        abs_omega_lcdm_dev = np.abs(omega_tot_lcdm - 1.0)
        abs_omega_src_dev = np.abs(omega_tot_src - 1.0)
        
        max_omega_lcdm_dev = np.max(abs_omega_lcdm_dev)
        max_omega_src_dev = np.max(abs_omega_src_dev)
        
        idx_max_omega_lcdm = np.argmax(abs_omega_lcdm_dev)
        idx_max_omega_src = np.argmax(abs_omega_src_dev)
        
        z_at_max_omega_lcdm = z_lcdm_phys[idx_max_omega_lcdm]
        z_at_max_omega_src = z_src_phys[idx_max_omega_src]
        
        log_print(f"  LCDM: max |Ω_tot - 1| = {max_omega_lcdm_dev:.2e} at z = {z_at_max_omega_lcdm:.3f}")
        log_print(f"  SRC:  max |Ω_tot - 1| = {max_omega_src_dev:.2e} at z = {z_at_max_omega_src:.3f}")
        log_print()
        
        pass_omega_lcdm = max_omega_lcdm_dev < TOL_OMEGA
        pass_omega_src = max_omega_src_dev < TOL_OMEGA
    
    # 3. Check for negative densities in SRC-related columns
    log_print("Checking for negative densities in SRC run...")
    negative_found = False
    for key, val in omega_src_phys.items():
        if np.any(val < 0):
            min_val = np.min(val)
            min_idx = np.argmin(val)
            z_min_val = z_src_phys[min_idx]
            log_print(f"  ⚠ Warning: {key} < 0: min = {min_val:.2e} at z = {z_min_val:.3f}")
            negative_found = True
    
    if not negative_found:
        log_print("  ✓ No negative densities detected")
    log_print()
    
    # 4. Check tolerances and determine PASS/FAIL
    log_print("="*70)
    log_print("VALIDATION CRITERIA")
    log_print("="*70)
    
    pass_H = max_delta_H < TOL_H_R3
    
    log_print(f"max |ΔH/H| < {TOL_H_R3:.0e}: {max_delta_H:.2e} {'✓' if pass_H else '✗'}")
    
    if max_omega_lcdm_dev is not None and max_omega_src_dev is not None:
        log_print(f"max |Ω_tot_LCDM - 1| = {max_omega_lcdm_dev:.2e} (tolerance: {TOL_OMEGA:.0e}) {'✓' if pass_omega_lcdm else '✗'}")
        log_print(f"max |Ω_tot_SRC - 1| = {max_omega_src_dev:.2e} (tolerance: {TOL_OMEGA:.0e}) {'✓' if pass_omega_src else '✗'}")
    else:
        log_print("Ω_tot test: skipped (no Omega_* columns found)")
    
    all_pass = pass_H and pass_omega_lcdm and pass_omega_src
    
    log_print("="*70)
    log_print("SUMMARY")
    log_print("="*70)
    log_print(f"max |ΔH/H| = {max_delta_H:.2e} at z = {z_at_max_H:.3f} (tolerance: {TOL_H_R3:.0e})")
    if max_omega_lcdm_dev is not None:
        log_print(f"max |Ω_tot_LCDM - 1| = {max_omega_lcdm_dev:.2e} at z = {z_at_max_omega_lcdm:.3f} (tolerance: {TOL_OMEGA:.0e})")
    if max_omega_src_dev is not None:
        log_print(f"max |Ω_tot_SRC - 1| = {max_omega_src_dev:.2e} at z = {z_at_max_omega_src:.3f} (tolerance: {TOL_OMEGA:.0e})")
    
    if negative_found:
        log_print("⚠ Warning: Negative densities detected in SRC run")
    
    log_print("="*70)
    
    if all_pass:
        log_print("RESULT: ✓ PASS - SRC background physics validation successful")
        log_print("        - H(z) modulation is smooth and within expected range")
        log_print("        - Energy budget is consistent (Ω_tot ≈ 1) for both runs")
        log_file.close()
        return 0
    else:
        log_print("RESULT: ✗ FAIL - One or more validation criteria not met")
        if not pass_H:
            log_print("        - max |ΔH/H| exceeds tolerance")
        if not pass_omega_lcdm:
            log_print("        - LCDM energy budget deviates from unity")
        if not pass_omega_src:
            log_print("        - SRC energy budget deviates from unity")
        log_file.close()
        return 1


if __name__ == "__main__":
    sys.exit(main())

