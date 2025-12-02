#!/usr/bin/env python3
"""
Compare SRC vs LCDM perturbation spectra for R4 physics validation.

R4 Perturbation Spectra Validation: Verify that:
1. Matter power spectra P(k,z) differ smoothly from ΛCDM with reasonable amplitude
2. CMB angular power spectra C_ℓ^TT differ smoothly from ΛCDM
3. Differences are O(1-5%) in the physical window (k ≤ 1 h/Mpc, ℓ ≤ 2000)

Run from: ~/dev/cosmo-class/
"""

import numpy as np
import sys
import os
from pathlib import Path


# Numerical tolerances for R4 validation
TOL_Pk_phys = 0.02      # 2% tolerance on P(k) for k <= 1 h/Mpc (physical window)
TOL_Pk_full = 0.10      # 10% tolerance on P(k) for full range (diagnostic only)
TOL_Cl_phys = 0.02      # 2% tolerance on C_ℓ for ℓ <= 2000 (physical window)
TOL_Cl_full = 0.05      # 5% tolerance on C_ℓ for full range (diagnostic only)


def read_class_pk(filename):
    """
    Read CLASS power spectrum file (single redshift per file).
    
    Parameters:
    -----------
    filename : str
        Path to CLASS pk file
    
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


def read_class_cl_lensed(filename):
    """
    Read CLASS lensed CMB power spectrum file.
    
    Parameters:
    -----------
    filename : str
        Path to CLASS cl_lensed file
    
    Returns:
    --------
    ell : ndarray
        Multipole array
    Cl : ndarray
        C_ℓ^TT values (lensed)
    """
    if not os.path.exists(filename):
        return None, None
    
    # CLASS cl_lensed files have format:
    # Header lines starting with '#'
    # Data line with column names: "1:l  2:TT  3:EE  4:TE ..."
    # Data: ell in column 0, Cl_TT in column 1
    
    with open(filename, 'r') as f:
        # Find the header line with column names
        spectrum_idx = None
        for line in f:
            if line.startswith('#') and ':' in line and 'TT' in line.upper():
                # Parse column names
                parts = line.strip().split()
                for i, part in enumerate(parts):
                    if 'TT' in part.upper() and ':' in part:
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


def compute_relative_diff(x_ref, y_ref, x_test, y_test, return_position=False):
    """
    Compute relative difference between two functions on potentially different grids.
    
    Parameters:
    -----------
    x_ref, y_ref : ndarray
        Reference function (x_ref, y_ref)
    x_test, y_test : ndarray
        Test function (x_test, y_test)
    return_position : bool
        If True, also return the x position where the maximum relative difference occurs.
    
    Returns:
    --------
    max_rel_diff : float
        Maximum relative difference |y_test - y_ref| / |y_ref|
    x_at_max : float (only if return_position=True)
        x position where maximum relative difference occurs
    """
    # Find common range
    x_min = max(x_ref.min(), x_test.min())
    x_max = min(x_ref.max(), x_test.max())
    
    if x_min >= x_max:
        if return_position:
            return np.inf, None
        return np.inf
    
    # Create common grid
    x_common = np.linspace(x_min, x_max, min(len(x_ref), len(x_test), 1000))
    
    # Interpolate both functions
    y_ref_interp = np.interp(x_common, x_ref, y_ref)
    y_test_interp = np.interp(x_common, x_test, y_test)
    
    # Compute relative difference
    mask = np.abs(y_ref_interp) > 0
    if not np.any(mask):
        if return_position:
            return np.inf, None
        return np.inf
    
    rel_diff = np.abs((y_test_interp[mask] - y_ref_interp[mask]) / y_ref_interp[mask])
    max_idx = np.argmax(rel_diff)
    max_rel_diff = rel_diff[max_idx]
    
    if return_position:
        x_at_max = x_common[mask][max_idx]
        return max_rel_diff, x_at_max
    
    return max_rel_diff


def main():
    """Main comparison function."""
    
    # Get script directory
    script_dir = Path(__file__).parent.absolute()
    os.chdir(script_dir)
    
    # Redirect output to log file
    log_file = open("R4_compare.log", "w")
    
    def log_print(*args, **kwargs):
        """Print to both console and log file."""
        print(*args, **kwargs)
        print(*args, file=log_file, **kwargs)
    
    log_print("="*70)
    log_print("R4 Perturbation Spectra Validation: SRC vs LCDM")
    log_print("="*70)
    log_print(f"Working directory: {script_dir}")
    log_print()
    
    # Define z_pk values (from ini files: 0.,0.5,1.,2.)
    z_pk_values = [0.0, 0.5, 1.0, 2.0]
    z_pk_labels = ['z1', 'z2', 'z3', 'z4']  # CLASS names files as z1_pk.dat, z2_pk.dat, etc.
    
    # Define file paths
    base_lcdm = Path("CLASS_clean/output/R4_lcdm")
    base_src = Path("CLASS/output/R4_src")
    
    # Check if directories exist
    if not base_lcdm.exists():
        log_print(f"Error: {base_lcdm} not found")
        log_print("Please run ./run_R4_perturbation_comparison.sh first")
        log_file.close()
        sys.exit(1)
    
    if not base_src.exists():
        log_print(f"Error: {base_src} not found")
        log_print("Please run ./run_R4_perturbation_comparison.sh first")
        log_file.close()
        sys.exit(1)
    
    results = {}
    
    # 1. Compare P(k) at each redshift
    log_print("Comparing P(k,z) at multiple redshifts...")
    log_print()
    
    pk_phys_all_pass = True
    
    for z_val, z_label in zip(z_pk_values, z_pk_labels):
        log_print(f"P(k) at z = {z_val:.1f}...")
        
        pk_lcdm_file = base_lcdm / f"R4_lcdm_{z_label}_pk.dat"
        pk_src_file = base_src / f"R4_src_{z_label}_pk.dat"
        
        k_lcdm, Pk_lcdm = read_class_pk(str(pk_lcdm_file))
        k_src, Pk_src = read_class_pk(str(pk_src_file))
        
        if k_lcdm is None or k_src is None:
            log_print(f"  ✗ Error: Could not read P(k) files at z = {z_val:.1f}")
            results[f'Pk_z{z_val}'] = (False, np.inf, None, None, None, None)
            pk_phys_all_pass = False
            continue
        
        # Full range comparison
        max_delta_Pk_full, k_max_full = compute_relative_diff(
            k_lcdm, Pk_lcdm, k_src, Pk_src, return_position=True
        )
        pass_Pk_full = max_delta_Pk_full < TOL_Pk_full
        
        # Physical window: k <= 1 h/Mpc
        k_max_phys = 1.0
        mask_lcdm = k_lcdm <= k_max_phys
        mask_src = k_src <= k_max_phys
        
        if not np.any(mask_lcdm) or not np.any(mask_src):
            log_print(f"  ✗ Error: No points with k <= 1 h/Mpc at z = {z_val:.1f}")
            results[f'Pk_z{z_val}'] = (False, np.inf, None, None, None, None)
            pk_phys_all_pass = False
            continue
        
        max_delta_Pk_phys, k_max_phys_pos = compute_relative_diff(
            k_lcdm[mask_lcdm], Pk_lcdm[mask_lcdm],
            k_src[mask_src], Pk_src[mask_src],
            return_position=True
        )
        pass_Pk_phys = max_delta_Pk_phys < TOL_Pk_phys
        
        results[f'Pk_z{z_val}'] = (
            pass_Pk_phys, max_delta_Pk_phys, k_max_phys_pos,
            pass_Pk_full, max_delta_Pk_full, k_max_full
        )
        
        status_phys = "✓" if pass_Pk_phys else "✗"
        status_full = "✓" if pass_Pk_full else "✗"
        
        log_print(f"  Physical (k <= 1 h/Mpc): max ΔP/P = {max_delta_Pk_phys:.2e} at k = {k_max_phys_pos:.6e} h/Mpc (tolerance: {TOL_Pk_phys:.0e}) {status_phys}")
        log_print(f"  Full range: max ΔP/P = {max_delta_Pk_full:.2e} at k = {k_max_full:.6e} h/Mpc (tolerance: {TOL_Pk_full:.0e}) {status_full}")
        log_print()
        
        if not pass_Pk_phys:
            pk_phys_all_pass = False
    
    # 2. Compare C_ℓ^TT (lensed)
    log_print("Comparing C_ℓ^TT (lensed)...")
    
    cl_lcdm_file = base_lcdm / "R4_lcdm_cl_lensed.dat"
    cl_src_file = base_src / "R4_src_cl_lensed.dat"
    
    ell_lcdm, Cl_lcdm = read_class_cl_lensed(str(cl_lcdm_file))
    ell_src, Cl_src = read_class_cl_lensed(str(cl_src_file))
    
    if ell_lcdm is None or ell_src is None:
        log_print("  ✗ Error: Could not read CMB files")
        results['Cl'] = (False, np.inf, None, None, None, None)
        cl_phys_pass = False
    else:
        # Full range comparison
        max_delta_Cl_full, ell_max_full = compute_relative_diff(
            ell_lcdm, Cl_lcdm, ell_src, Cl_src, return_position=True
        )
        pass_Cl_full = max_delta_Cl_full < TOL_Cl_full
        
        # Physical window: ℓ <= 2000
        ell_max_phys = 2000.0
        mask_lcdm_ell = ell_lcdm <= ell_max_phys
        mask_src_ell = ell_src <= ell_max_phys
        
        if not np.any(mask_lcdm_ell) or not np.any(mask_src_ell):
            log_print("  ✗ Error: No multipoles with ℓ <= 2000")
            results['Cl'] = (False, np.inf, None, None, None, None)
            cl_phys_pass = False
        else:
            max_delta_Cl_phys, ell_max_phys_pos = compute_relative_diff(
                ell_lcdm[mask_lcdm_ell], Cl_lcdm[mask_lcdm_ell],
                ell_src[mask_src_ell], Cl_src[mask_src_ell],
                return_position=True
            )
            pass_Cl_phys = max_delta_Cl_phys < TOL_Cl_phys
            
            results['Cl'] = (
                pass_Cl_phys, max_delta_Cl_phys, ell_max_phys_pos,
                pass_Cl_full, max_delta_Cl_full, ell_max_full
            )
            
            status_phys = "✓" if pass_Cl_phys else "✗"
            status_full = "✓" if pass_Cl_full else "✗"
            
            log_print(f"  Physical (ℓ <= 2000): max ΔC_ℓ/C_ℓ = {max_delta_Cl_phys:.2e} at ℓ = {ell_max_phys_pos:.1f} (tolerance: {TOL_Cl_phys:.0e}) {status_phys}")
            log_print(f"  Full range: max ΔC_ℓ/C_ℓ = {max_delta_Cl_full:.2e} at ℓ = {ell_max_full:.1f} (tolerance: {TOL_Cl_full:.0e}) {status_full}")
            cl_phys_pass = pass_Cl_phys
    
    log_print()
    
    # Summary
    log_print("="*70)
    log_print("SUMMARY")
    log_print("="*70)
    
    log_print("P(k) physical window (k <= 1 h/Mpc):")
    for z_val, z_label in zip(z_pk_values, z_pk_labels):
        if f'Pk_z{z_val}' in results:
            r = results[f'Pk_z{z_val}']
            if r[0] is not False and r[1] is not np.inf:
                status = "✓" if r[0] else "✗"
                log_print(f"  z = {z_val:.1f}: max ΔP/P = {r[1]:.2e} at k = {r[2]:.6e} h/Mpc (tolerance: {TOL_Pk_phys:.0e}) {status}")
    
    log_print()
    log_print("C_ℓ^TT physical window (ℓ <= 2000):")
    if 'Cl' in results:
        r = results['Cl']
        if r[0] is not False and r[1] is not np.inf:
            status = "✓" if r[0] else "✗"
            log_print(f"  max ΔC_ℓ/C_ℓ = {r[1]:.2e} at ℓ = {r[2]:.1f} (tolerance: {TOL_Cl_phys:.0e}) {status}")
    
    log_print("="*70)
    
    # PASS/FAIL logic: require physical windows to pass
    all_pass = pk_phys_all_pass and cl_phys_pass
    
    if all_pass:
        log_print("RESULT: ✓ PASS - SRC perturbation spectra consistent with ΛCDM at O(1–2%) in the physical window")
        log_print("        - P(k) differences are smooth and within 2% tolerance for k ≤ 1 h/Mpc at all redshifts")
        log_print("        - C_ℓ differences are smooth and within 2% tolerance for ℓ ≤ 2000")
        log_file.close()
        return 0
    else:
        log_print("RESULT: ✗ FAIL - One or more validation criteria not met")
        if not pk_phys_all_pass:
            log_print("        - P(k) physical window exceeds tolerance at one or more redshifts")
        if not cl_phys_pass:
            log_print("        - C_ℓ physical window exceeds tolerance")
        log_file.close()
        return 1


if __name__ == "__main__":
    sys.exit(main())

