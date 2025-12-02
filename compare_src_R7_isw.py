#!/usr/bin/env python3
"""
R7 ISW & low-ℓ CMB Diagnostics: Validate that SRC does not spoil low-ℓ TT power.

R7 uses existing R4/R4b CMB spectra to check ISW-sensitive regime (ℓ ≤ 50).

Run from: ~/dev/cosmo-class/
"""

import numpy as np
import sys
import os
from pathlib import Path


# Numerical tolerances for R7 validation
TOL_LOW_MAX = 2e-2   # 2% max deviation in 2≤ℓ≤50
TOL_LOW_RMS = 1e-2   # 1% RMS in 2≤ℓ≤50
TOL_VLOW_MAX = 3e-2  # 3% max deviation in 2≤ℓ≤10 (more lenient due to cosmic variance)


def read_class_cl_lensed(filename, spectrum='TT'):
    """
    Read CLASS lensed CMB power spectrum file.
    
    Parameters:
    -----------
    filename : str
        Path to CLASS cl_lensed file
    spectrum : str
        Spectrum to extract: 'TT', 'EE', 'TE', or 'PP' (lensing potential)
    
    Returns:
    --------
    ell : ndarray
        Multipole array
    Cl : ndarray
        C_ℓ values for the requested spectrum
    """
    if not os.path.exists(filename):
        return None, None
    
    with open(filename, 'r') as f:
        # Find the header line with column names
        spectrum_idx = None
        for line in f:
            if line.startswith('#') and ':' in line:
                # Parse column names
                parts = line.strip().split()
                for part in parts:
                    if ':' in part:
                        try:
                            col_num, col_name = part.split(':', 1)
                            # Check if this is the spectrum we want
                            if spectrum.upper() in col_name.upper():
                                col_idx = int(col_num) - 1  # Convert to 0-based index
                                spectrum_idx = col_idx
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


def interpolate_cl(ell_ref, cl_ref, ell_test, cl_test):
    """
    Interpolate test Cl spectrum onto reference ell grid.
    
    Parameters:
    -----------
    ell_ref : ndarray
        Reference multipole grid
    cl_ref : ndarray
        Reference Cl values (not used, for consistency with interface)
    ell_test : ndarray
        Test multipole grid
    cl_test : ndarray
        Test Cl values
    
    Returns:
    --------
    cl_test_interp : ndarray
        Test Cl interpolated onto ell_ref grid
    """
    cl_test_interp = np.interp(ell_ref, ell_test, cl_test)
    return cl_test_interp


def rel_diff(ref, test):
    """
    Compute relative difference (test - ref) / ref.
    
    Parameters:
    -----------
    ref : ndarray
        Reference values
    test : ndarray
        Test values
    
    Returns:
    --------
    delta : ndarray
        Relative difference array
    """
    mask = np.abs(ref) > 0
    delta = np.zeros_like(ref)
    delta[mask] = (test[mask] - ref[mask]) / ref[mask]
    return delta


def log_print(msg=""):
    """Print to both stdout and log file."""
    print(msg)
    if hasattr(log_print, 'log_file'):
        log_print.log_file.write(msg + '\n')


def main():
    """Main driver for R7 ISW diagnostics."""
    
    # Set up working directory
    work_dir = Path.cwd()
    log_print("="*70)
    log_print("R7 ISW & low-ℓ TT Diagnostics: SRC vs LCDM")
    log_print("="*70)
    log_print(f"Working directory: {work_dir}")
    log_print()
    
    # Base directories
    base_lcdm = Path("CLASS_clean/output/R4_lcdm")
    base_src = Path("CLASS/output/R4_src")
    
    # Input files
    cl_file_lcdm = base_lcdm / "R4_lcdm_cl_lensed.dat"
    cl_file_src = base_src / "R4_src_cl_lensed.dat"
    
    # Check files exist
    if not cl_file_lcdm.exists():
        log_print(f"ERROR: Could not find LCDM Cl file: {cl_file_lcdm}")
        sys.exit(1)
    if not cl_file_src.exists():
        log_print(f"ERROR: Could not find SRC Cl file: {cl_file_src}")
        sys.exit(1)
    
    log_print("Loading lensed Cl_TT spectra...")
    log_print(f"  LCDM: {cl_file_lcdm}")
    log_print(f"  SRC:  {cl_file_src}")
    log_print()
    
    # Load Cl spectra
    ell_lcdm, Cl_lcdm = read_class_cl_lensed(str(cl_file_lcdm), spectrum='TT')
    ell_src, Cl_src = read_class_cl_lensed(str(cl_file_src), spectrum='TT')
    
    if ell_lcdm is None or Cl_lcdm is None:
        log_print("ERROR: Failed to load LCDM Cl_TT spectrum")
        sys.exit(1)
    if ell_src is None or Cl_src is None:
        log_print("ERROR: Failed to load SRC Cl_TT spectrum")
        sys.exit(1)
    
    log_print(f"Loaded LCDM: {len(ell_lcdm)} multipoles, ℓ ∈ [{ell_lcdm.min():.0f}, {ell_lcdm.max():.0f}]")
    log_print(f"Loaded SRC:  {len(ell_src)} multipoles, ℓ ∈ [{ell_src.min():.0f}, {ell_src.max():.0f}]")
    log_print()
    
    # Interpolate SRC onto LCDM ell grid
    Cl_src_interp = interpolate_cl(ell_lcdm, Cl_lcdm, ell_src, Cl_src)
    
    # Use LCDM ell grid as reference
    ell = ell_lcdm
    
    # Compute relative differences
    delta = rel_diff(Cl_lcdm, Cl_src_interp)
    
    log_print("="*70)
    log_print("Low-ℓ TT Diagnostics (lensed Cl_TT)")
    log_print("="*70)
    log_print()
    
    # Define windows
    mask_2_10 = (ell >= 2) & (ell <= 10)
    mask_2_50 = (ell >= 2) & (ell <= 50)
    mask_2_30 = (ell >= 2) & (ell <= 30)
    
    # Window 2-10
    if np.any(mask_2_10):
        delta_2_10 = delta[mask_2_10]
        ell_2_10 = ell[mask_2_10]
        max_abs_2_10 = np.max(np.abs(delta_2_10))
        idx_max_2_10 = np.argmax(np.abs(delta_2_10))
        ell_at_max_2_10 = ell_2_10[idx_max_2_10]
        rms_2_10 = np.sqrt(np.mean(delta_2_10**2))
    else:
        max_abs_2_10 = None
        ell_at_max_2_10 = None
        rms_2_10 = None
    
    # Window 2-50
    if np.any(mask_2_50):
        delta_2_50 = delta[mask_2_50]
        ell_2_50 = ell[mask_2_50]
        max_abs_2_50 = np.max(np.abs(delta_2_50))
        idx_max_2_50 = np.argmax(np.abs(delta_2_50))
        ell_at_max_2_50 = ell_2_50[idx_max_2_50]
        rms_2_50 = np.sqrt(np.mean(delta_2_50**2))
    else:
        max_abs_2_50 = None
        ell_at_max_2_50 = None
        rms_2_50 = None
    
    # Print diagnostics
    log_print("ℓ in [2,10]:")
    if max_abs_2_10 is not None:
        log_print(f"  max |ΔC_ℓ/C_ℓ| = {max_abs_2_10:.6e} at ℓ = {ell_at_max_2_10:.0f}")
        log_print(f"  RMS |ΔC_ℓ/C_ℓ| = {rms_2_10:.6e}")
    else:
        log_print("  No data in this window")
    log_print()
    
    log_print("ℓ in [2,50]:")
    if max_abs_2_50 is not None:
        log_print(f"  max |ΔC_ℓ/C_ℓ| = {max_abs_2_50:.6e} at ℓ = {ell_at_max_2_50:.0f}")
        log_print(f"  RMS |ΔC_ℓ/C_ℓ| = {rms_2_50:.6e}")
    else:
        log_print("  No data in this window")
    log_print()
    
    # ISW amplitude proxy (2 ≤ ℓ ≤ 30)
    if np.any(mask_2_30):
        ell_2_30 = ell[mask_2_30]
        Cl_lcdm_2_30 = Cl_lcdm[mask_2_30]
        Cl_src_2_30 = Cl_src_interp[mask_2_30]
        
        # Compute D_ℓ = ℓ(ℓ+1)C_ℓ / (2π)
        Dl_lcdm = ell_2_30 * (ell_2_30 + 1) * Cl_lcdm_2_30 / (2 * np.pi)
        Dl_src = ell_2_30 * (ell_2_30 + 1) * Cl_src_2_30 / (2 * np.pi)
        
        # ISW amplitude ratio
        A_isw = np.sum(Dl_src) / np.sum(Dl_lcdm)
        
        log_print("ISW amplitude proxy (2 ≤ ℓ ≤ 30):")
        log_print(f"  A_ISW = Σ D_ℓ^SRC / Σ D_ℓ^LCDM = {A_isw:.6e}")
        log_print()
    else:
        A_isw = None
        log_print("ISW amplitude proxy: No data in 2≤ℓ≤30 window")
        log_print()
    
    # Check tolerances
    pass_max_2_50 = (max_abs_2_50 is not None) and (max_abs_2_50 <= TOL_LOW_MAX)
    pass_rms_2_50 = (rms_2_50 is not None) and (rms_2_50 <= TOL_LOW_RMS)
    pass_max_2_10 = (max_abs_2_10 is not None) and (max_abs_2_10 <= TOL_VLOW_MAX)
    
    all_pass = pass_max_2_50 and pass_rms_2_50 and pass_max_2_10
    
    # Summary
    log_print("="*70)
    log_print("SUMMARY")
    log_print("="*70)
    if max_abs_2_10 is not None:
        status_2_10 = '✓' if pass_max_2_10 else '✗'
        log_print(f"max |ΔC_ℓ/C_ℓ| (2≤ℓ≤10)  = {max_abs_2_10:.6e} (tolerance: {TOL_VLOW_MAX:.0e}) {status_2_10}")
    if max_abs_2_50 is not None:
        status_2_50 = '✓' if pass_max_2_50 else '✗'
        log_print(f"max |ΔC_ℓ/C_ℓ| (2≤ℓ≤50)  = {max_abs_2_50:.6e} (tolerance: {TOL_LOW_MAX:.0e}) {status_2_50}")
    if rms_2_50 is not None:
        status_rms = '✓' if pass_rms_2_50 else '✗'
        log_print(f"RMS |ΔC_ℓ/C_ℓ| (2≤ℓ≤50) = {rms_2_50:.6e} (tolerance: {TOL_LOW_RMS:.0e}) {status_rms}")
    if A_isw is not None:
        log_print(f"A_ISW (2≤ℓ≤30)           = {A_isw:.6e}")
    log_print("="*70)
    log_print()
    
    # Result
    if all_pass:
        log_print("RESULT: ✓ PASS - SRC low-ℓ TT power remains within tolerances")
        log_print("        of ΛCDM in the ISW-sensitive regime.")
    else:
        log_print("RESULT: ✗ FAIL - One or more validation criteria not met.")
        if not pass_max_2_50:
            log_print(f"        - max |ΔC_ℓ/C_ℓ| (2≤ℓ≤50) = {max_abs_2_50:.6e} > {TOL_LOW_MAX:.0e}")
        if not pass_rms_2_50:
            log_print(f"        - RMS |ΔC_ℓ/C_ℓ| (2≤ℓ≤50) = {rms_2_50:.6e} > {TOL_LOW_RMS:.0e}")
        if not pass_max_2_10:
            log_print(f"        - max |ΔC_ℓ/C_ℓ| (2≤ℓ≤10) = {max_abs_2_10:.6e} > {TOL_VLOW_MAX:.0e}")
    log_print()
    
    # Export data file
    log_print("="*70)
    log_print("Exporting data file...")
    log_print("="*70)
    
    output_dir = Path("R7_outputs")
    output_dir.mkdir(exist_ok=True)
    
    data_file = output_dir / "R7_deltaCl_lowL.dat"
    
    # Restrict to low-ℓ range for export (ℓ ≤ 100 for Volume 2.5)
    mask_export = ell <= 100
    ell_export = ell[mask_export]
    Cl_lcdm_export = Cl_lcdm[mask_export]
    Cl_src_export = Cl_src_interp[mask_export]
    delta_export = delta[mask_export]
    
    header = (
        "# ell  Cl_lcdm  Cl_src  delta_rel\n"
        "# Relative difference: delta_rel = (Cl_src - Cl_lcdm) / Cl_lcdm\n"
        "# R7 ISW & low-ℓ TT diagnostics (SRC vs LCDM)"
    )
    
    data_array = np.column_stack([ell_export, Cl_lcdm_export, Cl_src_export, delta_export])
    np.savetxt(data_file, data_array, header=header, fmt="%.10e")
    
    log_print(f"Saved to: {data_file}")
    log_print()
    
    # Exit code
    sys.exit(0 if all_pass else 1)


if __name__ == "__main__":
    main()

