#!/usr/bin/env python3
"""
Compare SRC base precision vs high-accuracy precision outputs.

R2 Precision-Stability Validation: Verify that increasing precision
settings does not change physical results beyond expected numerical
convergence limits.

Run from: ~/dev/cosmo-class/
"""

import numpy as np
import sys
import os
from pathlib import Path


# Numerical tolerances for R2 validation
# These are more relaxed than R1 since we're comparing different precision settings
TOL_BG = 1e-3      # max ΔH/H or relative background differences

# Full-range tolerances (diagnostics only, not blockers for PASS)
TOL_PK = 1e-1      # 10% for the full k-range (sanity check only)
TOL_CL = 1e-2      # 1% for the full ℓ-range (sanity check only)

# Physical window tolerances (these determine PASS/FAIL)
TOL_PK_PHYS = 2e-2   # 2% tolerance on P(k) for k <= 1 h/Mpc
TOL_CL_PHYS = 1e-2   # 1% tolerance on C_ℓ for ℓ <= 2000


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
    log_file = open("R2_compare.log", "w")
    
    def log_print(*args, **kwargs):
        """Print to both console and log file."""
        print(*args, **kwargs)
        print(*args, file=log_file, **kwargs)
    
    log_print("="*70)
    log_print("R2 Precision-Stability Validation: Base vs High-Accuracy SRC")
    log_print("="*70)
    log_print(f"Working directory: {script_dir}")
    log_print()
    
    # Define file paths
    base_base = Path("CLASS/output/R2_base")
    base_hiacc = Path("CLASS/output/R2_hiacc")
    
    # Check if directories exist
    if not base_base.exists():
        log_print(f"Error: {base_base} not found")
        log_print("Please run ./run_R2_precision_stability.sh first")
        log_file.close()
        sys.exit(1)
    
    if not base_hiacc.exists():
        log_print(f"Error: {base_hiacc} not found")
        log_print("Please run ./run_R2_precision_stability.sh first")
        log_file.close()
        sys.exit(1)
    
    # File paths
    bg_base = base_base / "R2_base_background.dat"
    bg_hiacc = base_hiacc / "R2_hiacc_background.dat"
    
    # Note: CLASS names files as z1_pk.dat for first z_pk value, z2_pk.dat for second
    # With z_pk = 0.,1.: z1=0, z2=1
    pk_base_z0 = base_base / "R2_base_z1_pk.dat"  # z1 corresponds to first z_pk value (0.0)
    pk_hiacc_z0 = base_hiacc / "R2_hiacc_z1_pk.dat"
    
    pk_base_z1 = base_base / "R2_base_z2_pk.dat"  # z2 corresponds to second z_pk value (1.0)
    pk_hiacc_z1 = base_hiacc / "R2_hiacc_z2_pk.dat"
    
    cl_base = base_base / "R2_base_cl_lensed.dat"
    cl_hiacc = base_hiacc / "R2_hiacc_cl_lensed.dat"
    
    results = {}
    
    # 1. Compare H(z)
    log_print("Comparing H(z) from background files...")
    z_base, H_base = read_class_background(str(bg_base))
    z_hiacc, H_hiacc = read_class_background(str(bg_hiacc))
    
    if z_base is None or z_hiacc is None:
        log_print("  ✗ Error: Could not read background files")
        results['H'] = (False, np.inf, "File read error")
    else:
        max_delta_H = compute_relative_diff(z_base, H_base, z_hiacc, H_hiacc)
        pass_H = max_delta_H < TOL_BG
        results['H'] = (pass_H, max_delta_H, f"max ΔH/H = {max_delta_H:.2e}")
        status = "✓" if pass_H else "✗"
        log_print(f"  {status} {results['H'][2]} (tolerance: {TOL_BG:.0e})")
    
    # 2. Compare P(k) at z=0
    log_print("\nComparing P(k) at z=0...")
    k_base_z0, Pk_base_z0 = read_class_pk(str(pk_base_z0))
    k_hiacc_z0, Pk_hiacc_z0 = read_class_pk(str(pk_hiacc_z0))
    
    if k_base_z0 is None or k_hiacc_z0 is None:
        log_print("  ✗ Error: Could not read P(k) files at z=0")
        results['Pk_z0'] = (False, np.inf, "File read error", None)
    else:
        max_delta_Pk_z0, k_max_z0 = compute_relative_diff(k_base_z0, Pk_base_z0, k_hiacc_z0, Pk_hiacc_z0, return_position=True)
        pass_Pk_z0 = max_delta_Pk_z0 < TOL_PK
        results['Pk_z0'] = (pass_Pk_z0, max_delta_Pk_z0, f"max ΔP/P = {max_delta_Pk_z0:.2e}", k_max_z0)
        status = "✓" if pass_Pk_z0 else "✗"
        log_print(f"  {status} {results['Pk_z0'][2]} (tolerance: {TOL_PK:.0e})")
        log_print(f"    Maximum occurs at k = {k_max_z0:.6e} h/Mpc")
    
    # 2b. Physical window for P(k): k <= 1 h/Mpc
    log_print("\nP(k) physical window (k <= 1 h/Mpc)...")
    if k_base_z0 is None or k_hiacc_z0 is None:
        log_print("  ✗ Error: P(k) data missing, cannot compute physical window")
        results['Pk_phys'] = (False, np.inf, "File read error", None)
    else:
        # Select k <= 1 h/Mpc in both runs
        k_max_phys = 1.0
        mask_base = k_base_z0 <= k_max_phys
        mask_hiacc = k_hiacc_z0 <= k_max_phys
        
        if not np.any(mask_base) or not np.any(mask_hiacc):
            log_print("  ✗ Error: No points with k <= 1 h/Mpc in P(k) files")
            results['Pk_phys'] = (False, np.inf, "No k <= 1 h/Mpc", None)
        else:
            max_delta_Pk_phys, k_max_phys_pos = compute_relative_diff(
                k_base_z0[mask_base], Pk_base_z0[mask_base],
                k_hiacc_z0[mask_hiacc], Pk_hiacc_z0[mask_hiacc],
                return_position=True
            )
            pass_Pk_phys = max_delta_Pk_phys < TOL_PK_PHYS
            results['Pk_phys'] = (pass_Pk_phys, max_delta_Pk_phys,
                                  f"max ΔP/P (k <= 1 h/Mpc) = {max_delta_Pk_phys:.2e}",
                                  k_max_phys_pos)
            status = "✓" if pass_Pk_phys else "✗"
            log_print(f"  {status} {results['Pk_phys'][2]} (tolerance: {TOL_PK_PHYS:.0e})")
            log_print(f"      at k = {k_max_phys_pos:.6e} h/Mpc")
    
    # 3. Compare P(k) at z=1
    log_print("\nComparing P(k) at z=1...")
    k_base_z1, Pk_base_z1 = read_class_pk(str(pk_base_z1))
    k_hiacc_z1, Pk_hiacc_z1 = read_class_pk(str(pk_hiacc_z1))
    
    if k_base_z1 is None or k_hiacc_z1 is None:
        log_print("  ✗ Error: Could not read P(k) files at z=1")
        results['Pk_z1'] = (False, np.inf, "File read error", None)
    else:
        max_delta_Pk_z1, k_max_z1 = compute_relative_diff(k_base_z1, Pk_base_z1, k_hiacc_z1, Pk_hiacc_z1, return_position=True)
        pass_Pk_z1 = max_delta_Pk_z1 < TOL_PK
        results['Pk_z1'] = (pass_Pk_z1, max_delta_Pk_z1, f"max ΔP/P = {max_delta_Pk_z1:.2e}", k_max_z1)
        status = "✓" if pass_Pk_z1 else "✗"
        log_print(f"  {status} {results['Pk_z1'][2]} (tolerance: {TOL_PK:.0e})")
        log_print(f"    Maximum occurs at k = {k_max_z1:.6e} h/Mpc")
    
    # 4. Compare C_ℓ^TT (lensed)
    log_print("\nComparing C_ℓ^TT (lensed) from CMB files...")
    ell_base, Cl_base = read_class_cl_lensed(str(cl_base))
    ell_hiacc, Cl_hiacc = read_class_cl_lensed(str(cl_hiacc))
    
    if ell_base is None or ell_hiacc is None:
        log_print("  ✗ Error: Could not read CMB files")
        results['Cl'] = (False, np.inf, "File read error", None)
    else:
        max_delta_Cl, ell_max = compute_relative_diff(ell_base, Cl_base, ell_hiacc, Cl_hiacc, return_position=True)
        pass_Cl = max_delta_Cl < TOL_CL
        results['Cl'] = (pass_Cl, max_delta_Cl, f"max ΔC_ℓ/C_ℓ = {max_delta_Cl:.2e}", ell_max)
        status = "✓" if pass_Cl else "✗"
        log_print(f"  {status} {results['Cl'][2]} (tolerance: {TOL_CL:.0e})")
        log_print(f"    Maximum occurs at ℓ = {ell_max:.1f}")
    
    # 4b. Physical window for C_ℓ: ℓ <= 2000
    log_print("\nC_ℓ^TT physical window (ℓ <= 2000)...")
    if ell_base is None or ell_hiacc is None:
        log_print("  ✗ Error: CMB data missing, cannot compute physical window")
        results['Cl_phys'] = (False, np.inf, "File read error", None)
    else:
        ell_max_phys = 2000.0
        mask_base_ell = ell_base <= ell_max_phys
        mask_hiacc_ell = ell_hiacc <= ell_max_phys
        
        if not np.any(mask_base_ell) or not np.any(mask_hiacc_ell):
            log_print("  ✗ Error: No multipoles with ℓ <= 2000 in CMB files")
            results['Cl_phys'] = (False, np.inf, "No ℓ <= 2000", None)
        else:
            max_delta_Cl_phys, ell_max_phys_pos = compute_relative_diff(
                ell_base[mask_base_ell], Cl_base[mask_base_ell],
                ell_hiacc[mask_hiacc_ell], Cl_hiacc[mask_hiacc_ell],
                return_position=True
            )
            pass_Cl_phys = max_delta_Cl_phys < TOL_CL_PHYS
            results['Cl_phys'] = (pass_Cl_phys, max_delta_Cl_phys,
                                  f"max ΔC_ℓ/C_ℓ (ℓ <= 2000) = {max_delta_Cl_phys:.2e}",
                                  ell_max_phys_pos)
            status = "✓" if pass_Cl_phys else "✗"
            log_print(f"  {status} {results['Cl_phys'][2]} (tolerance: {TOL_CL_PHYS:.0e})")
            log_print(f"      at ℓ = {ell_max_phys_pos:.1f}")
    
    # Summary
    log_print("\n" + "="*70)
    log_print("SUMMARY")
    log_print("="*70)
    
    # PASS/FAIL logic: require H, physical P(k), and physical C_ℓ to pass
    # Full-range comparisons are diagnostics only
    pass_H = results['H'][0]
    pass_Pk_phys = results.get('Pk_phys', (False, np.inf, "Missing", None))[0]
    pass_Cl_phys = results.get('Cl_phys', (False, np.inf, "Missing", None))[0]
    
    all_pass = pass_H and pass_Pk_phys and pass_Cl_phys
    
    log_print(f"max ΔH/H        = {results['H'][1]:.2e} (tolerance: {TOL_BG:.0e}) {'✓' if results['H'][0] else '✗'}")
    log_print(f"max ΔP/P (z=0)  = {results['Pk_z0'][1]:.2e} (tolerance: {TOL_PK:.0e}) {'✓' if results['Pk_z0'][0] else '✗'}")
    if len(results['Pk_z0']) > 3 and results['Pk_z0'][3] is not None:
        log_print(f"                at k = {results['Pk_z0'][3]:.6e} h/Mpc")
    log_print(f"max ΔP/P (z=1)  = {results['Pk_z1'][1]:.2e} (tolerance: {TOL_PK:.0e}) {'✓' if results['Pk_z1'][0] else '✗'}")
    if len(results['Pk_z1']) > 3 and results['Pk_z1'][3] is not None:
        log_print(f"                at k = {results['Pk_z1'][3]:.6e} h/Mpc")
    log_print(f"max ΔC_ℓ/C_ℓ    = {results['Cl'][1]:.2e} (tolerance: {TOL_CL:.0e}) {'✓' if results['Cl'][0] else '✗'}")
    if len(results['Cl']) > 3 and results['Cl'][3] is not None:
        log_print(f"                at ℓ = {results['Cl'][3]:.1f}")
    
    # New: physical window summary
    if 'Pk_phys' in results:
        log_print(f"max ΔP/P (phys, k<=1)  = {results['Pk_phys'][1]:.2e} (tolerance: {TOL_PK_PHYS:.0e}) {'✓' if results['Pk_phys'][0] else '✗'}")
        if results['Pk_phys'][3] is not None:
            log_print(f"                       at k = {results['Pk_phys'][3]:.6e} h/Mpc")
    if 'Cl_phys' in results:
        log_print(f"max ΔC_ℓ/C_ℓ (phys, ℓ<=2000) = {results['Cl_phys'][1]:.2e} (tolerance: {TOL_CL_PHYS:.0e}) {'✓' if results['Cl_phys'][0] else '✗'}")
        if results['Cl_phys'][3] is not None:
            log_print(f"                              at ℓ = {results['Cl_phys'][3]:.1f}")
    
    log_print("="*70)
    
    if all_pass:
        log_print("RESULT: ✓ PASS - SRC model is stable under precision changes")
        log_print("        Differences are within expected numerical convergence limits")
        log_file.close()
        return 0
    else:
        log_print("RESULT: ✗ FAIL - Differences exceed tolerances")
        log_print("        This may indicate precision-related issues or implementation problems")
        log_file.close()
        return 1


if __name__ == "__main__":
    sys.exit(main())

