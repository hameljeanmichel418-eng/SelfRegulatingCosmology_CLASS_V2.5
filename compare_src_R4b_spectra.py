#!/usr/bin/env python3
"""
Extended perturbation/spectra diagnostics for R4b validation.

R4b Extended Spectra Diagnostics: Provide deeper diagnostics and exportable
data for Volume 2.5, using the same R4 outputs.

Run from: ~/dev/cosmo-class/
"""

import numpy as np
import sys
import os
from pathlib import Path


# Numerical tolerances for R4b validation
TOL_PK_PHYS = 2e-2      # 2% tolerance on P(k) for k <= 1 h/Mpc (physical window)
TOL_PK_FULL = 0.10     # 10% tolerance on P(k) for full range (diagnostic only)
TOL_CL_PHYS = 2e-2     # 2% tolerance on C_ℓ for ℓ <= 2000 (physical window)
TOL_CL_FULL = 0.05     # 5% tolerance on C_ℓ for full range (diagnostic only)


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
    
    data = np.loadtxt(filename, comments="#")
    
    k = data[:, 0]
    Pk = data[:, 1]
    
    return k, Pk


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


def compute_rms(x, y, x_min, x_max):
    """
    Compute RMS of y over the range x_min <= x <= x_max.
    
    Parameters:
    -----------
    x : ndarray
        x values
    y : ndarray
        y values
    x_min, x_max : float
        Range limits
    
    Returns:
    --------
    rms : float
        RMS value of y in the specified range
    """
    mask = (x >= x_min) & (x <= x_max)
    if not np.any(mask):
        return None
    
    y_range = y[mask]
    rms = np.sqrt(np.mean(y_range**2))
    return rms


def main():
    """Main comparison function."""
    
    # Get script directory
    script_dir = Path(__file__).parent.absolute()
    os.chdir(script_dir)
    
    # Create output directory for data files
    output_dir = Path("R4b_outputs")
    output_dir.mkdir(exist_ok=True)
    
    # Redirect output to log file
    log_file = open("R4b_compare.log", "w")
    
    def log_print(*args, **kwargs):
        """Print to both console and log file."""
        print(*args, **kwargs)
        print(*args, file=log_file, **kwargs)
    
    log_print("="*70)
    log_print("R4b Extended Spectra Diagnostics: SRC vs LCDM")
    log_print("="*70)
    log_print(f"Working directory: {script_dir}")
    log_print()
    
    # Define z_pk values (from R4 ini files: 0.,0.5,1.,2.)
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
    log_print("="*70)
    log_print("P(k,z) Diagnostics")
    log_print("="*70)
    log_print()
    
    pk_phys_all_pass = True
    
    for z_val, z_label in zip(z_pk_values, z_pk_labels):
        log_print(f"z = {z_val:.1f}:")
        
        pk_lcdm_file = base_lcdm / f"R4_lcdm_{z_label}_pk.dat"
        pk_src_file = base_src / f"R4_src_{z_label}_pk.dat"
        
        k_lcdm, Pk_lcdm = read_class_pk(str(pk_lcdm_file))
        k_src, Pk_src = read_class_pk(str(pk_src_file))
        
        if k_lcdm is None or k_src is None:
            log_print(f"  ✗ Error: Could not read P(k) files at z = {z_val:.1f}")
            results[f'Pk_z{z_val}'] = (False, None, None, None, None, None)
            pk_phys_all_pass = False
            continue
        
        # Create common k-grid for interpolation
        k_min = max(k_lcdm.min(), k_src.min())
        k_max = min(k_lcdm.max(), k_src.max())
        k_common = np.linspace(k_min, k_max, min(len(k_lcdm), len(k_src), 2000))
        
        Pk_lcdm_interp = np.interp(k_common, k_lcdm, Pk_lcdm)
        Pk_src_interp = np.interp(k_common, k_src, Pk_src)
        
        # Compute ΔP/P
        deltaP_over_P = (Pk_src_interp - Pk_lcdm_interp) / Pk_lcdm_interp
        abs_deltaP = np.abs(deltaP_over_P)
        
        # Save data file
        data_file = output_dir / f"R4b_deltaPk_z{z_val:.1f}.dat"
        data_array = np.column_stack([k_common, Pk_lcdm_interp, Pk_src_interp, deltaP_over_P])
        header = f"# k [h/Mpc]  P_LCDM [(Mpc/h)^3]  P_SRC [(Mpc/h)^3]  deltaP_over_P\n# z = {z_val:.1f}"
        np.savetxt(data_file, data_array, header=header, fmt="%.10e")
        log_print(f"  Saved data to: {data_file}")
        
        # Physical window: k <= 1 h/Mpc
        k_max_phys = 1.0
        mask_phys = k_common <= k_max_phys
        
        if not np.any(mask_phys):
            log_print(f"  ✗ Error: No points with k <= 1 h/Mpc at z = {z_val:.1f}")
            results[f'Pk_z{z_val}'] = (False, None, None, None, None, None)
            pk_phys_all_pass = False
            continue
        
        max_deltaP_phys = np.max(abs_deltaP[mask_phys])
        idx_max_phys = np.argmax(abs_deltaP[mask_phys])
        k_max_phys_pos = k_common[mask_phys][idx_max_phys]
        
        # RMS in linear/quasi-linear band: 0.01 <= k <= 0.5 h/Mpc
        rms_deltaP_lin = compute_rms(k_common, abs_deltaP, 0.01, 0.5)
        
        # Full range
        max_deltaP_full = np.max(abs_deltaP)
        idx_max_full = np.argmax(abs_deltaP)
        k_max_full_pos = k_common[idx_max_full]
        
        pass_Pk_phys = max_deltaP_phys < TOL_PK_PHYS
        pass_Pk_full = max_deltaP_full < TOL_PK_FULL
        
        results[f'Pk_z{z_val}'] = (
            pass_Pk_phys, max_deltaP_phys, k_max_phys_pos, rms_deltaP_lin,
            max_deltaP_full, k_max_full_pos
        )
        
        status_phys = "✓" if pass_Pk_phys else "✗"
        status_full = "✓" if pass_Pk_full else "✗"
        
        log_print(f"  Physical window (k <= 1 h/Mpc):")
        log_print(f"    max |ΔP/P| = {max_deltaP_phys:.2e} at k = {k_max_phys_pos:.6e} h/Mpc (tolerance: {TOL_PK_PHYS:.0e}) {status_phys}")
        if rms_deltaP_lin is not None:
            log_print(f"    RMS |ΔP/P| (0.01 <= k <= 0.5 h/Mpc) = {rms_deltaP_lin:.2e}")
        log_print(f"  Full range (diagnostic):")
        log_print(f"    max |ΔP/P| = {max_deltaP_full:.2e} at k = {k_max_full_pos:.6e} h/Mpc (tolerance: {TOL_PK_FULL:.0e}) {status_full}")
        log_print()
        
        if not pass_Pk_phys:
            pk_phys_all_pass = False
    
    # 2. Compare CMB spectra (TT, EE, TE, and PP if available)
    log_print("="*70)
    log_print("C_ℓ Diagnostics")
    log_print("="*70)
    log_print()
    
    cl_spectra = ['TT', 'EE', 'TE']
    cl_spectra_validation = ['TT', 'EE']  # Only TT and EE are used for PASS/FAIL
    cl_lcdm_file = base_lcdm / "R4_lcdm_cl_lensed.dat"
    cl_src_file = base_src / "R4_src_cl_lensed.dat"
    
    # Check if lensing potential (PP) is available
    ell_test, _ = read_class_cl_lensed(str(cl_lcdm_file), spectrum='PP')
    if ell_test is not None:
        cl_spectra.append('PP')
    
    cl_phys_all_pass = True
    
    for spectrum in cl_spectra:
        log_print(f"C_ℓ^{spectrum}:")
        
        ell_lcdm, Cl_lcdm = read_class_cl_lensed(str(cl_lcdm_file), spectrum=spectrum)
        ell_src, Cl_src = read_class_cl_lensed(str(cl_src_file), spectrum=spectrum)
        
        if ell_lcdm is None or ell_src is None:
            log_print(f"  ✗ Error: Could not read C_ℓ^{spectrum} data")
            results[f'Cl_{spectrum}'] = (False, None, None, None, None, None)
            cl_phys_all_pass = False
            continue
        
        # Create common ell-grid for interpolation
        ell_min = max(ell_lcdm.min(), ell_src.min())
        ell_max = min(ell_lcdm.max(), ell_src.max())
        ell_common = np.linspace(ell_min, ell_max, min(len(ell_lcdm), len(ell_src), 2000))
        
        Cl_lcdm_interp = np.interp(ell_common, ell_lcdm, Cl_lcdm)
        Cl_src_interp = np.interp(ell_common, ell_src, Cl_src)
        
        # Compute ΔC_ℓ/C_ℓ
        deltaC_over_C = (Cl_src_interp - Cl_lcdm_interp) / Cl_lcdm_interp
        abs_deltaC = np.abs(deltaC_over_C)
        
        # For TE, also compute absolute differences normalized by peak |TE|
        is_te = (spectrum == 'TE')
        if is_te:
            abs_deltaC_abs = np.abs(Cl_src_interp - Cl_lcdm_interp)
            peak_te = np.max(np.abs(Cl_lcdm_interp))
            if peak_te > 0:
                abs_deltaC_normalized = abs_deltaC_abs / peak_te
            else:
                abs_deltaC_normalized = None
        
        # Save data file (skip PP for now as it might not be standard)
        if spectrum != 'PP':
            data_file = output_dir / f"R4b_deltaCl_{spectrum}.dat"
            data_array = np.column_stack([ell_common, Cl_lcdm_interp, Cl_src_interp, deltaC_over_C])
            header = f"# ell  Cl_LCDM  Cl_SRC  deltaC_over_C\n# Spectrum: {spectrum}"
            np.savetxt(data_file, data_array, header=header, fmt="%.10e")
            log_print(f"  Saved data to: {data_file}")
        
        # Physical window: ℓ <= 2000
        ell_max_phys = 2000.0
        mask_phys_ell = ell_common <= ell_max_phys
        
        if not np.any(mask_phys_ell):
            log_print(f"  ✗ Error: No multipoles with ℓ <= 2000")
            results[f'Cl_{spectrum}'] = (False, None, None, None, None, None, is_te)
            if spectrum in cl_spectra_validation:
                cl_phys_all_pass = False
            continue
        
        max_deltaC_phys = np.max(abs_deltaC[mask_phys_ell])
        idx_max_phys = np.argmax(abs_deltaC[mask_phys_ell])
        ell_max_phys_pos = ell_common[mask_phys_ell][idx_max_phys]
        
        # RMS over 30 <= ℓ <= 2000
        rms_deltaC = compute_rms(ell_common, abs_deltaC, 30.0, 2000.0)
        
        # Full range
        max_deltaC_full = np.max(abs_deltaC)
        idx_max_full = np.argmax(abs_deltaC)
        ell_max_full_pos = ell_common[idx_max_full]
        
        # For validation spectra, check tolerance; TE is diagnostic only
        if spectrum in cl_spectra_validation:
            pass_Cl_phys = max_deltaC_phys < TOL_CL_PHYS
            pass_Cl_full = max_deltaC_full < TOL_CL_FULL
        else:
            # TE and PP are diagnostic only
            pass_Cl_phys = True  # Don't fail on TE
            pass_Cl_full = True
        
        results[f'Cl_{spectrum}'] = (
            pass_Cl_phys, max_deltaC_phys, ell_max_phys_pos, rms_deltaC,
            max_deltaC_full, ell_max_full_pos, is_te
        )
        
        status_phys = "✓" if pass_Cl_phys else "✗"
        status_full = "✓" if pass_Cl_full else "✗"
        
        if is_te:
            log_print(f"  Physical window (ℓ <= 2000) [DIAGNOSTIC ONLY]:")
            log_print(f"    max |ΔC_ℓ/C_ℓ| = {max_deltaC_phys:.2e} at ℓ = {ell_max_phys_pos:.1f} (tolerance: {TOL_CL_PHYS:.0e}) {status_phys}")
            if rms_deltaC is not None:
                log_print(f"    RMS |ΔC_ℓ/C_ℓ| (30 <= ℓ <= 2000) = {rms_deltaC:.2e}")
            if abs_deltaC_normalized is not None:
                max_abs_norm = np.max(abs_deltaC_normalized[mask_phys_ell])
                log_print(f"    max |ΔC_ℓ|/peak|TE| = {max_abs_norm:.2e} (absolute metric)")
            if max_deltaC_phys >= TOL_CL_PHYS:
                log_print(f"    ⚠ Warning: TE relative differences are large due to zero crossings; TE is used for diagnostics only, not for validation.")
        else:
            log_print(f"  Physical window (ℓ <= 2000):")
            log_print(f"    max |ΔC_ℓ/C_ℓ| = {max_deltaC_phys:.2e} at ℓ = {ell_max_phys_pos:.1f} (tolerance: {TOL_CL_PHYS:.0e}) {status_phys}")
            if rms_deltaC is not None:
                log_print(f"    RMS |ΔC_ℓ/C_ℓ| (30 <= ℓ <= 2000) = {rms_deltaC:.2e}")
        
        log_print(f"  Full range (diagnostic):")
        log_print(f"    max |ΔC_ℓ/C_ℓ| = {max_deltaC_full:.2e} at ℓ = {ell_max_full_pos:.1f} (tolerance: {TOL_CL_FULL:.0e}) {status_full}")
        log_print()
        
        # Only check validation spectra for PASS/FAIL
        if spectrum in cl_spectra_validation and not pass_Cl_phys:
            cl_phys_all_pass = False
    
    # Summary
    log_print("="*70)
    log_print("SUMMARY")
    log_print("="*70)
    
    log_print("P(k) physical window (k <= 1 h/Mpc):")
    for z_val in z_pk_values:
        if f'Pk_z{z_val}' in results:
            r = results[f'Pk_z{z_val}']
            if r[0] is not False and r[1] is not None:
                status = "✓" if r[0] else "✗"
                log_print(f"  z = {z_val:.1f}: max |ΔP/P| = {r[1]:.2e} at k = {r[2]:.6e} h/Mpc (tolerance: {TOL_PK_PHYS:.0e}) {status}")
                if r[3] is not None:
                    log_print(f"            RMS |ΔP/P| (0.01 <= k <= 0.5) = {r[3]:.2e}")
    
    log_print()
    log_print("C_ℓ physical window (ℓ <= 2000):")
    for spectrum in cl_spectra:
        if f'Cl_{spectrum}' in results:
            r = results[f'Cl_{spectrum}']
            if r[0] is not False and r[1] is not None:
                is_te = r[6] if len(r) > 6 else False
                if is_te:
                    log_print(f"  {spectrum}: max |ΔC_ℓ/C_ℓ| = {r[1]:.2e} at ℓ = {r[2]:.1f} [DIAGNOSTIC ONLY]")
                    if r[3] is not None:
                        log_print(f"         RMS |ΔC_ℓ/C_ℓ| (30 <= ℓ <= 2000) = {r[3]:.2e}")
                else:
                    status = "✓" if r[0] else "✗"
                    log_print(f"  {spectrum}: max |ΔC_ℓ/C_ℓ| = {r[1]:.2e} at ℓ = {r[2]:.1f} (tolerance: {TOL_CL_PHYS:.0e}) {status}")
                    if r[3] is not None:
                        log_print(f"         RMS |ΔC_ℓ/C_ℓ| (30 <= ℓ <= 2000) = {r[3]:.2e}")
    
    log_print("="*70)
    
    # PASS/FAIL logic: require physical windows to pass
    all_pass = pk_phys_all_pass and cl_phys_all_pass
    
    if all_pass:
        log_print("RESULT: ✓ PASS - SRC perturbation spectra consistent with ΛCDM at O(1–2%) in the physical window")
        log_print("        - P(k) differences are smooth and within 2% tolerance for k ≤ 1 h/Mpc at all redshifts")
        log_print("        - C_ℓ^TT and C_ℓ^EE differences are smooth and within 2% tolerance for ℓ ≤ 2000")
        # Check if TE has large differences
        if 'Cl_TE' in results:
            r_te = results['Cl_TE']
            if r_te[1] is not None and r_te[1] >= TOL_CL_PHYS:
                log_print("        - Note: TE relative differences are large due to zero crossings; TE is diagnostic only, not used in validation")
        log_print(f"        - Data files saved to: {output_dir}/")
        log_file.close()
        return 0
    else:
        log_print("RESULT: ✗ FAIL - One or more validation criteria not met")
        if not pk_phys_all_pass:
            log_print("        - P(k) physical window exceeds tolerance at one or more redshifts")
        if not cl_phys_all_pass:
            log_print("        - C_ℓ physical window exceeds tolerance for TT or EE")
        log_file.close()
        return 1


if __name__ == "__main__":
    sys.exit(main())

