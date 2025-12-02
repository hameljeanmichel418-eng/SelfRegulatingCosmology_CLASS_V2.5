#!/usr/bin/env python3
"""
Late-time H(z) & distance diagnostics for R5 validation.

R5 Late-time H(z) & Distance Diagnostics: Quantify the late-time ΔH/H modulation
of the SRC model and its impact on cosmological distances.

Run from: ~/dev/cosmo-class/
"""

import numpy as np
import sys
import os
from pathlib import Path


# Tolerances for R5 validation
Z_MIN = 0.0
Z_MAX = 2.0
TOL_DH_MIN = 1e-4      # Minimum expected modulation (non-trivial)
TOL_DH_MAX = 5e-2      # Maximum allowed modulation
TOL_DCHI = 5e-2        # Maximum allowed comoving distance deviation


def read_class_background_with_header(filename):
    """
    Read CLASS background file and extract z, H(z), and comoving distance.
    
    Parameters:
    -----------
    filename : str
        Path to CLASS background file
    
    Returns:
    --------
    z : ndarray
        Redshift array (sorted by increasing z)
    H : ndarray
        Hubble parameter H(z) [1/Mpc] (sorted)
    chi : ndarray
        Comoving distance χ(z) [Mpc] (sorted)
    """
    if not os.path.exists(filename):
        return None, None, None
    
    # Read header to find column indices
    z_col = 0  # Default: z is always column 0
    H_col = None
    chi_col = None
    
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('#'):
                # Parse header lines like "# 1:z  2:a  3:H[1/Mpc]  4:comov. dist. [Mpc] ..."
                if ':' in line:
                    parts = line.strip().split()
                    for part in parts:
                        if ':' in part:
                            try:
                                col_num, col_name = part.split(':', 1)
                                col_idx = int(col_num) - 1  # Convert to 0-based
                                
                                # Find H column
                                if 'H' in col_name.upper() and ('1/Mpc' in col_name or 'H[' in col_name):
                                    H_col = col_idx
                                
                                # Find comoving distance column
                                if 'comov' in col_name.lower() or 'chi' in col_name.lower():
                                    chi_col = col_idx
                            except (ValueError, IndexError):
                                continue
            else:
                break
    
    # Defaults if not found in header
    if H_col is None:
        H_col = 3  # Typical default
    if chi_col is None:
        # Try to find it by checking common column positions
        # CLASS often puts comoving distance in column 4 or 5
        chi_col = 4  # Try column 4 first
    
    # Read data
    data = np.loadtxt(filename, comments="#")
    
    if len(data) == 0:
        return None, None, None
    
    # Extract columns
    z = data[:, z_col]
    H = data[:, H_col]
    
    # Try to get chi; if column doesn't exist, try alternative
    try:
        chi = data[:, chi_col]
    except IndexError:
        # Try other common positions
        if data.shape[1] > 5:
            chi_col = 5
            try:
                chi = data[:, chi_col]
            except IndexError:
                return None, None, None
    
    # Sort by increasing z
    sort_idx = np.argsort(z)
    z_sorted = z[sort_idx]
    H_sorted = H[sort_idx]
    chi_sorted = chi[sort_idx]
    
    return z_sorted, H_sorted, chi_sorted


def compute_relative_diff(x_ref, y_ref, x_test, y_test, z_min, z_max):
    """
    Compute relative difference between two functions restricted to z_min <= z <= z_max.
    
    Parameters:
    -----------
    x_ref, y_ref : ndarray
        Reference function (x_ref, y_ref) - x should be redshift z
    x_test, y_test : ndarray
        Test function (x_test, y_test) - x should be redshift z
    z_min, z_max : float
        Redshift range limits
    
    Returns:
    --------
    max_abs_rel : float
        Maximum of |Δ/Ref| = |(y_test - y_ref) / y_ref|
    z_at_max : float
        Redshift where maximum occurs
    """
    # Restrict to physical window
    mask_ref = (x_ref >= z_min) & (x_ref <= z_max)
    mask_test = (x_test >= z_min) & (x_test <= z_max)
    
    if not np.any(mask_ref) or not np.any(mask_test):
        return np.inf, None
    
    z_ref_phys = x_ref[mask_ref]
    y_ref_phys = y_ref[mask_ref]
    z_test_phys = x_test[mask_test]
    y_test_phys = y_test[mask_test]
    
    # Create common z-grid
    z_min_common = max(z_ref_phys.min(), z_test_phys.min())
    z_max_common = min(z_ref_phys.max(), z_test_phys.max())
    z_common = np.linspace(z_min_common, z_max_common, 2000)
    
    # Interpolate both functions
    y_ref_interp = np.interp(z_common, z_ref_phys, y_ref_phys)
    y_test_interp = np.interp(z_common, z_test_phys, y_test_phys)
    
    # Compute relative difference
    mask = np.abs(y_ref_interp) > 0
    if not np.any(mask):
        return np.inf, None
    
    rel_diff = np.abs((y_test_interp[mask] - y_ref_interp[mask]) / y_ref_interp[mask])
    max_idx = np.argmax(rel_diff)
    max_abs_rel = rel_diff[max_idx]
    z_at_max = z_common[mask][max_idx]
    
    return max_abs_rel, z_at_max


def distance_table(z_ref, chi_lcdm, chi_src, z_targets):
    """
    Compute luminosity distance ratios at target redshifts.
    
    Parameters:
    -----------
    z_ref : ndarray
        Redshift array for interpolation
    chi_lcdm : ndarray
        Comoving distance for LCDM [Mpc]
    chi_src : ndarray
        Comoving distance for SRC [Mpc]
    z_targets : list
        List of target redshifts
    
    Returns:
    --------
    table : list
        List of tuples (z, delta_D_L_over_D_L)
    """
    table = []
    
    for z_t in z_targets:
        # Interpolate comoving distances
        chi_lcdm_t = np.interp(z_t, z_ref, chi_lcdm)
        chi_src_t = np.interp(z_t, z_ref, chi_src)
        
        # Luminosity distance D_L = (1+z) * chi
        D_L_lcdm = (1.0 + z_t) * chi_lcdm_t
        D_L_src = (1.0 + z_t) * chi_src_t
        
        # Relative difference
        delta_D_L = (D_L_src - D_L_lcdm) / D_L_lcdm
        
        table.append((z_t, delta_D_L))
    
    return table


def main():
    """Main comparison function."""
    
    # Get script directory
    script_dir = Path(__file__).parent.absolute()
    os.chdir(script_dir)
    
    # Create output directory for data files
    output_dir = Path("R5_outputs")
    output_dir.mkdir(exist_ok=True)
    
    # Redirect output to log file
    log_file = open("R5_compare.log", "w")
    
    def log_print(*args, **kwargs):
        """Print to both console and log file."""
        print(*args, **kwargs)
        print(*args, file=log_file, **kwargs)
    
    log_print("="*70)
    log_print("R5 Late-time H(z) & Distance Diagnostics: SRC vs LCDM")
    log_print("="*70)
    log_print(f"Working directory: {script_dir}")
    log_print()
    
    # Define file paths
    bg_lcdm = Path("CLASS_clean/output/R5_lcdm/R5_lcdm_background.dat")
    bg_src = Path("CLASS/output/R5_src/R5_src_background.dat")
    
    # Check if files exist
    if not bg_lcdm.exists():
        log_print(f"Error: {bg_lcdm} not found")
        log_print("Please run ./run_R5_late_time_observables.sh first")
        log_file.close()
        sys.exit(1)
    
    if not bg_src.exists():
        log_print(f"Error: {bg_src} not found")
        log_print("Please run ./run_R5_late_time_observables.sh first")
        log_file.close()
        sys.exit(1)
    
    # Read background files
    log_print("Reading background files...")
    z_lcdm, H_lcdm, chi_lcdm = read_class_background_with_header(str(bg_lcdm))
    z_src, H_src, chi_src = read_class_background_with_header(str(bg_src))
    
    if z_lcdm is None or z_src is None:
        log_print("  ✗ Error: Could not read background files")
        log_file.close()
        sys.exit(1)
    
    if chi_lcdm is None or chi_src is None:
        log_print("  ⚠ Warning: Could not extract comoving distance column")
        log_print("    Attempting to continue with H(z) only...")
        chi_lcdm = np.zeros_like(z_lcdm)
        chi_src = np.zeros_like(z_src)
    
    log_print(f"  LCDM: {len(z_lcdm)} points, z ∈ [{z_lcdm.min():.3f}, {z_lcdm.max():.3f}]")
    log_print(f"  SRC:  {len(z_src)} points, z ∈ [{z_src.min():.3f}, {z_src.max():.3f}]")
    log_print()
    
    # Late-time physical window: 0 <= z <= 2
    log_print(f"Late-time physical window: {Z_MIN} <= z <= {Z_MAX}")
    log_print()
    
    # 1. Compute ΔH/H(z)
    log_print("Computing ΔH/H(z) = (H_SRC - H_LCDM) / H_LCDM...")
    max_abs_dH, z_at_max_dH = compute_relative_diff(z_lcdm, H_lcdm, z_src, H_src, Z_MIN, Z_MAX)
    
    if max_abs_dH == np.inf or z_at_max_dH is None:
        log_print("  ✗ Error: Could not compute ΔH/H")
        log_file.close()
        sys.exit(1)
    
    log_print(f"  max |ΔH/H| = {max_abs_dH:.2e} at z = {z_at_max_dH:.3f}")
    log_print()
    
    # 2. Compute Δχ/χ(z)
    log_print("Computing Δχ/χ(z) = (χ_SRC - χ_LCDM) / χ_LCDM...")
    max_abs_dchi, z_at_max_dchi = compute_relative_diff(z_lcdm, chi_lcdm, z_src, chi_src, Z_MIN, Z_MAX)
    
    if max_abs_dchi == np.inf or z_at_max_dchi is None:
        log_print("  ⚠ Warning: Could not compute Δχ/χ (comoving distance may not be available)")
        max_abs_dchi = 0.0
        z_at_max_dchi = 0.0
    else:
        log_print(f"  max |Δχ/χ| = {max_abs_dchi:.2e} at z = {z_at_max_dchi:.3f}")
    log_print()
    
    # 3. Compute distance ratios at sample redshifts
    log_print("Computing luminosity distance ratios at sample redshifts...")
    z_targets = [0.1, 0.25, 0.5, 1.0, 2.0]
    dist_table = distance_table(z_lcdm, chi_lcdm, chi_src, z_targets)
    
    log_print("  z      ΔD_L/D_L")
    for z_t, delta_D in dist_table:
        log_print(f"  {z_t:4.2f}   {delta_D:+.6e}")
    log_print()
    
    # 4. Export data file
    log_print("Exporting data file...")
    
    # Create common z-grid for export (restricted to physical window)
    mask_lcdm = (z_lcdm >= Z_MIN) & (z_lcdm <= Z_MAX)
    mask_src = (z_src >= Z_MIN) & (z_src <= Z_MAX)
    
    z_lcdm_phys = z_lcdm[mask_lcdm]
    H_lcdm_phys = H_lcdm[mask_lcdm]
    chi_lcdm_phys = chi_lcdm[mask_lcdm] if chi_lcdm is not None else np.zeros_like(z_lcdm_phys)
    
    z_src_phys = z_src[mask_src]
    H_src_phys = H_src[mask_src]
    chi_src_phys = chi_src[mask_src] if chi_src is not None else np.zeros_like(z_src_phys)
    
    # Create common grid
    z_min_common = max(z_lcdm_phys.min(), z_src_phys.min())
    z_max_common = min(z_lcdm_phys.max(), z_src_phys.max())
    z_export = np.linspace(z_min_common, z_max_common, 2000)
    
    # Interpolate
    H_lcdm_interp = np.interp(z_export, z_lcdm_phys, H_lcdm_phys)
    H_src_interp = np.interp(z_export, z_src_phys, H_src_phys)
    chi_lcdm_interp = np.interp(z_export, z_lcdm_phys, chi_lcdm_phys)
    chi_src_interp = np.interp(z_export, z_src_phys, chi_src_phys)
    
    # Compute relative differences
    delta_H_over_H = (H_src_interp - H_lcdm_interp) / H_lcdm_interp
    delta_chi_over_chi = (chi_src_interp - chi_lcdm_interp) / chi_lcdm_interp
    
    # Save data file
    data_file = output_dir / "R5_delta_background.dat"
    data_array = np.column_stack([
        z_export, H_lcdm_interp, H_src_interp, delta_H_over_H,
        chi_lcdm_interp, chi_src_interp, delta_chi_over_chi
    ])
    header = (
        "# z  H_LCDM[1/Mpc]  H_SRC[1/Mpc]  delta_H_over_H  "
        "chi_LCDM[Mpc]  chi_SRC[Mpc]  delta_chi_over_chi\n"
        "# Late-time H(z) & distance diagnostics (0 <= z <= 2)"
    )
    np.savetxt(data_file, data_array, header=header, fmt="%.10e")
    log_print(f"  Saved to: {data_file}")
    log_print()
    
    # 5. Check tolerances and determine PASS/FAIL
    log_print("="*70)
    log_print("VALIDATION CRITERIA")
    log_print("="*70)
    
    pass_dH = (max_abs_dH >= TOL_DH_MIN) and (max_abs_dH <= TOL_DH_MAX)
    pass_dchi = max_abs_dchi <= TOL_DCHI
    
    log_print(f"max |ΔH/H| in [{TOL_DH_MIN:.0e}, {TOL_DH_MAX:.0e}]: {max_abs_dH:.2e} at z = {z_at_max_dH:.3f} {'✓' if pass_dH else '✗'}")
    log_print(f"max |Δχ/χ| < {TOL_DCHI:.0e}: {max_abs_dchi:.2e} at z = {z_at_max_dchi:.3f} {'✓' if pass_dchi else '✗'}")
    
    all_pass = pass_dH and pass_dchi
    
    log_print("="*70)
    log_print("SUMMARY")
    log_print("="*70)
    log_print(f"max |ΔH/H| = {max_abs_dH:.2e} at z = {z_at_max_dH:.3f}")
    log_print(f"  (tolerance window: [{TOL_DH_MIN:.0e}, {TOL_DH_MAX:.0e}])")
    log_print(f"max |Δχ/χ| = {max_abs_dchi:.2e} at z = {z_at_max_dchi:.3f}")
    log_print(f"  (tolerance: {TOL_DCHI:.0e})")
    log_print()
    log_print("Luminosity distance ratios:")
    log_print("  z      ΔD_L/D_L")
    for z_t, delta_D in dist_table:
        log_print(f"  {z_t:4.2f}   {delta_D:+.6e}")
    log_print()
    log_print("="*70)
    
    if all_pass:
        log_print("RESULT: ✓ PASS - SRC late-time modulation is non-trivial and within expected range")
        log_print("        - H(z) modulation is detectable (≥ 1×10⁻⁴) and moderate (≤ 5×10⁻²)")
        log_print("        - Distance deviations remain within 5% tolerance")
        log_print(f"        - Data file saved to: {data_file}")
        log_file.close()
        return 0
    else:
        log_print("RESULT: ✗ FAIL - One or more validation criteria not met")
        if not pass_dH:
            if max_abs_dH < TOL_DH_MIN:
                log_print("        - H(z) modulation is too small (< 1×10⁻⁴)")
            else:
                log_print("        - H(z) modulation exceeds maximum tolerance (> 5×10⁻²)")
        if not pass_dchi:
            log_print("        - Comoving distance deviation exceeds tolerance (> 5×10⁻²)")
        log_file.close()
        return 1


if __name__ == "__main__":
    sys.exit(main())

