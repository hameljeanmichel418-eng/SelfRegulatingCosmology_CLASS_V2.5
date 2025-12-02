#!/usr/bin/env python3
"""
Linear growth & fσ8 diagnostics for R6 validation.

R6 Linear Growth & fσ8 Diagnostics: Compare growth factor D(z), growth rate f(z),
σ8(z), and fσ8(z) between SRC and ΛCDM.

Run from: ~/dev/cosmo-class/
"""

import numpy as np
import sys
import os
from pathlib import Path


# Tolerances for R6 validation
TOL_D = 0.02      # 2% tolerance on ΔD/D
TOL_F = 0.03      # 3% tolerance on Δf/f
TOL_SIGMA8 = 0.03 # 3% tolerance on Δσ8/σ8
TOL_FSIGMA8 = 0.03 # 3% tolerance on Δ(fσ8)/(fσ8) for z ≤ 1.5

# Reference scale for growth factor
K0 = 0.1  # h/Mpc

# Top-hat filter radius
R8 = 8.0  # h⁻¹ Mpc


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


def compute_sigma8(k, Pk, R=R8):
    """
    Compute σ_R from P(k) using top-hat window function.
    
    Parameters:
    -----------
    k : ndarray
        Wavenumber array [h/Mpc]
    Pk : ndarray
        Power spectrum P(k) [(Mpc/h)^3]
    R : float
        Top-hat radius [h⁻¹ Mpc], default 8.0
    
    Returns:
    --------
    sigma_R : float
        σ_R value
    """
    # Window function: W(x) = 3 (sin x - x cos x) / x^3
    # where x = k*R
    x = k * R
    
    # Avoid division by zero at x=0
    mask = x > 1e-10
    W_sq = np.zeros_like(x)
    W_sq[mask] = (3.0 * (np.sin(x[mask]) - x[mask] * np.cos(x[mask])) / x[mask]**3)**2
    W_sq[~mask] = 1.0  # W(0) = 1
    
    # Integrate: σ_R^2 = (1/2π^2) ∫ dk k^2 P(k) W^2(kR)
    # Use trapezoidal rule
    integrand = k**2 * Pk * W_sq
    dk = np.diff(k)
    
    # Handle non-uniform k-grid
    if len(dk) > 0:
        # Use trapezoidal rule
        sigma_sq = (1.0 / (2.0 * np.pi**2)) * np.trapz(integrand, k)
    else:
        sigma_sq = 0.0
    
    sigma_R = np.sqrt(max(0.0, sigma_sq))
    return sigma_R


def compute_growth_rate(z_array, D_array):
    """
    Compute the linear growth rate f(z) = d ln D / d ln a.

    Parameters
    ----------
    z_array : array-like
        Redshift values (monotonically increasing).
    D_array : array-like
        Growth factor D(z), same length as z_array.

    Returns
    -------
    f_array : ndarray
        Growth rate f(z) evaluated at the same z points.
    """
    # Convert to numpy arrays
    z_array = np.asarray(z_array, dtype=float)
    D_array = np.asarray(D_array, dtype=float)

    # Scale factor a = 1 / (1+z) → ln a = -ln(1+z)
    ln_a = -np.log1p(z_array)

    # Ensure monotonic ordering
    sort_idx = np.argsort(ln_a)
    ln_a_sorted = ln_a[sort_idx]
    ln_D_sorted = np.log(D_array[sort_idx])

    # f = d ln D / d ln a
    f_sorted = np.gradient(ln_D_sorted, ln_a_sorted)

    # Return in original z order
    inv_idx = np.argsort(sort_idx)
    return f_sorted[inv_idx]


def main():
    """Main comparison function."""
    
    # Get script directory
    script_dir = Path(__file__).parent.absolute()
    os.chdir(script_dir)
    
    # Create output directory for data files
    output_dir = Path("R6_outputs")
    output_dir.mkdir(exist_ok=True)
    
    # Redirect output to log file
    log_file = open("R6_compare.log", "w")
    
    def log_print(*args, **kwargs):
        """Print to both console and log file."""
        print(*args, **kwargs)
        print(*args, file=log_file, **kwargs)
    
    log_print("="*70)
    log_print("R6 Linear Growth & fσ8 Diagnostics: SRC vs LCDM")
    log_print("="*70)
    log_print(f"Working directory: {script_dir}")
    log_print()
    
    # Define z_pk values (from R4: 0.,0.5,1.,2.)
    z_pk_values = [0.0, 0.5, 1.0, 2.0]
    z_pk_labels = ['z1', 'z2', 'z3', 'z4']  # CLASS names files as z1_pk.dat, etc.
    
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
    
    # Read P(k,z) files for all redshifts
    log_print("Reading P(k,z) files from R4 outputs...")
    log_print()
    
    Pk_lcdm_dict = {}
    Pk_src_dict = {}
    k_dict = {}
    
    for z_val, z_label in zip(z_pk_values, z_pk_labels):
        pk_lcdm_file = base_lcdm / f"R4_lcdm_{z_label}_pk.dat"
        pk_src_file = base_src / f"R4_src_{z_label}_pk.dat"
        
        k_lcdm, Pk_lcdm = read_class_pk(str(pk_lcdm_file))
        k_src, Pk_src = read_class_pk(str(pk_src_file))
        
        if k_lcdm is None or k_src is None:
            log_print(f"✗ Error: Could not read P(k) files at z = {z_val:.1f}")
            log_file.close()
            sys.exit(1)
        
        # Use LCDM k-grid as reference (should be the same for both)
        k_dict[z_val] = k_lcdm
        Pk_lcdm_dict[z_val] = Pk_lcdm
        Pk_src_dict[z_val] = np.interp(k_lcdm, k_src, Pk_src)  # Interpolate SRC onto LCDM grid
    
    log_print(f"Loaded P(k,z) for z = {z_pk_values}")
    log_print()
    
    # 1. Compute linear growth factor D(z)
    log_print("="*70)
    log_print("1. Linear Growth Factor D(z)")
    log_print("="*70)
    log_print(f"Using reference scale k0 = {K0} h/Mpc")
    log_print()
    
    D_lcdm = []
    D_src = []
    
    for z_val in z_pk_values:
        k = k_dict[z_val]
        Pk_lcdm = Pk_lcdm_dict[z_val]
        Pk_src = Pk_src_dict[z_val]
        
        # Interpolate P(k) at k0
        Pk_lcdm_k0 = np.interp(K0, k, Pk_lcdm)
        Pk_src_k0 = np.interp(K0, k, Pk_src)
        
        # D(z) ∝ sqrt(P(k0,z))
        D_lcdm.append(np.sqrt(Pk_lcdm_k0))
        D_src.append(np.sqrt(Pk_src_k0))
    
    D_lcdm = np.array(D_lcdm)
    D_src = np.array(D_src)
    
    # Normalize so that D(z=0) = 1
    D_lcdm = D_lcdm / D_lcdm[0]
    D_src = D_src / D_src[0]
    
    D_ratio = D_src / D_lcdm
    delta_D_over_D = D_ratio - 1.0
    
    log_print("z      D_LCDM    D_SRC     D_ratio   ΔD/D")
    for i, z_val in enumerate(z_pk_values):
        log_print(f"{z_val:4.1f}   {D_lcdm[i]:.6f}  {D_src[i]:.6f}  {D_ratio[i]:.6f}  {delta_D_over_D[i]:+.6e}")
    log_print()
    
    max_abs_delta_D = np.max(np.abs(delta_D_over_D))
    idx_max_D = np.argmax(np.abs(delta_D_over_D))
    z_at_max_D = z_pk_values[idx_max_D]
    
    log_print(f"max |ΔD/D| = {max_abs_delta_D:.2e} at z = {z_at_max_D:.1f} (tolerance: {TOL_D:.0e})")
    log_print()
    
    # Check monotonicity
    if np.any(np.diff(D_lcdm) > 0) or np.any(np.diff(D_src) > 0):
        log_print("⚠ Warning: D(z) is not monotonically decreasing (should decrease as z increases)")
        log_print()
    
    # 2. Compute growth rate f(z)
    log_print("="*70)
    log_print("2. Growth Rate f(z) = d ln D / d ln a")
    log_print("="*70)
    log_print()
    
    f_lcdm = compute_growth_rate(z_pk_values, D_lcdm)
    f_src = compute_growth_rate(z_pk_values, D_src)
    
    f_ratio = f_src / f_lcdm
    delta_f_over_f = f_ratio - 1.0
    
    log_print("z      f_LCDM    f_SRC     f_ratio   Δf/f")
    for i, z_val in enumerate(z_pk_values):
        log_print(f"{z_val:4.1f}   {f_lcdm[i]:.6f}  {f_src[i]:.6f}  {f_ratio[i]:.6f}  {delta_f_over_f[i]:+.6e}")
    log_print()
    
    max_abs_delta_f = np.max(np.abs(delta_f_over_f))
    idx_max_f = np.argmax(np.abs(delta_f_over_f))
    z_at_max_f = z_pk_values[idx_max_f]
    
    log_print(f"max |Δf/f| = {max_abs_delta_f:.2e} at z = {z_at_max_f:.1f} (tolerance: {TOL_F:.0e})")
    log_print()
    
    # 3. Compute σ8(z)
    log_print("="*70)
    log_print("3. σ8(z) from Top-Hat Filter")
    log_print("="*70)
    log_print(f"Using top-hat radius R = {R8} h⁻¹ Mpc")
    log_print()
    
    sigma8_lcdm = []
    sigma8_src = []
    
    for z_val in z_pk_values:
        k = k_dict[z_val]
        Pk_lcdm = Pk_lcdm_dict[z_val]
        Pk_src = Pk_src_dict[z_val]
        
        sigma8_lcdm.append(compute_sigma8(k, Pk_lcdm, R=R8))
        sigma8_src.append(compute_sigma8(k, Pk_src, R=R8))
    
    sigma8_lcdm = np.array(sigma8_lcdm)
    sigma8_src = np.array(sigma8_src)
    
    sigma8_ratio = sigma8_src / sigma8_lcdm
    delta_sigma8_over_sigma8 = sigma8_ratio - 1.0
    
    log_print("z      σ8_LCDM   σ8_SRC    σ8_ratio  Δσ8/σ8")
    for i, z_val in enumerate(z_pk_values):
        log_print(f"{z_val:4.1f}   {sigma8_lcdm[i]:.6f}  {sigma8_src[i]:.6f}  {sigma8_ratio[i]:.6f}  {delta_sigma8_over_sigma8[i]:+.6e}")
    log_print()
    
    max_abs_delta_sigma8 = np.max(np.abs(delta_sigma8_over_sigma8))
    idx_max_sigma8 = np.argmax(np.abs(delta_sigma8_over_sigma8))
    z_at_max_sigma8 = z_pk_values[idx_max_sigma8]
    
    log_print(f"max |Δσ8/σ8| = {max_abs_delta_sigma8:.2e} at z = {z_at_max_sigma8:.1f} (tolerance: {TOL_SIGMA8:.0e})")
    log_print()
    
    # 4. Compute fσ8(z)
    log_print("="*70)
    log_print("4. fσ8(z) = f(z) × σ8(z)")
    log_print("="*70)
    log_print()
    
    fsigma8_lcdm = f_lcdm * sigma8_lcdm
    fsigma8_src = f_src * sigma8_src
    
    fsigma8_ratio = fsigma8_src / fsigma8_lcdm
    delta_fsigma8_over_fsigma8 = fsigma8_ratio - 1.0
    
    log_print("z      fσ8_LCDM  fσ8_SRC   fσ8_ratio Δ(fσ8)/(fσ8)")
    for i, z_val in enumerate(z_pk_values):
        log_print(f"{z_val:4.1f}   {fsigma8_lcdm[i]:.6f}  {fsigma8_src[i]:.6f}  {fsigma8_ratio[i]:.6f}  {delta_fsigma8_over_fsigma8[i]:+.6e}")
    log_print()
    
    max_abs_delta_fsigma8 = np.max(np.abs(delta_fsigma8_over_fsigma8))
    idx_max_fsigma8 = np.argmax(np.abs(delta_fsigma8_over_fsigma8))
    z_at_max_fsigma8 = z_pk_values[idx_max_fsigma8]
    
    log_print(f"max |Δ(fσ8)/(fσ8)| = {max_abs_delta_fsigma8:.2e} at z = {z_at_max_fsigma8:.1f} (tolerance: {TOL_FSIGMA8:.0e})")
    
    # For fσ8, check z ≤ 1.5 separately
    z_pk_array = np.array(z_pk_values)
    mask_z15 = z_pk_array <= 1.5
    
    if np.any(mask_z15):
        # Restrict to z ≤ 1.5
        delta_fsigma8_z15 = delta_fsigma8_over_fsigma8[mask_z15]
        
        max_abs_delta_fsigma8_z15 = np.max(np.abs(delta_fsigma8_z15))
        idx_max_fsigma8_z15 = np.argmax(np.abs(delta_fsigma8_z15))
        z_at_max_fsigma8_z15 = z_pk_array[mask_z15][idx_max_fsigma8_z15]
        log_print(f"  For z ≤ 1.5: max |Δ(fσ8)/(fσ8)| = {max_abs_delta_fsigma8_z15:.2e} at z = {z_at_max_fsigma8_z15:.1f}")
    else:
        max_abs_delta_fsigma8_z15 = None
        z_at_max_fsigma8_z15 = None
    
    log_print()
    
    # 5. Save data file
    log_print("="*70)
    log_print("Exporting data file...")
    log_print("="*70)
    
    data_file = output_dir / "R6_growth_table.dat"
    header = (
        "# z  D_lcdm  D_src  D_ratio  f_lcdm  f_src  f_ratio  "
        "sigma8_lcdm  sigma8_src  sigma8_ratio  fsigma8_lcdm  fsigma8_src  fsigma8_ratio\n"
        "# Linear growth & fσ8 diagnostics (SRC vs LCDM)"
    )
    
    data_array = np.column_stack([
        z_pk_values, D_lcdm, D_src, D_ratio,
        f_lcdm, f_src, f_ratio,
        sigma8_lcdm, sigma8_src, sigma8_ratio,
        fsigma8_lcdm, fsigma8_src, fsigma8_ratio
    ])
    
    np.savetxt(data_file, data_array, header=header, fmt="%.10e")
    log_print(f"Saved to: {data_file}")
    log_print()
    
    # 6. Check tolerances and determine PASS/FAIL
    log_print("="*70)
    log_print("VALIDATION CRITERIA")
    log_print("="*70)
    
    pass_D = max_abs_delta_D <= TOL_D
    pass_f = max_abs_delta_f <= TOL_F
    pass_sigma8 = max_abs_delta_sigma8 <= TOL_SIGMA8
    pass_fsigma8 = (max_abs_delta_fsigma8_z15 <= TOL_FSIGMA8) if (max_abs_delta_fsigma8_z15 is not None) else True
    
    log_print(f"max |ΔD/D| ≤ {TOL_D:.0e}: {max_abs_delta_D:.2e} at z = {z_at_max_D:.1f} {'✓' if pass_D else '✗'}")
    log_print(f"max |Δf/f| ≤ {TOL_F:.0e}: {max_abs_delta_f:.2e} at z = {z_at_max_f:.1f} {'✓' if pass_f else '✗'}")
    log_print(f"max |Δσ8/σ8| ≤ {TOL_SIGMA8:.0e}: {max_abs_delta_sigma8:.2e} at z = {z_at_max_sigma8:.1f} {'✓' if pass_sigma8 else '✗'}")
    if max_abs_delta_fsigma8_z15 is not None:
        log_print(f"max |Δ(fσ8)/(fσ8)| (z ≤ 1.5) ≤ {TOL_FSIGMA8:.0e}: {max_abs_delta_fsigma8_z15:.2e} at z = {z_at_max_fsigma8_z15:.1f} {'✓' if pass_fsigma8 else '✗'}")
    
    all_pass = pass_D and pass_f and pass_sigma8 and pass_fsigma8
    
    log_print("="*70)
    log_print("SUMMARY")
    log_print("="*70)
    log_print(f"max |ΔD/D| = {max_abs_delta_D:.2e} at z = {z_at_max_D:.1f} (tolerance: {TOL_D:.0e})")
    log_print(f"max |Δf/f| = {max_abs_delta_f:.2e} at z = {z_at_max_f:.1f} (tolerance: {TOL_F:.0e})")
    log_print(f"max |Δσ8/σ8| = {max_abs_delta_sigma8:.2e} at z = {z_at_max_sigma8:.1f} (tolerance: {TOL_SIGMA8:.0e})")
    if max_abs_delta_fsigma8_z15 is not None:
        log_print(f"max |Δ(fσ8)/(fσ8)| (z ≤ 1.5) = {max_abs_delta_fsigma8_z15:.2e} at z = {z_at_max_fsigma8_z15:.1f} (tolerance: {TOL_FSIGMA8:.0e})")
    log_print()
    log_print("Growth table:")
    log_print("  z      ΔD/D      Δf/f      Δσ8/σ8    Δ(fσ8)/(fσ8)")
    for i, z_val in enumerate(z_pk_values):
        log_print(f"  {z_val:4.1f}   {delta_D_over_D[i]:+.6e}  {delta_f_over_f[i]:+.6e}  {delta_sigma8_over_sigma8[i]:+.6e}  {delta_fsigma8_over_fsigma8[i]:+.6e}")
    log_print("="*70)
    
    if all_pass:
        log_print("RESULT: ✓ PASS - SRC linear growth and fσ8 differ from ΛCDM at the O(1–3%) level")
        log_print("        in the late-time window, consistent with a small dissipative modulation.")
        log_print(f"        - Data file saved to: {data_file}")
        log_file.close()
        return 0
    else:
        log_print("RESULT: ✗ FAIL - One or more validation criteria not met")
        if not pass_D:
            log_print("        - Growth factor D(z) deviation exceeds tolerance")
        if not pass_f:
            log_print("        - Growth rate f(z) deviation exceeds tolerance")
        if not pass_sigma8:
            log_print("        - σ8(z) deviation exceeds tolerance")
        if not pass_fsigma8:
            log_print("        - fσ8(z) deviation exceeds tolerance for z ≤ 1.5")
        log_file.close()
        return 1


if __name__ == "__main__":
    sys.exit(main())

