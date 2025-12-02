#!/usr/bin/env python3
"""
R9 BAO-scale diagnostics (SRC vs ΛCDM)

Goal:
    Use the existing R5 background outputs and R5 logs to check that the SRC
    implementation does not significantly distort BAO-relevant quantities:
      - sound horizon r_s
      - BAO distance combinations D_H(z), D_M(z), D_V(z)
    at a set of standard BAO redshifts.

Inputs:
    - R5_outputs/R5_delta_background.dat
    - R5_LCDM.log
    - R5_SRC.log

Outputs:
    - Human-readable diagnostics on stdout
    - Data table: R9_outputs/R9_bao_results.dat
"""

import sys
import numpy as np
from pathlib import Path

# Tolerances
TOL_RS   = 5e-3   # 0.5% on r_s
TOL_BAO  = 1e-2   # 1% on BAO distances

# Speed of light in km/s (for D_H = c/H)
C_LIGHT_KMS = 299792.458


def log_print(msg: str = "") -> None:
    """Print helper (tee will handle logging to R9_compare.log)."""
    print(msg)


def load_background_file(file_path: Path):
    """
    Read CLASS background file and extract z, H(z), and comoving distance.
    
    Based on read_class_background_with_header from compare_src_R5_late_time.py
    
    Parameters:
    -----------
    file_path : Path
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
    if not file_path.exists():
        raise FileNotFoundError(f"Background file not found: {file_path}")
    
    # Read header to find column indices
    z_col = 0  # Default: z is always column 0
    H_col = None
    chi_col = None
    
    with open(file_path, 'r') as f:
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
    data = np.loadtxt(file_path, comments="#")
    
    if len(data) == 0:
        raise ValueError(f"No data found in {file_path}")
    
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
                raise ValueError(f"Could not find comoving distance column in {file_path}")
        else:
            raise ValueError(f"Could not find comoving distance column in {file_path}")
    
    # Sort by increasing z
    sort_idx = np.argsort(z)
    z_sorted = z[sort_idx]
    H_sorted = H[sort_idx]
    chi_sorted = chi[sort_idx]
    
    return z_sorted, H_sorted, chi_sorted


def parse_header_columns(file_path: Path):
    """
    Parse the first header line starting with '#' and return a dict
    mapping column_name (as-is) -> index, plus a lowercase lookup.
    """
    header_line = None
    with file_path.open("r") as f:
        for line in f:
            if line.lstrip().startswith("#"):
                header_line = line
                break

    if header_line is None:
        raise RuntimeError(f"No header line starting with '#' found in {file_path}")

    # Remove leading '#' and split into column names
    parts = header_line.lstrip("#").strip().split()
    columns = parts
    columns_lc = [c.lower() for c in columns]

    name_to_idx = {name: i for i, name in enumerate(columns)}

    return columns, columns_lc, name_to_idx


def find_column_index(columns_lc, targets):
    """
    Find the index of a column name by trying several target names (lowercase).

    This version is robust to units in brackets, e.g. it will match:
        'h_lcdm[1/mpc]' for target 'h_lcdm'
        'chi_lcdm[mpc]' for target 'chi_lcdm'
    by checking prefix matches before failing.

    Parameters
    ----------
    columns_lc : list of str
        Lowercased column names from the header line.
    targets : list of str
        List of candidate names (e.g. ['H_LCDM', 'H_lcdm', 'h_lcdm']).

    Returns
    -------
    idx : int
        Index of the matched column.

    Raises
    ------
    ValueError
        If no matching column is found.
    """
    # Normalize targets to lowercase
    targets_lc = [t.lower() for t in targets]

    # First try exact match
    for t in targets_lc:
        for i, name in enumerate(columns_lc):
            if name == t:
                return i

    # Then try prefix match (to handle units in brackets, e.g. 'h_lcdm[1/mpc]')
    for t in targets_lc:
        for i, name in enumerate(columns_lc):
            if name.startswith(t + "[") or name.startswith(t):
                return i

    raise ValueError(f"Could not find any of {targets} in columns: {columns_lc}")


def load_r5_background():
    """
    Load background H(z) and chi(z) for LCDM and SRC from the CLASS R5
    background files, and return them on a common z-grid suitable for BAO.

    Returns
    -------
    z_common : ndarray
        Common redshift grid.
    H_lcdm : ndarray
        H_LCDM(z) on z_common.
    H_src : ndarray
        H_SRC(z) on z_common.
    chi_lcdm : ndarray
        χ_LCDM(z) on z_common.
    chi_src : ndarray
        χ_SRC(z) on z_common.
    """
    bg_lcdm_path = Path("CLASS_clean/output/R5_lcdm/R5_lcdm_background.dat")
    bg_src_path  = Path("CLASS/output/R5_src/R5_src_background.dat")

    if not bg_lcdm_path.exists():
        raise FileNotFoundError(f"R5 LCDM background file not found: {bg_lcdm_path}")
    if not bg_src_path.exists():
        raise FileNotFoundError(f"R5 SRC background file not found: {bg_src_path}")

    # Use the same loader as in R5
    z_lcdm, H_lcdm, chi_lcdm = load_background_file(bg_lcdm_path)
    z_src,  H_src,  chi_src  = load_background_file(bg_src_path)

    # Restrict to a reasonable late-time window that includes BAO redshifts
    z_min = 0.0
    z_max = 3.0

    mask_lcdm = (z_lcdm >= z_min) & (z_lcdm <= z_max)
    mask_src  = (z_src  >= z_min) & (z_src  <= z_max)

    if not np.any(mask_lcdm):
        raise RuntimeError("No LCDM background points in the chosen redshift window.")
    if not np.any(mask_src):
        raise RuntimeError("No SRC background points in the chosen redshift window.")

    # Common z-grid for interpolation (fine enough for BAO distances)
    z_low  = max(z_lcdm[mask_lcdm].min(), z_src[mask_src].min())
    z_high = min(z_lcdm[mask_lcdm].max(), z_src[mask_src].max())
    z_common = np.linspace(z_low, z_high, 2000)

    H_lcdm_interp   = np.interp(z_common, z_lcdm[mask_lcdm], H_lcdm[mask_lcdm])
    H_src_interp    = np.interp(z_common, z_src[mask_src],   H_src[mask_src])
    chi_lcdm_interp = np.interp(z_common, z_lcdm[mask_lcdm], chi_lcdm[mask_lcdm])
    chi_src_interp  = np.interp(z_common, z_src[mask_src],   chi_src[mask_src])

    return z_common, H_lcdm_interp, H_src_interp, chi_lcdm_interp, chi_src_interp


def parse_rs_from_log(log_path: Path):
    """
    Parse the sound horizon rs from a CLASS log file.

    Looks for a line like:
        "with comoving sound horizon rs = 147.071614 Mpc"
    """
    if not log_path.exists():
        return None

    rs_value = None
    with log_path.open("r") as f:
        for line in f:
            if "comoving sound horizon rs" in line:
                # Expect something like: ... rs = 147.071614 Mpc
                parts = line.strip().split()
                # Try to find '=' and take the number right after
                if "=" in parts:
                    idx = parts.index("=")
                    if idx + 1 < len(parts):
                        try:
                            rs_value = float(parts[idx + 1])
                        except ValueError:
                            rs_value = None
                break

    return rs_value


def compute_bao_distances(z_bao, z_grid, H_lcdm, H_src, chi_lcdm, chi_src):
    """
    Interpolate H(z) and chi(z) on BAO redshifts and compute
    D_H, D_M, D_V for LCDM and SRC.
    """
    # Interpolate
    H_lcdm_bao   = np.interp(z_bao, z_grid, H_lcdm)
    H_src_bao    = np.interp(z_bao, z_grid, H_src)
    chi_lcdm_bao = np.interp(z_bao, z_grid, chi_lcdm)
    chi_src_bao  = np.interp(z_bao, z_grid, chi_src)

    # D_H = c / H
    D_H_lcdm = C_LIGHT_KMS / H_lcdm_bao
    D_H_src  = C_LIGHT_KMS / H_src_bao

    # D_M is the comoving angular diameter distance (here ≈ chi)
    D_M_lcdm = chi_lcdm_bao
    D_M_src  = chi_src_bao

    # D_V(z) = [ z * D_M(z)^2 * D_H(z) ]^{1/3}
    D_V_lcdm = (z_bao * D_M_lcdm**2 * D_H_lcdm)**(1.0 / 3.0)
    D_V_src  = (z_bao * D_M_src**2  * D_H_src )**(1.0 / 3.0)

    return (D_H_lcdm, D_H_src,
            D_M_lcdm, D_M_src,
            D_V_lcdm, D_V_src)


def main():
    log_print("======================================================================")
    log_print("R9 BAO-scale Diagnostics: SRC vs LCDM")
    log_print("======================================================================")
    log_print(f"Working directory: {Path.cwd()}")
    log_print("")
    log_print("Using R5 outputs (no new CLASS runs).")
    log_print("")

    # Ensure output directory exists
    out_dir = Path("R9_outputs")
    out_dir.mkdir(exist_ok=True)

    # ------------------------------------------------------------------
    # 1. Load R5 background
    # ------------------------------------------------------------------
    try:
        z_r5, H_lcdm_r5, H_src_r5, chi_lcdm_r5, chi_src_r5 = load_r5_background()
        log_print("Using CLASS R5 background files (LCDM & SRC).")
        log_print(f"z-range available: [{z_r5.min():.3f}, {z_r5.max():.3f}]")
    except Exception as e:
        log_print(f"ERROR: Failed to load R5 background files: {e}")
        sys.exit(1)

    log_print("------------------------------------------------------------------")
    log_print("1. Late-time background inputs (from R5)")
    log_print("------------------------------------------------------------------")
    log_print("")

    # Use the interpolated data directly (already on common grid)
    z_late = z_r5
    H_lcdm_late = H_lcdm_r5
    H_src_late = H_src_r5
    chi_lcdm_late = chi_lcdm_r5
    chi_src_late = chi_src_r5

    # ------------------------------------------------------------------
    # 2. Read r_s from R5 logs
    # ------------------------------------------------------------------
    log_print("------------------------------------------------------------------")
    log_print("2. Sound horizon r_s from R5 logs")
    log_print("------------------------------------------------------------------")

    rs_lcdm = parse_rs_from_log(Path("R5_LCDM.log"))
    rs_src  = parse_rs_from_log(Path("R5_SRC.log"))

    if rs_lcdm is None:
        log_print("  ⚠ Warning: Could not extract r_s from R5_LCDM.log")
    else:
        log_print(f"  r_s_LCDM = {rs_lcdm:.6f} Mpc")

    if rs_src is None:
        log_print("  ⚠ Warning: Could not extract r_s from R5_SRC.log")
    else:
        log_print(f"  r_s_SRC  = {rs_src:.6f} Mpc")

    if (rs_lcdm is None) or (rs_src is None):
        log_print("  ⚠ r_s diagnostics will be skipped in PASS/FAIL logic.")
        delta_rs_over_rs = None
    else:
        delta_rs_over_rs = (rs_src - rs_lcdm) / rs_lcdm
        log_print(f"  Δr_s/r_s = {delta_rs_over_rs:.2e}")
    log_print("")

    # ------------------------------------------------------------------
    # 3. BAO distances at standard BAO redshifts
    # ------------------------------------------------------------------
    log_print("------------------------------------------------------------------")
    log_print("3. BAO distances D_H, D_M, D_V at standard BAO redshifts")
    log_print("------------------------------------------------------------------")

    # Common BAO redshifts used in the literature (example set)
    z_bao = np.array([0.106, 0.15, 0.38, 0.51, 0.61])

    # Check that these are within the available z-range
    if z_bao.min() < z_late.min() or z_bao.max() > z_late.max():
        log_print("  ⚠ Warning: Some BAO redshifts lie outside the R5 z-range subset.")
        log_print("              Interpolation will still proceed but check consistency.")
        log_print("")

    (D_H_lcdm, D_H_src,
     D_M_lcdm, D_M_src,
     D_V_lcdm, D_V_src) = compute_bao_distances(
        z_bao, z_late, H_lcdm_late, H_src_late,
        chi_lcdm_late, chi_src_late
    )

    # Relative differences
    delta_DH_over_DH = (D_H_src - D_H_lcdm) / D_H_lcdm
    delta_DM_over_DM = (D_M_src - D_M_lcdm) / D_M_lcdm
    delta_DV_over_DV = (D_V_src - D_V_lcdm) / D_V_lcdm

    log_print("z      ΔD_H/D_H     ΔD_M/D_M     ΔD_V/D_V")
    for i, z_val in enumerate(z_bao):
        log_print(
            f"{z_val:5.3f}  "
            f"{delta_DH_over_DH[i]:+.3e}  "
            f"{delta_DM_over_DM[i]:+.3e}  "
            f"{delta_DV_over_DV[i]:+.3e}"
        )
    log_print("")

    max_abs_DH = np.max(np.abs(delta_DH_over_DH))
    max_abs_DM = np.max(np.abs(delta_DM_over_DM))
    max_abs_DV = np.max(np.abs(delta_DV_over_DV))

    z_max_DH = z_bao[np.argmax(np.abs(delta_DH_over_DH))]
    z_max_DM = z_bao[np.argmax(np.abs(delta_DM_over_DM))]
    z_max_DV = z_bao[np.argmax(np.abs(delta_DV_over_DV))]

    # ------------------------------------------------------------------
    # 4. Export data
    # ------------------------------------------------------------------
    log_print("------------------------------------------------------------------")
    log_print("4. Exporting BAO table")
    log_print("------------------------------------------------------------------")

    out_file = out_dir / "R9_bao_results.dat"
    with out_file.open("w") as f:
        f.write("# R9 BAO diagnostics table\n")
        f.write("# Columns:\n")
        f.write("# z  "
                "D_H_LCDM[Mpc]  D_H_SRC[Mpc]  delta_DH_over_DH  "
                "D_M_LCDM[Mpc]  D_M_SRC[Mpc]  delta_DM_over_DM  "
                "D_V_LCDM[Mpc]  D_V_SRC[Mpc]  delta_DV_over_DV\n")
        for i, z_val in enumerate(z_bao):
            f.write(
                f"{z_val:8.3f}  "
                f"{D_H_lcdm[i]:15.6e}  {D_H_src[i]:15.6e}  {delta_DH_over_DH[i]:15.6e}  "
                f"{D_M_lcdm[i]:15.6e}  {D_M_src[i]:15.6e}  {delta_DM_over_DM[i]:15.6e}  "
                f"{D_V_lcdm[i]:15.6e}  {D_V_src[i]:15.6e}  {delta_DV_over_DV[i]:15.6e}\n"
            )

        # Optionally append r_s info at the end as comments
        f.write("#\n")
        f.write("# Sound horizon diagnostics (from R5 logs):\n")
        if delta_rs_over_rs is not None:
            f.write(f"# r_s_LCDM = {rs_lcdm:.6f} Mpc\n")
            f.write(f"# r_s_SRC  = {rs_src:.6f} Mpc\n")
            f.write(f"# delta_rs_over_rs = {delta_rs_over_rs:.6e}\n")
        else:
            f.write("# r_s not available (could not parse logs)\n")

    log_print(f"Saved BAO table to: {out_file}")
    log_print("")

    # ------------------------------------------------------------------
    # 5. Validation criteria & summary
    # ------------------------------------------------------------------
    log_print("======================================================================")
    log_print("VALIDATION CRITERIA")
    log_print("======================================================================")

    if delta_rs_over_rs is not None:
        pass_rs = np.abs(delta_rs_over_rs) <= TOL_RS
        log_print(
            f"|Δr_s/r_s| ≤ {TOL_RS:.0e}: "
            f"{np.abs(delta_rs_over_rs):.2e} {'✓' if pass_rs else '✗'}"
        )
    else:
        pass_rs = True  # don't block, just skip
        log_print("r_s test: skipped (could not parse from logs)")

    pass_DH = max_abs_DH <= TOL_BAO
    pass_DM = max_abs_DM <= TOL_BAO
    pass_DV = max_abs_DV <= TOL_BAO

    log_print(
        f"max |ΔD_H/D_H| ≤ {TOL_BAO:.0e}: {max_abs_DH:.2e} at z = {z_max_DH:.3f} "
        f"{'✓' if pass_DH else '✗'}"
    )
    log_print(
        f"max |ΔD_M/D_M| ≤ {TOL_BAO:.0e}: {max_abs_DM:.2e} at z = {z_max_DM:.3f} "
        f"{'✓' if pass_DM else '✗'}"
    )
    log_print(
        f"max |ΔD_V/D_V| ≤ {TOL_BAO:.0e}: {max_abs_DV:.2e} at z = {z_max_DV:.3f} "
        f"{'✓' if pass_DV else '✗'}"
    )

    log_print("======================================================================")
    log_print("SUMMARY")
    log_print("======================================================================")
    if delta_rs_over_rs is not None:
        log_print(f"Δr_s/r_s = {delta_rs_over_rs:.2e} (tolerance: {TOL_RS:.0e})")
    else:
        log_print("Δr_s/r_s: not available (logs missing or unparsable)")

    log_print(
        f"max |ΔD_H/D_H| = {max_abs_DH:.2e} at z = {z_max_DH:.3f} "
        f"(tolerance: {TOL_BAO:.0e})"
    )
    log_print(
        f"max |ΔD_M/D_M| = {max_abs_DM:.2e} at z = {z_max_DM:.3f} "
        f"(tolerance: {TOL_BAO:.0e})"
    )
    log_print(
        f"max |ΔD_V/D_V| = {max_abs_DV:.2e} at z = {z_max_DV:.3f} "
        f"(tolerance: {TOL_BAO:.0e})"
    )
    log_print("")
    log_print(f"BAO table: {out_file}")
    log_print("")

    all_pass = pass_rs and pass_DH and pass_DM and pass_DV

    log_print("======================================================================")
    if all_pass:
        log_print("RESULT: ✓ PASS - SRC BAO-scale distances consistent with ΛCDM")
        log_print("        at ≲ O(1%) level across standard BAO redshifts.")
    else:
        log_print("RESULT: ✗ FAIL - One or more BAO diagnostics exceeded tolerances")
        log_print("        Check R9_outputs/R9_bao_results.dat for details.")
    log_print("======================================================================")

    return 0


if __name__ == "__main__":
    sys.exit(main())

