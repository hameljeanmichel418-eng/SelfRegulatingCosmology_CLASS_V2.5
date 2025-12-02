#!/usr/bin/env python3
"""
R8 Global Observables Summary: Consolidated diagnostics for SRC vs ΛCDM.

R8 reuses outputs from R4b, R5, R6, and R7 to summarize:
- Late-time background modulation ΔH/H and Δχ/χ (R5)
- Growth, σ8(z), and fσ8(z) (R6)
- Linear P(k,z) and C_ell (TT,EE) spectra differences (R4b)
- Low-ℓ TT and ISW amplitude proxy (R7)

No new CLASS runs are performed. R8 only reads the existing data files.

Run from: ~/dev/cosmo-class/
"""

import numpy as np
import sys
import os
from pathlib import Path


# Numerical tolerances (reused from previous validation blocks)
TOL_DH = 5e-2         # max |ΔH/H|
TOL_DCHI = 5e-2       # max |Δχ/χ|
TOL_D = 2e-2          # max |ΔD/D|
TOL_F = 3e-2          # max |Δf/f|
TOL_SIGMA8 = 3e-2     # max |Δσ8/σ8|
TOL_FSIGMA8 = 3e-2    # max |Δ(fσ8)/(fσ8)| (z <= 1.5)
TOL_PK_PHYS = 2e-2    # P(k) physical window
TOL_CL_PHYS = 2e-2    # TT/EE physical window ℓ <= 2000
TOL_LOWELL_MAX = 3e-2 # low-ℓ TT max (2<=ℓ<=10)
TOL_LOWELL_RMS = 1e-2 # low-ℓ TT RMS (2<=ℓ<=50)
TOL_AISW_REL = 5e-2   # |A_ISW - 1|


def log_print(msg=""):
    """Print to both stdout and log file."""
    print(msg)
    if hasattr(log_print, 'log_file'):
        log_print.log_file.write(msg + '\n')


def main():
    """Main driver for R8 global observables summary."""
    
    # Set up working directory and log file
    work_dir = Path.cwd()
    log_file = open("R8_compare.log", "w")
    log_print.log_file = log_file
    
    log_print("="*70)
    log_print("R8 Global Observables Summary: SRC vs ΛCDM")
    log_print("="*70)
    log_print(f"Working directory: {work_dir}")
    log_print()
    log_print("Reusing outputs from R4b, R5, R6, and R7 (no new CLASS runs).")
    log_print()
    
    # Initialize results dictionary
    results = {}
    
    # =======================================================================
    # Section 1: Background (R5)
    # =======================================================================
    log_print("="*70)
    log_print("1. Late-time Background & Distances (R5)")
    log_print("="*70)
    log_print()
    
    # Check for R5 file (try both possible names)
    r5_file = Path("R5_outputs/R5_late_time_observables.dat")
    if not r5_file.exists():
        r5_file_alt = Path("R5_outputs/R5_delta_background.dat")
        if r5_file_alt.exists():
            r5_file = r5_file_alt
            log_print(f"Note: Using R5 file: {r5_file}")
        else:
            log_print(f"ERROR: Could not find R5 file: {r5_file} or {r5_file_alt}")
            log_print("Please run ./run_R5_late_time_observables.sh first")
            log_file.close()
            sys.exit(1)
    
    try:
        data_r5 = np.loadtxt(str(r5_file), comments="#")
        # Expected columns: z, H_LCDM, H_SRC, delta_H_over_H, chi_LCDM, chi_SRC, delta_chi_over_chi
        z_r5 = data_r5[:, 0]
        delta_H_over_H_r5 = data_r5[:, 3]  # column 3 (0-indexed)
        delta_chi_over_chi_r5 = data_r5[:, 6]  # column 6 (0-indexed)
        
        max_abs_dH = np.max(np.abs(delta_H_over_H_r5))
        idx_max_dH = np.argmax(np.abs(delta_H_over_H_r5))
        z_at_max_dH = z_r5[idx_max_dH]
        
        max_abs_dchi = np.max(np.abs(delta_chi_over_chi_r5))
        idx_max_dchi = np.argmax(np.abs(delta_chi_over_chi_r5))
        z_at_max_dchi = z_r5[idx_max_dchi]
        
        # Extract D_L sample points if available (or compute from R5 data)
        # For now, we'll note that D_L diagnostics were computed in R5
        # but may need to be read from a separate section of the file
        
        log_print(f"max |ΔH/H| = {max_abs_dH:.6e} at z = {z_at_max_dH:.3f}")
        log_print(f"max |Δχ/χ| = {max_abs_dchi:.6e} at z = {z_at_max_dchi:.3f}")
        log_print()
        
        results['bg'] = {
            'max_dH': max_abs_dH,
            'z_max_dH': z_at_max_dH,
            'max_dchi': max_abs_dchi,
            'z_max_dchi': z_at_max_dchi,
            'z_array': z_r5,
            'delta_H': delta_H_over_H_r5,
            'delta_chi': delta_chi_over_chi_r5
        }
        
    except Exception as e:
        log_print(f"ERROR: Failed to read R5 file: {e}")
        log_file.close()
        sys.exit(1)
    
    # =======================================================================
    # Section 2: Growth & fσ8 (R6)
    # =======================================================================
    log_print("="*70)
    log_print("2. Linear Growth & fσ8 (R6)")
    log_print("="*70)
    log_print()
    
    r6_file = Path("R6_outputs/R6_growth_table.dat")
    if not r6_file.exists():
        log_print(f"ERROR: Could not find R6 file: {r6_file}")
        log_print("Please run ./run_R6_growth_diagnostics.sh first")
        log_file.close()
        sys.exit(1)
    
    try:
        data_r6 = np.loadtxt(str(r6_file), comments="#")
        # Expected columns: z, D_lcdm, D_src, D_ratio, f_lcdm, f_src, f_ratio,
        #                   sigma8_lcdm, sigma8_src, sigma8_ratio,
        #                   fsigma8_lcdm, fsigma8_src, fsigma8_ratio
        z_r6 = data_r6[:, 0]
        D_ratio = data_r6[:, 3]  # D_src/D_lcdm
        f_ratio = data_r6[:, 6]  # f_src/f_lcdm
        sigma8_ratio = data_r6[:, 9]  # sigma8_src/sigma8_lcdm
        fsigma8_ratio = data_r6[:, 12]  # fsigma8_src/fsigma8_lcdm
        
        # Compute relative differences
        delta_D_over_D = D_ratio - 1.0
        delta_f_over_f = f_ratio - 1.0
        delta_sigma8_over_sigma8 = sigma8_ratio - 1.0
        delta_fsigma8_over_fsigma8 = fsigma8_ratio - 1.0
        
        max_abs_dD = np.max(np.abs(delta_D_over_D))
        idx_max_dD = np.argmax(np.abs(delta_D_over_D))
        z_at_max_dD = z_r6[idx_max_dD]
        
        max_abs_df = np.max(np.abs(delta_f_over_f))
        idx_max_df = np.argmax(np.abs(delta_f_over_f))
        z_at_max_df = z_r6[idx_max_df]
        
        max_abs_dsigma8 = np.max(np.abs(delta_sigma8_over_sigma8))
        idx_max_dsigma8 = np.argmax(np.abs(delta_sigma8_over_sigma8))
        z_at_max_dsigma8 = z_r6[idx_max_dsigma8]
        
        # For fσ8, check z <= 1.5 separately
        mask_z15 = z_r6 <= 1.5
        if np.any(mask_z15):
            delta_fsigma8_z15 = delta_fsigma8_over_fsigma8[mask_z15]
            z_r6_z15 = z_r6[mask_z15]
            max_abs_dfsigma8_z15 = np.max(np.abs(delta_fsigma8_z15))
            idx_max_dfsigma8_z15 = np.argmax(np.abs(delta_fsigma8_z15))
            z_at_max_dfsigma8_z15 = z_r6_z15[idx_max_dfsigma8_z15]
        else:
            max_abs_dfsigma8_z15 = None
            z_at_max_dfsigma8_z15 = None
        
        log_print(f"max |ΔD/D| = {max_abs_dD:.6e} at z = {z_at_max_dD:.1f}")
        log_print(f"max |Δf/f| = {max_abs_df:.6e} at z = {z_at_max_df:.1f}")
        log_print(f"max |Δσ8/σ8| = {max_abs_dsigma8:.6e} at z = {z_at_max_dsigma8:.1f}")
        if max_abs_dfsigma8_z15 is not None:
            log_print(f"max |Δ(fσ8)/(fσ8)| (z ≤ 1.5) = {max_abs_dfsigma8_z15:.6e} at z = {z_at_max_dfsigma8_z15:.1f}")
        log_print()
        
        results['growth'] = {
            'max_dD': max_abs_dD,
            'z_max_dD': z_at_max_dD,
            'max_df': max_abs_df,
            'z_max_df': z_at_max_df,
            'max_dsigma8': max_abs_dsigma8,
            'z_max_dsigma8': z_at_max_dsigma8,
            'max_dfsigma8_z15': max_abs_dfsigma8_z15,
            'z_max_dfsigma8_z15': z_at_max_dfsigma8_z15,
            'z_array': z_r6,
            'delta_D': delta_D_over_D,
            'delta_f': delta_f_over_f,
            'delta_sigma8': delta_sigma8_over_sigma8,
            'delta_fsigma8': delta_fsigma8_over_fsigma8
        }
        
    except Exception as e:
        log_print(f"ERROR: Failed to read R6 file: {e}")
        log_file.close()
        sys.exit(1)
    
    # =======================================================================
    # Section 3: Spectra (P(k), Cℓ) (R4b)
    # =======================================================================
    log_print("="*70)
    log_print("3. Spectra: P(k,z) & C_ℓ (R4b)")
    log_print("="*70)
    log_print()
    
    # P(k) files
    z_pk_values = [0.0, 0.5, 1.0, 2.0]
    # R4b actually creates files with format: R4b_deltaPk_z0.0.dat, R4b_deltaPk_z0.5.dat, etc.
    z_pk_file_suffixes = ['z0.0', 'z0.5', 'z1.0', 'z2.0']
    
    max_abs_dPk_phys_all = []
    z_at_max_dPk_phys_all = []
    k_at_max_dPk_phys_all = []
    
    for z_val, z_suffix in zip(z_pk_values, z_pk_file_suffixes):
        pk_file = Path(f"R4b_outputs/R4b_deltaPk_{z_suffix}.dat")
        if not pk_file.exists():
            log_print(f"WARNING: Could not find P(k) file at z={z_val}: {pk_file}")
            continue
        
        try:
            data_pk = np.loadtxt(str(pk_file), comments="#")
            # Expected columns: k, P_LCDM, P_SRC, deltaP_over_P
            k_pk = data_pk[:, 0]
            deltaP_over_P_pk = data_pk[:, 3]
            
            # Physical window: k <= 1 h/Mpc
            mask_phys = k_pk <= 1.0
            if np.any(mask_phys):
                abs_deltaP_phys = np.abs(deltaP_over_P_pk[mask_phys])
                k_phys = k_pk[mask_phys]
                max_abs_dPk_phys = np.max(abs_deltaP_phys)
                idx_max = np.argmax(abs_deltaP_phys)
                k_at_max = k_phys[idx_max]
                
                max_abs_dPk_phys_all.append(max_abs_dPk_phys)
                z_at_max_dPk_phys_all.append(z_val)
                k_at_max_dPk_phys_all.append(k_at_max)
        except Exception as e:
            log_print(f"WARNING: Failed to read P(k) file at z={z_val}: {e}")
            continue
    
    if max_abs_dPk_phys_all:
        global_max_dPk_phys = max(max_abs_dPk_phys_all)
        idx_global_max = max_abs_dPk_phys_all.index(global_max_dPk_phys)
        z_at_max_dPk_phys = z_at_max_dPk_phys_all[idx_global_max]
        k_at_max_dPk_phys = k_at_max_dPk_phys_all[idx_global_max]
        
        log_print(f"max |ΔP/P| (k ≤ 1 h/Mpc, all z) = {global_max_dPk_phys:.6e}")
        log_print(f"  at z = {z_at_max_dPk_phys:.1f}, k = {k_at_max_dPk_phys:.6e} h/Mpc")
        log_print()
        
        results['pk'] = {
            'max_dPk_phys': global_max_dPk_phys,
            'z_max_dPk_phys': z_at_max_dPk_phys,
            'k_max_dPk_phys': k_at_max_dPk_phys
        }
    else:
        log_print("ERROR: No valid P(k) data found")
        results['pk'] = {'max_dPk_phys': np.inf}
    
    # C_ℓ files (TT, EE, TE)
    max_abs_dTT = None
    ell_at_max_dTT = None
    max_abs_dEE = None
    ell_at_max_dEE = None
    
    # TT
    cl_tt_file = Path("R4b_outputs/R4b_deltaCl_TT.dat")
    if cl_tt_file.exists():
        try:
            data_tt = np.loadtxt(str(cl_tt_file), comments="#")
            # Expected columns: ell, Cl_LCDM, Cl_SRC, delta_rel
            ell_tt = data_tt[:, 0]
            delta_tt = data_tt[:, 3]
            
            # Physical window: ℓ <= 2000
            mask_phys_tt = ell_tt <= 2000
            if np.any(mask_phys_tt):
                abs_delta_tt_phys = np.abs(delta_tt[mask_phys_tt])
                ell_tt_phys = ell_tt[mask_phys_tt]
                max_abs_dTT = np.max(abs_delta_tt_phys)
                idx_max = np.argmax(abs_delta_tt_phys)
                ell_at_max_dTT = ell_tt_phys[idx_max]
                
                log_print(f"max |ΔC_ℓ^TT/C_ℓ^TT| (ℓ ≤ 2000) = {max_abs_dTT:.6e} at ℓ = {ell_at_max_dTT:.0f}")
        except Exception as e:
            log_print(f"WARNING: Failed to read C_ℓ^TT file: {e}")
    else:
        log_print(f"WARNING: Could not find C_ℓ^TT file: {cl_tt_file}")
    
    # EE
    cl_ee_file = Path("R4b_outputs/R4b_deltaCl_EE.dat")
    if cl_ee_file.exists():
        try:
            data_ee = np.loadtxt(str(cl_ee_file), comments="#")
            ell_ee = data_ee[:, 0]
            delta_ee = data_ee[:, 3]
            
            mask_phys_ee = ell_ee <= 2000
            if np.any(mask_phys_ee):
                abs_delta_ee_phys = np.abs(delta_ee[mask_phys_ee])
                ell_ee_phys = ell_ee[mask_phys_ee]
                max_abs_dEE = np.max(abs_delta_ee_phys)
                idx_max = np.argmax(abs_delta_ee_phys)
                ell_at_max_dEE = ell_ee_phys[idx_max]
                
                log_print(f"max |ΔC_ℓ^EE/C_ℓ^EE| (ℓ ≤ 2000) = {max_abs_dEE:.6e} at ℓ = {ell_at_max_dEE:.0f}")
        except Exception as e:
            log_print(f"WARNING: Failed to read C_ℓ^EE file: {e}")
    else:
        log_print(f"WARNING: Could not find C_ℓ^EE file: {cl_ee_file}")
    
    # TE (diagnostic only)
    cl_te_file = Path("R4b_outputs/R4b_deltaCl_TE.dat")
    if cl_te_file.exists():
        log_print("(TE spectrum available for diagnostics, not used in PASS/FAIL)")
    
    log_print()
    
    results['cl'] = {
        'max_dTT': max_abs_dTT,
        'ell_max_dTT': ell_at_max_dTT,
        'max_dEE': max_abs_dEE,
        'ell_max_dEE': ell_at_max_dEE
    }
    
    # =======================================================================
    # Section 4: Low-ℓ TT & ISW (R7)
    # =======================================================================
    log_print("="*70)
    log_print("4. Low-ℓ TT & ISW Amplitude Proxy (R7)")
    log_print("="*70)
    log_print()
    
    r7_file = Path("R7_outputs/R7_deltaCl_lowL.dat")
    if not r7_file.exists():
        log_print(f"ERROR: Could not find R7 file: {r7_file}")
        log_print("Please run ./run_R7_isw_diagnostics.sh first")
        log_file.close()
        sys.exit(1)
    
    try:
        data_r7 = np.loadtxt(str(r7_file), comments="#")
        # Expected columns: ell, Cl_lcdm, Cl_src, delta_rel
        ell_r7 = data_r7[:, 0]
        delta_r7 = data_r7[:, 3]
        
        # Window 2 <= ℓ <= 10
        mask_2_10 = (ell_r7 >= 2) & (ell_r7 <= 10)
        if np.any(mask_2_10):
            abs_delta_2_10 = np.abs(delta_r7[mask_2_10])
            ell_2_10 = ell_r7[mask_2_10]
            max_lowL_2_10 = np.max(abs_delta_2_10)
            idx_max = np.argmax(abs_delta_2_10)
            ell_at_max_2_10 = ell_2_10[idx_max]
        else:
            max_lowL_2_10 = None
            ell_at_max_2_10 = None
        
        # Window 2 <= ℓ <= 50 (for RMS)
        mask_2_50 = (ell_r7 >= 2) & (ell_r7 <= 50)
        if np.any(mask_2_50):
            abs_delta_2_50 = np.abs(delta_r7[mask_2_50])
            rms_lowL_2_50 = np.sqrt(np.mean(delta_r7[mask_2_50]**2))
        else:
            rms_lowL_2_50 = None
        
        # ISW amplitude proxy (2 <= ℓ <= 30)
        mask_2_30 = (ell_r7 >= 2) & (ell_r7 <= 30)
        if np.any(mask_2_30):
            ell_2_30 = ell_r7[mask_2_30]
            Cl_lcdm_2_30 = data_r7[mask_2_30, 1]
            Cl_src_2_30 = data_r7[mask_2_30, 2]
            
            # D_ℓ = ℓ(ℓ+1)C_ℓ / (2π)
            Dl_lcdm = ell_2_30 * (ell_2_30 + 1) * Cl_lcdm_2_30 / (2 * np.pi)
            Dl_src = ell_2_30 * (ell_2_30 + 1) * Cl_src_2_30 / (2 * np.pi)
            
            A_ISW = np.sum(Dl_src) / np.sum(Dl_lcdm)
        else:
            A_ISW = None
        
        if max_lowL_2_10 is not None:
            log_print(f"max |ΔC_ℓ/C_ℓ| (2≤ℓ≤10) = {max_lowL_2_10:.6e} at ℓ = {ell_at_max_2_10:.0f}")
        if rms_lowL_2_50 is not None:
            log_print(f"RMS |ΔC_ℓ/C_ℓ| (2≤ℓ≤50) = {rms_lowL_2_50:.6e}")
        if A_ISW is not None:
            log_print(f"A_ISW (2≤ℓ≤30) = {A_ISW:.6e}")
        log_print()
        
        results['lowL'] = {
            'max_2_10': max_lowL_2_10,
            'ell_max_2_10': ell_at_max_2_10,
            'rms_2_50': rms_lowL_2_50,
            'A_ISW': A_ISW
        }
        
    except Exception as e:
        log_print(f"ERROR: Failed to read R7 file: {e}")
        log_file.close()
        sys.exit(1)
    
    # =======================================================================
    # Combined PASS/FAIL Logic
    # =======================================================================
    log_print("="*70)
    log_print("COMBINED VALIDATION CRITERIA")
    log_print("="*70)
    log_print()
    
    pass_bg = (results['bg']['max_dH'] <= TOL_DH) and (results['bg']['max_dchi'] <= TOL_DCHI)
    
    pass_growth = (
        results['growth']['max_dD'] <= TOL_D and
        results['growth']['max_df'] <= TOL_F and
        results['growth']['max_dsigma8'] <= TOL_SIGMA8 and
        (results['growth']['max_dfsigma8_z15'] is not None) and
        (results['growth']['max_dfsigma8_z15'] <= TOL_FSIGMA8)
    )
    
    pass_pk = results['pk']['max_dPk_phys'] <= TOL_PK_PHYS
    
    pass_cl = (
        (results['cl']['max_dTT'] is not None) and (results['cl']['max_dTT'] <= TOL_CL_PHYS) and
        (results['cl']['max_dEE'] is not None) and (results['cl']['max_dEE'] <= TOL_CL_PHYS)
    )
    
    pass_lowL = (
        (results['lowL']['max_2_10'] is not None) and (results['lowL']['max_2_10'] <= TOL_LOWELL_MAX) and
        (results['lowL']['rms_2_50'] is not None) and (results['lowL']['rms_2_50'] <= TOL_LOWELL_RMS) and
        (results['lowL']['A_ISW'] is not None) and (np.abs(results['lowL']['A_ISW'] - 1.0) <= TOL_AISW_REL)
    )
    
    overall_pass = pass_bg and pass_growth and pass_pk and pass_cl and pass_lowL
    
    # =======================================================================
    # Summary
    # =======================================================================
    log_print("="*70)
    log_print("SUMMARY")
    log_print("="*70)
    log_print()
    
    # Background
    status_bg = '✓' if pass_bg else '✗'
    log_print(f"Background (R5): {status_bg}")
    log_print(f"  max |ΔH/H| = {results['bg']['max_dH']:.6e} (tolerance: {TOL_DH:.0e})")
    log_print(f"  max |Δχ/χ| = {results['bg']['max_dchi']:.6e} (tolerance: {TOL_DCHI:.0e})")
    log_print()
    
    # Growth
    status_growth = '✓' if pass_growth else '✗'
    log_print(f"Growth (R6): {status_growth}")
    log_print(f"  max |ΔD/D| = {results['growth']['max_dD']:.6e} (tolerance: {TOL_D:.0e})")
    log_print(f"  max |Δf/f| = {results['growth']['max_df']:.6e} (tolerance: {TOL_F:.0e})")
    log_print(f"  max |Δσ8/σ8| = {results['growth']['max_dsigma8']:.6e} (tolerance: {TOL_SIGMA8:.0e})")
    if results['growth']['max_dfsigma8_z15'] is not None:
        log_print(f"  max |Δ(fσ8)/(fσ8)| (z≤1.5) = {results['growth']['max_dfsigma8_z15']:.6e} (tolerance: {TOL_FSIGMA8:.0e})")
    log_print()
    
    # P(k)
    status_pk = '✓' if pass_pk else '✗'
    log_print(f"P(k) spectra (R4b): {status_pk}")
    log_print(f"  max |ΔP/P| (k≤1 h/Mpc) = {results['pk']['max_dPk_phys']:.6e} (tolerance: {TOL_PK_PHYS:.0e})")
    log_print()
    
    # C_ℓ
    status_cl = '✓' if pass_cl else '✗'
    log_print(f"C_ℓ spectra (R4b): {status_cl}")
    if results['cl']['max_dTT'] is not None:
        log_print(f"  max |ΔC_ℓ^TT/C_ℓ^TT| (ℓ≤2000) = {results['cl']['max_dTT']:.6e} (tolerance: {TOL_CL_PHYS:.0e})")
    if results['cl']['max_dEE'] is not None:
        log_print(f"  max |ΔC_ℓ^EE/C_ℓ^EE| (ℓ≤2000) = {results['cl']['max_dEE']:.6e} (tolerance: {TOL_CL_PHYS:.0e})")
    log_print()
    
    # Low-ℓ
    status_lowL = '✓' if pass_lowL else '✗'
    log_print(f"Low-ℓ TT & ISW (R7): {status_lowL}")
    if results['lowL']['max_2_10'] is not None:
        log_print(f"  max |ΔC_ℓ/C_ℓ| (2≤ℓ≤10) = {results['lowL']['max_2_10']:.6e} (tolerance: {TOL_LOWELL_MAX:.0e})")
    if results['lowL']['rms_2_50'] is not None:
        log_print(f"  RMS |ΔC_ℓ/C_ℓ| (2≤ℓ≤50) = {results['lowL']['rms_2_50']:.6e} (tolerance: {TOL_LOWELL_RMS:.0e})")
    if results['lowL']['A_ISW'] is not None:
        log_print(f"  |A_ISW - 1| = {np.abs(results['lowL']['A_ISW'] - 1.0):.6e} (tolerance: {TOL_AISW_REL:.0e})")
    log_print()
    
    # Result
    log_print("="*70)
    if overall_pass:
        log_print("RESULT: ✓ PASS - SRC global observables consistent with ΛCDM")
        log_print("        at O(1–2%) in the physical window")
    else:
        log_print("RESULT: ✗ FAIL - One or more global diagnostics exceed tolerances")
        failed = []
        if not pass_bg:
            failed.append("Background")
        if not pass_growth:
            failed.append("Growth")
        if not pass_pk:
            failed.append("P(k)")
        if not pass_cl:
            failed.append("C_ℓ")
        if not pass_lowL:
            failed.append("Low-ℓ TT & ISW")
        log_print(f"        Failed groups: {', '.join(failed)}")
    log_print("="*70)
    log_print()
    
    # =======================================================================
    # Export Combined Table
    # =======================================================================
    log_print("="*70)
    log_print("Exporting combined table...")
    log_print("="*70)
    log_print()
    
    output_dir = Path("R8_outputs")
    output_dir.mkdir(exist_ok=True)
    
    # Use R6 z grid as main grid (0, 0.5, 1, 2)
    z_main = results['growth']['z_array']
    
    # Interpolate R5 data onto R6 grid
    delta_H_interp = np.interp(z_main, results['bg']['z_array'], results['bg']['delta_H'])
    delta_chi_interp = np.interp(z_main, results['bg']['z_array'], results['bg']['delta_chi'])
    
    # Build combined table
    data_array = np.column_stack([
        z_main,
        delta_H_interp,
        delta_chi_interp,
        results['growth']['delta_D'],
        results['growth']['delta_f'],
        results['growth']['delta_sigma8'],
        results['growth']['delta_fsigma8']
    ])
    
    header = (
        "# z  delta_H_over_H  delta_chi_over_chi  delta_D_over_D  "
        "delta_f_over_f  delta_sigma8_over_sigma8  delta_fsigma8_over_fsigma8\n"
        "# Global observables summary (SRC vs LCDM)\n"
        "# All relative differences: (SRC - LCDM) / LCDM\n"
        "# Data aggregated from R5 (background), R6 (growth) validation blocks"
    )
    
    output_file = output_dir / "R8_global_observables.dat"
    np.savetxt(output_file, data_array, header=header, fmt="%.10e")
    log_print(f"Saved to: {output_file}")
    log_print()
    
    log_file.close()
    sys.exit(0 if overall_pass else 1)


if __name__ == "__main__":
    main()

