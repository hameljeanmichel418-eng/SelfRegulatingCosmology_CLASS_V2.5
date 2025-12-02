# Self-Regulating Cosmology (SRC) Implementation in CLASS

## Project Overview

This repository contains the complete implementation and validation of the Self-Regulating Cosmology (SRC) model within the CLASS (Cosmic Linear Anisotropy Solving System) framework. The SRC model introduces a late-time energy transfer mechanism from primordial black hole (PBH) evaporation into a light scalar field, resulting in a small, controlled modulation of the Hubble expansion history H(z).

**Version:** 2.5 (Final)  
**Author:** Jm Hamel  
**Year:** 2025

## Project Structure

This project spans three volumes:

- **Volume 1:** Theoretical foundations and conceptual consistency checks (A1–A10)
- **Volume 2:** CLASS implementation of the SRC background modification
- **Volume 2.5:** Complete numerical validation pipeline (R1–R9)

## Model Description

The SRC model adds a late-time energy transfer term Q(a) to the fluid component:

```
Q(a) = gamma0 * exp(-0.5 * (ln(a/a_peak) / width)^2)
```

This contributes to the background evolution of the fluid:

```
dot(rho_fld) = -3H(1+w) rho_fld + Q(a)
```

The implementation is **background-only** in Volume 2.5; no perturbation-level coupling is included.

## Validation Pipeline

### Section A — Conceptual Baseline (Volume 1)

The A1–A10 checks establish the physical and mathematical foundations:

- **A1–A4:** Physical consistency, thermodynamics, energy transfer, stability
- **A5–A7:** Scaling predictions, fine-tuning avoidance, naturalness
- **A8–A10:** Causality, conservation laws, observational signatures

All A1–A10 checks **PASS** (see `VERIFICATION_SUMMARY.md` for details).

### Section R — Numerical Validation Pipeline (Volume 2.5)

The R1–R9 suite provides comprehensive CLASS-based validation:

- **R1:** Baseline numerical stability (CLASS compilation, background monotonicity, closure)
- **R2:** Early ΔH/H & Δχ/χ inspection
- **R3:** Raw ΔP(k) inspection
- **R4:** Preliminary CMB spectra inspection
- **R4b:** Strict ΔP(k) and ΔC_ℓ validation with physical windows
- **R5:** Late-time background & distances (H(z), χ(z), D_L(z))
- **R6:** Linear growth, σ8, fσ8 diagnostics
- **R7:** Low-ℓ TT & ISW test
- **R8:** Global observable summary (consolidated from R4b, R5, R6, R7)
- **R9:** BAO-scale distances (r_s, D_H, D_M, D_V)

All R1–R9 validation blocks **PASS** (see `VERIFICATION_SUMMARY.md` for numerical results).

## Compilation Instructions

### Prerequisites

- C compiler (gcc recommended)
- Make
- Python 3 with NumPy (for analysis scripts)
- WSL/Linux environment (for shell scripts)

### Building CLASS

1. **Clean build (recommended for first-time setup):**
   ```bash
   cd CLASS
   make clean
   make -j4
   ```

2. **Build clean upstream CLASS (for R1 validation):**
   ```bash
   cd CLASS_clean
   make clean
   make -j4
   ```

3. **Verify compilation:**
   ```bash
   ./CLASS/class --version
   ./CLASS_clean/class --version
   ```

## Running Validation Tests

### Quick Start

All validation scripts are located in the project root (`~/dev/cosmo-class/`).

### Individual R-Tests

**R1 – Baseline numerical stability:**
```bash
./run_R1_validation.sh
```

**R2 – Precision-stability validation:**
```bash
./run_R2_precision_stability.sh
```

**R3 – Background physics validation:**
```bash
./run_R3_background_comparison.sh
```

**R4 – Perturbation spectra validation:**
```bash
./run_R4_perturbation_comparison.sh
```

**R4b – Extended spectra diagnostics:**
```bash
./run_R4b_spectra_diagnostics.sh
```

**R5 – Late-time H(z) & distance diagnostics:**
```bash
./run_R5_late_time_observables.sh
```

**R6 – Linear growth & fσ8 diagnostics:**
```bash
./run_R6_growth_diagnostics.sh
```

**R7 – ISW & low-ℓ CMB diagnostics:**
```bash
./run_R7_isw_diagnostics.sh
```

**R8 – Global observables summary:**
```bash
./run_R8_global_summary.sh
```

**R9 – BAO-scale diagnostics:**
```bash
./run_R9_bao_diagnostics.sh
```

### Running All Tests

To run the complete validation suite sequentially:

```bash
./run_R1_validation.sh
./run_R2_precision_stability.sh
./run_R3_background_comparison.sh
./run_R4_perturbation_comparison.sh
./run_R4b_spectra_diagnostics.sh
./run_R5_late_time_observables.sh
./run_R6_growth_diagnostics.sh
./run_R7_isw_diagnostics.sh
./run_R8_global_summary.sh
./run_R9_bao_diagnostics.sh
```

**Note:** R4b, R6, R7, R8, and R9 reuse outputs from earlier tests (R4, R5) and do not require new CLASS runs.

## Output Files and Paths

### CLASS Outputs

- **R1:** `CLASS_clean/output/`, `CLASS/output/src_gamma0/`
- **R2:** `CLASS/output/R2_base/`, `CLASS/output/R2_hiacc/`
- **R3:** `CLASS_clean/output/R3_lcdm/`, `CLASS/output/R3_src/`
- **R4:** `CLASS_clean/output/R4_lcdm/`, `CLASS/output/R4_src/`
- **R5:** `CLASS_clean/output/R5_lcdm/`, `CLASS/output/R5_src/`

### Analysis Outputs

- **R4b:** `R4b_outputs/` (ΔP/P and ΔC_ℓ data files)
- **R5:** `R5_outputs/R5_delta_background.dat`
- **R6:** `R6_outputs/R6_growth_table.dat`
- **R7:** `R7_outputs/R7_deltaCl_lowL.dat`
- **R8:** `R8_outputs/R8_global_observables.dat`
- **R9:** `R9_outputs/R9_bao_results.dat`

### Log Files

Each validation script produces a log file:
- `R1_compare.log`, `R2_compare.log`, `R3_compare.log`, etc.

### Comparison Scripts

All comparison scripts are in the project root:
- `compare_lcdm_clean_vs_src_gamma0.py` (R1)
- `compare_src_R2_precision.py` (R2)
- `compare_src_R3_background.py` (R3)
- `compare_src_R4_perturbations.py` (R4)
- `compare_src_R4b_spectra.py` (R4b)
- `compare_src_R5_late_time.py` (R5)
- `compare_src_R6_growth.py` (R6)
- `compare_src_R7_isw.py` (R7)
- `compare_src_R8_summary.py` (R8)
- `compare_src_R9_bao.py` (R9)

## Key Files

### Documentation

- **`VERIFICATION_SUMMARY.md`** — Complete verification record with all A1–A10 and R1–R9 results
- **`README_final.md`** — This file

### Configuration Files

- **`src_science_LCDM.ini`** — Pure ΛCDM reference configuration
- **`src_science_SRC.ini`** — SRC science configuration (src_gamma0 = 0.02)
- **`src_science_SRC_gamma0.ini`** — SRC with gamma0 = 0.0 (for R1 regression test)
- **`src_R2_base.ini`, `src_R2_hiacc.ini`** — R2 precision tests
- **`src_R3_LCDM.ini`, `src_R3_SRC.ini`** — R3 background validation
- **`src_R4_LCDM.ini`, `src_R4_SRC.ini`** — R4 perturbation validation
- **`src_R5_LCDM.ini`, `src_R5_SRC.ini`** — R5 late-time observables

### CLASS Source Modifications

The SRC implementation modifies:
- **`CLASS/source/background.c`** — Adds Q(a) term to fluid continuity equation
- **`CLASS/source/background.h`** — SRC parameter declarations
- **`CLASS/source/input.c`** — SRC parameter parsing from .ini files

**Important:** The implementation is background-only; no perturbation modules are modified.

## Validation Results Summary

All validation tests pass. Key results:

- **Background:** max |ΔH/H| = 2.04×10⁻³, max |Δχ/χ| = 3.99×10⁻⁵
- **Growth:** max |ΔD/D| = 4.49×10⁻⁴, max |Δf/f| = 1.68×10⁻³, max |Δσ8/σ8| = 4.47×10⁻⁴
- **Spectra:** max |ΔP/P| = 8.99×10⁻⁴, max |ΔC_ℓ/C_ℓ| = 2.17×10⁻³
- **BAO:** max |ΔD_V/D_V| = 5.72×10⁻⁴, Δr_s/r_s = 0.00

All differences are at the 10⁻³–10⁻² level, consistent with the predicted Volume 1 theoretical window.

See `VERIFICATION_SUMMARY.md` for complete numerical results and validation criteria.

## Notes for Zenodo/Figshare Upload

### Recommended Archive Contents

- All CLASS source modifications (`CLASS/source/`, `CLASS/include/`)
- All validation scripts (`compare_src_R*.py`, `run_R*.sh`)
- All configuration files (`src_*.ini`)
- All analysis outputs (`R*_outputs/`, log files)
- Documentation (`VERIFICATION_SUMMARY.md`, `README_final.md`, `LICENSE`)
- Clean upstream CLASS reference (`CLASS_clean/`) for reproducibility

### Archive Metadata

- **Title:** Self-Regulating Cosmology (SRC) Implementation in CLASS v2.5
- **Description:** Complete implementation and validation of late-time PBH→φ energy transfer model
- **Keywords:** cosmology, CLASS, primordial black holes, scalar fields, late-time physics
- **Version:** 2.5 (Final)
- **License:** MIT

### Citation

If you use this code, please cite:
- The CLASS codebase (Lesgourgues 2011)
- This SRC implementation (Hamel 2025)
- The verification summary document

## Contact and Support

For questions or issues related to this implementation, please refer to the verification summary and log files for diagnostic information.

## License

This project is licensed under the MIT License. See `LICENSE` file for details.

---

**Last Updated:** 2025  
**Verification Status:** All A1–A10 and R1–R9 tests PASS  
**Documentation:** See `VERIFICATION_SUMMARY.md` for complete validation record

