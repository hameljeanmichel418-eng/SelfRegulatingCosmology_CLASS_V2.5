#!/bin/bash
# R7 ISW & low-ℓ TT Diagnostics
# Run from: ~/dev/cosmo-class/

set -e

echo "=== R7 ISW & low-ℓ TT Diagnostics ==="
echo "Working directory: $(pwd)"
echo

# R7 uses R4/R4b outputs; no new CLASS runs.
echo "Using R4/R4b CMB spectra (no rebuild)."
echo

python3 compare_src_R7_isw.py | tee R7_compare.log

echo
echo "=== R7 Diagnostics Complete ==="
echo
echo "Outputs:"
echo "  Log:    R7_compare.log"
echo "  Data:   R7_outputs/R7_deltaCl_lowL.dat"

