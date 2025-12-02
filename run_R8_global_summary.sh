#!/bin/bash
# R8 Global Observables Summary (no new CLASS runs)
# Run from: ~/dev/cosmo-class/

set -e

echo "=== R8 Global Observables Summary ==="
echo "Working directory: $(pwd)"
echo
echo "Reusing outputs from R4b, R5, R6, and R7 (no new CLASS runs)."
echo

python3 compare_src_R8_summary.py | tee R8_compare.log

echo
echo "=== R8 Summary Complete ==="
echo "Outputs:"
echo "  Log:    R8_compare.log"
echo "  Data:   R8_outputs/R8_global_observables.dat"

