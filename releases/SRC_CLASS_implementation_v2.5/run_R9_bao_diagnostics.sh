#!/bin/bash
# R9 BAO-scale Diagnostics (SRC vs LCDM)
# Run from: ~/dev/cosmo-class

set -e

echo "=== R9 BAO-scale Diagnostics ==="
echo "Working directory: $(pwd)"
echo
echo "Reusing outputs from R5 (no new CLASS runs)."
echo

python3 compare_src_R9_bao.py | tee R9_compare.log

echo
echo "=== R9 Diagnostics Complete ==="
echo
echo "Outputs:"
echo "  Log:    R9_compare.log"
echo "  Data:   R9_outputs/R9_bao_results.dat"

