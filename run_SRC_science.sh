#!/bin/bash
# Run SRC science simulation and generate plots
# 
# This script:
# 1. Runs LCDM simulation (use_src = no) → output/lcdm/
# 2. Runs SRC simulation (use_src = yes) → output/src/
# 3. Generates comparison plots using plot_SRC_results.py
#
# Usage: ./run_SRC_science.sh

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "=========================================================="
echo "SRC Science Run"
echo "=========================================================="
echo

# Check if CLASS is compiled
if [ ! -f "CLASS/class" ]; then
    echo "Error: CLASS executable not found. Please compile CLASS first:"
    echo "  cd CLASS && make clean && make -j4"
    exit 1
fi

# Create output directories
echo "Creating output directories..."
mkdir -p CLASS/output/lcdm
mkdir -p CLASS/output/src
mkdir -p plots

# Clean previous outputs
echo "Cleaning previous outputs..."
rm -f CLASS/output/lcdm/lcdm*.dat
rm -f CLASS/output/lcdm/lcdm*.ini
rm -f CLASS/output/src/src*.dat
rm -f CLASS/output/src/src*.ini

# Copy .ini files to CLASS directory
echo "Preparing input files..."
cp src_science_LCDM.ini CLASS/
cp src_science_SRC.ini CLASS/

# Run LCDM simulation (use_src = no)
echo
echo "=========================================================="
echo "Running LCDM simulation (use_src = no)..."
echo "=========================================================="
cd CLASS
./class src_science_LCDM.ini
cd ..

# Check if LCDM output was generated
if [ ! -f "CLASS/output/lcdm/lcdm_background.dat" ]; then
    echo "Error: CLASS did not generate LCDM output file"
    exit 1
fi

echo "✓ LCDM run completed"
echo

# Run SRC simulation (use_src = yes)
echo "=========================================================="
echo "Running SRC simulation (use_src = yes)..."
echo "=========================================================="
cd CLASS
./class src_science_SRC.ini
cd ..

# Check if SRC output was generated
if [ ! -f "CLASS/output/src/src_background.dat" ]; then
    echo "Error: CLASS did not generate SRC output file"
    exit 1
fi

echo "✓ SRC run completed"
echo

# Generate comparison plots
echo "=========================================================="
echo "Generating comparison plots..."
echo "=========================================================="
python3 plot_SRC_results.py \
    CLASS/output/lcdm/lcdm_background.dat \
    CLASS/output/src/src_background.dat \
    plots

echo
echo "=========================================================="
echo "SRC Science Run Complete"
echo "=========================================================="
echo
echo "Output files:"
echo "  LCDM:"
echo "    - CLASS/output/lcdm/lcdm_background.dat"
echo "    - CLASS/output/lcdm/lcdm_thermodynamics.dat"
echo "    - CLASS/output/lcdm/lcdm_primordial.dat"
echo "  SRC:"
echo "    - CLASS/output/src/src_background.dat"
echo "    - CLASS/output/src/src_thermodynamics.dat"
echo "    - CLASS/output/src/src_primordial.dat"
echo
if [ -d "plots" ] && [ -n "$(ls -A plots/*.png 2>/dev/null)" ]; then
    echo "Plots:"
    ls -1 plots/*.png 2>/dev/null | sed 's/^/  - /'
fi
echo
