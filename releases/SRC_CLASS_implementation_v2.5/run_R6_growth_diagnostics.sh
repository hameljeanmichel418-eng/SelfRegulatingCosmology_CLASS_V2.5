#!/bin/bash
# R6 Linear Growth & fσ8 Diagnostics: Using R4 outputs
# Run from: ~/dev/cosmo-class/

set -e

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo -e "${YELLOW}=== R6 Linear Growth & fσ8 Diagnostics ===${NC}"
echo "Working directory: $SCRIPT_DIR"
echo ""

# Check if R4 outputs exist
R4_LCDM_DIR="CLASS_clean/output/R4_lcdm"
R4_SRC_DIR="CLASS/output/R4_src"

if [ ! -d "$R4_LCDM_DIR" ] || [ ! -d "$R4_SRC_DIR" ]; then
    echo -e "${YELLOW}R4 outputs not found. Running R4 first...${NC}"
    echo ""
    ./run_R4_perturbation_comparison.sh
    if [ $? -ne 0 ]; then
        echo -e "${RED}✗ R4 run failed${NC}"
        exit 1
    fi
    echo ""
fi

# Check for required files
REQUIRED_FILES=(
    "$R4_LCDM_DIR/R4_lcdm_z1_pk.dat"
    "$R4_LCDM_DIR/R4_lcdm_z2_pk.dat"
    "$R4_LCDM_DIR/R4_lcdm_z3_pk.dat"
    "$R4_LCDM_DIR/R4_lcdm_z4_pk.dat"
    "$R4_SRC_DIR/R4_src_z1_pk.dat"
    "$R4_SRC_DIR/R4_src_z2_pk.dat"
    "$R4_SRC_DIR/R4_src_z3_pk.dat"
    "$R4_SRC_DIR/R4_src_z4_pk.dat"
)

MISSING_FILES=()
for file in "${REQUIRED_FILES[@]}"; do
    if [ ! -f "$file" ]; then
        MISSING_FILES+=("$file")
    fi
done

if [ ${#MISSING_FILES[@]} -gt 0 ]; then
    echo -e "${YELLOW}Some R4 output files are missing. Running R4 first...${NC}"
    echo ""
    ./run_R4_perturbation_comparison.sh
    if [ $? -ne 0 ]; then
        echo -e "${RED}✗ R4 run failed${NC}"
        exit 1
    fi
    echo ""
fi

# Run R6 comparison
echo -e "${YELLOW}Running R6 growth diagnostics...${NC}"
python3 compare_src_R6_growth.py | tee R6_compare.log

if [ $? -eq 0 ]; then
    echo ""
    echo -e "${GREEN}=== R6 Diagnostics Complete ===${NC}"
    echo ""
    echo "Output files:"
    echo "  Log: R6_compare.log"
    echo "  Data: R6_outputs/R6_growth_table.dat"
    echo ""
else
    echo ""
    echo -e "${RED}✗ R6 diagnostics failed (check R6_compare.log)${NC}"
    exit 1
fi

