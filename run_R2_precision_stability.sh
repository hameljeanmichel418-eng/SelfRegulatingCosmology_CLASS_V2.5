#!/bin/bash
# R2 Precision-Stability Validation: Compare base vs high-accuracy SRC runs
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

echo -e "${YELLOW}=== R2 Precision-Stability Validation ===${NC}"
echo "Working directory: $SCRIPT_DIR"
echo ""

# Check prerequisites
if [ ! -d "CLASS" ]; then
    echo -e "${RED}✗ Error: CLASS directory not found${NC}"
    exit 1
fi

if [ ! -f "CLASS/class" ]; then
    echo -e "${RED}✗ Error: CLASS/class executable not found${NC}"
    echo "   Please build CLASS first: cd CLASS && make -j4"
    exit 1
fi

# Step 1: Copy ini files to CLASS directory
echo -e "${YELLOW}[1/4] Setting up R2 ini files...${NC}"
cp src_R2_base.ini CLASS/
cp src_R2_hiacc.ini CLASS/
echo -e "${GREEN}✓ Copied ini files to CLASS/${NC}"

# Step 2: Create output directories
echo -e "${YELLOW}[2/4] Creating output directories...${NC}"
mkdir -p CLASS/output/R2_base
mkdir -p CLASS/output/R2_hiacc
echo -e "${GREEN}✓ Output directories created${NC}"

# Step 3: Run base precision
echo -e "${YELLOW}[3/4] Running base precision SRC...${NC}"
cd CLASS
./class src_R2_base.ini > ../R2_SRC_base.log 2>&1
if [ $? -eq 0 ]; then
    echo -e "${GREEN}✓ Base precision run completed${NC}"
else
    echo -e "${RED}✗ Base precision run failed (check R2_SRC_base.log)${NC}"
    exit 1
fi
cd ..

# Step 4: Run high-accuracy precision
echo -e "${YELLOW}[4/4] Running high-accuracy precision SRC...${NC}"
cd CLASS
./class src_R2_hiacc.ini > ../R2_SRC_hiacc.log 2>&1
if [ $? -eq 0 ]; then
    echo -e "${GREEN}✓ High-accuracy precision run completed${NC}"
else
    echo -e "${RED}✗ High-accuracy precision run failed (check R2_SRC_hiacc.log)${NC}"
    exit 1
fi
cd ..

echo ""
echo -e "${GREEN}=== R2 Runs Complete ===${NC}"
echo ""
echo "Output files generated:"
echo "  Base:    CLASS/output/R2_base/R2_base_*.dat"
echo "  HiAcc:   CLASS/output/R2_hiacc/R2_hiacc_*.dat"
echo ""
echo "Logs saved:"
echo "  R2_SRC_base.log"
echo "  R2_SRC_hiacc.log"
echo ""
echo "Next step: python3 compare_src_R2_precision.py"
echo ""

