#!/bin/bash
# R4 Perturbation Spectra Validation: Compare LCDM vs SRC
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

echo -e "${YELLOW}=== R4 Perturbation Spectra Validation: LCDM vs SRC ===${NC}"
echo "Working directory: $SCRIPT_DIR"
echo ""

# Check prerequisites
if [ ! -d "CLASS_clean" ]; then
    echo -e "${RED}✗ Error: CLASS_clean directory not found${NC}"
    exit 1
fi

if [ ! -d "CLASS" ]; then
    echo -e "${RED}✗ Error: CLASS directory not found${NC}"
    exit 1
fi

# Step 1: Build CLASS_clean if needed
echo -e "${YELLOW}[1/5] Building CLASS_clean...${NC}"
cd CLASS_clean
if [ ! -f "class" ] || [ "../src_R4_LCDM.ini" -nt "class" ]; then
    make clean > /dev/null 2>&1
    make -j4 > compile_clean.log 2>&1
    if [ $? -eq 0 ]; then
        echo -e "${GREEN}✓ CLASS_clean built successfully${NC}"
    else
        echo -e "${RED}✗ CLASS_clean build failed (check CLASS_clean/compile_clean.log)${NC}"
        exit 1
    fi
else
    echo -e "${GREEN}✓ CLASS_clean already built${NC}"
fi
cd ..

# Step 2: Build CLASS if needed
echo -e "${YELLOW}[2/5] Building CLASS (SRC-modified)...${NC}"
cd CLASS
if [ ! -f "class" ] || [ "../src_R4_SRC.ini" -nt "class" ]; then
    make clean > /dev/null 2>&1
    make -j4 > compile_src.log 2>&1
    if [ $? -eq 0 ]; then
        echo -e "${GREEN}✓ CLASS built successfully${NC}"
    else
        echo -e "${RED}✗ CLASS build failed (check CLASS/compile_src.log)${NC}"
        exit 1
    fi
else
    echo -e "${GREEN}✓ CLASS already built${NC}"
fi
cd ..

# Step 3: Create output directories
echo -e "${YELLOW}[3/5] Creating output directories...${NC}"
mkdir -p CLASS_clean/output/R4_lcdm
mkdir -p CLASS/output/R4_src
echo -e "${GREEN}✓ Output directories created${NC}"

# Step 4: Run LCDM baseline
echo -e "${YELLOW}[4/5] Running LCDM baseline (CLASS_clean)...${NC}"
cp src_R4_LCDM.ini CLASS_clean/
cd CLASS_clean
./class src_R4_LCDM.ini > ../R4_LCDM.log 2>&1
if [ $? -eq 0 ]; then
    echo -e "${GREEN}✓ LCDM run completed${NC}"
else
    echo -e "${RED}✗ LCDM run failed (check R4_LCDM.log)${NC}"
    exit 1
fi
cd ..

# Step 5: Run SRC science run
echo -e "${YELLOW}[5/5] Running SRC science run (CLASS)...${NC}"
cp src_R4_SRC.ini CLASS/
cd CLASS
./class src_R4_SRC.ini > ../R4_SRC.log 2>&1
if [ $? -eq 0 ]; then
    echo -e "${GREEN}✓ SRC run completed${NC}"
else
    echo -e "${RED}✗ SRC run failed (check R4_SRC.log)${NC}"
    exit 1
fi
cd ..

echo ""
echo -e "${GREEN}=== R4 Runs Complete ===${NC}"
echo ""
echo "Output files generated:"
echo "  LCDM:  CLASS_clean/output/R4_lcdm/R4_lcdm_*.dat"
echo "  SRC:   CLASS/output/R4_src/R4_src_*.dat"
echo ""
echo "Logs saved:"
echo "  R4_LCDM.log"
echo "  R4_SRC.log"
echo ""
echo "Next step: python3 compare_src_R4_perturbations.py | tee R4_compare.log"
echo ""

