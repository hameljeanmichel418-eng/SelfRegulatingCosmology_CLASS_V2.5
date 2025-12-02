#!/bin/bash
# R1 Validation: Compare CLASS_clean vs CLASS (SRC gamma0=0)
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

echo -e "${YELLOW}=== R1 Validation: CLASS_clean vs CLASS (SRC gamma0=0) ===${NC}"
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

# Step 1: Build CLASS_clean
echo -e "${YELLOW}[1/6] Building CLASS_clean...${NC}"
cd CLASS_clean
make clean > /dev/null 2>&1
make -j4 > compile_clean.log 2>&1
if [ $? -eq 0 ]; then
    echo -e "${GREEN}✓ CLASS_clean built successfully${NC}"
else
    echo -e "${RED}✗ CLASS_clean build failed (check CLASS_clean/compile_clean.log)${NC}"
    exit 1
fi
cd ..

# Step 2: Copy LCDM ini to CLASS_clean if needed
echo -e "${YELLOW}[2/6] Setting up CLASS_clean LCDM run...${NC}"
if [ ! -f "CLASS_clean/src_science_LCDM.ini" ]; then
    cp src_science_LCDM.ini CLASS_clean/
    echo -e "${GREEN}✓ Copied src_science_LCDM.ini to CLASS_clean/${NC}"
else
    echo -e "${GREEN}✓ src_science_LCDM.ini already exists in CLASS_clean/${NC}"
fi

# Step 3: Run CLASS_clean LCDM
echo -e "${YELLOW}[3/6] Running CLASS_clean LCDM...${NC}"
cd CLASS_clean
mkdir -p output/lcdm
./class src_science_LCDM.ini > run_clean.log 2>&1
if [ $? -eq 0 ]; then
    echo -e "${GREEN}✓ CLASS_clean LCDM run completed${NC}"
else
    echo -e "${RED}✗ CLASS_clean run failed (check CLASS_clean/run_clean.log)${NC}"
    exit 1
fi
cd ..

# Step 4: Build CLASS (SRC-modified)
echo -e "${YELLOW}[4/6] Building CLASS (SRC-modified)...${NC}"
cd CLASS
make clean > /dev/null 2>&1
make -j4 > compile_mod.log 2>&1
if [ $? -eq 0 ]; then
    echo -e "${GREEN}✓ CLASS built successfully${NC}"
else
    echo -e "${RED}✗ CLASS build failed (check CLASS/compile_mod.log)${NC}"
    exit 1
fi
cd ..

# Step 5: Copy SRC gamma0=0 ini to CLASS
echo -e "${YELLOW}[5/6] Setting up CLASS SRC gamma0=0 run...${NC}"
if [ ! -f "CLASS/src_science_SRC_gamma0.ini" ]; then
    cp src_science_SRC_gamma0.ini CLASS/
    echo -e "${GREEN}✓ Copied src_science_SRC_gamma0.ini to CLASS/${NC}"
else
    echo -e "${GREEN}✓ src_science_SRC_gamma0.ini already exists in CLASS/${NC}"
fi

# Step 6: Run CLASS SRC gamma0=0
echo -e "${YELLOW}[6/6] Running CLASS SRC gamma0=0...${NC}"
cd CLASS
mkdir -p output/src_gamma0
./class src_science_SRC_gamma0.ini > run_src_gamma0.log 2>&1
if [ $? -eq 0 ]; then
    echo -e "${GREEN}✓ CLASS SRC gamma0=0 run completed${NC}"
else
    echo -e "${RED}✗ CLASS run failed (check CLASS/run_src_gamma0.log)${NC}"
    exit 1
fi
cd ..

echo ""
echo -e "${GREEN}=== Builds and Runs Complete ===${NC}"
echo ""
echo "Next step: python3 compare_lcdm_clean_vs_src_gamma0.py"
echo ""

