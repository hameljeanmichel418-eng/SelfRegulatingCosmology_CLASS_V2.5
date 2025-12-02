#!/bin/bash
# Verification script for SRC implementation
# This script compiles both clean and modified CLASS versions
# Run from: ~/dev/cosmo-class/

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Get script directory (~/dev/cosmo-class/)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo -e "${YELLOW}=== SRC Implementation Verification ===${NC}"
echo "Working directory: $SCRIPT_DIR"
echo ""

# Check that modified CLASS directory exists
if [ ! -d "CLASS" ]; then
    echo -e "${RED}✗ Error: Modified CLASS not found at $SCRIPT_DIR/CLASS/${NC}"
    exit 1
fi

# Step 1: Extract clean CLASS if needed
echo -e "${YELLOW}[1/3] Setting up clean CLASS...${NC}"
if [ -d "CLASS_clean" ]; then
    echo "  CLASS_clean/ already exists, skipping extraction"
    echo -e "${GREEN}✓ Using existing CLASS_clean/${NC}"
elif [ -f "CLASS_clean_v3.3.3.tgz" ]; then
    echo "  Extracting CLASS_clean_v3.3.3.tgz to CLASS_clean/..."
    mkdir -p CLASS_clean
    tar -xzf CLASS_clean_v3.3.3.tgz -C CLASS_clean --strip-components=1
    if [ $? -eq 0 ]; then
        echo -e "${GREEN}✓ Clean CLASS extracted to CLASS_clean/${NC}"
    else
        echo -e "${RED}✗ Extraction failed${NC}"
        exit 1
    fi
else
    echo -e "${RED}✗ Error: CLASS_clean_v3.3.3.tgz not found in $SCRIPT_DIR${NC}"
    exit 1
fi

# Step 2: Compile clean CLASS
echo -e "${YELLOW}[2/3] Compiling clean CLASS...${NC}"
cd CLASS_clean
make clean > /dev/null 2>&1 || true
make -j4 > compile_clean.log 2>&1
if [ $? -eq 0 ]; then
    echo -e "${GREEN}✓ Clean CLASS compiled successfully${NC}"
else
    echo -e "${RED}✗ Clean CLASS compilation failed (check CLASS_clean/compile_clean.log)${NC}"
    exit 1
fi
cd ..

# Step 3: Compile modified CLASS
echo -e "${YELLOW}[3/3] Compiling modified CLASS...${NC}"
cd CLASS
make clean > /dev/null 2>&1 || true
make -j4 > compile_mod.log 2>&1
if [ $? -eq 0 ]; then
    echo -e "${GREEN}✓ Modified CLASS compiled successfully${NC}"
else
    echo -e "${RED}✗ Modified CLASS compilation failed (check CLASS/compile_mod.log)${NC}"
    exit 1
fi
cd ..

# Step 4: Create output directories
echo -e "${YELLOW}Creating output directories...${NC}"
mkdir -p CLASS_clean/output
mkdir -p CLASS/output
echo -e "${GREEN}✓ Output directories created${NC}"

echo ""
echo -e "${GREEN}=== Compilation Complete ===${NC}"
echo ""
echo "Next steps:"
echo "  1. Run: ./run_comparison.sh"
echo "  2. Compare results: python3 compare_outputs.py"

