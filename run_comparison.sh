#!/bin/bash
# Run both CLASS versions with the same .ini file and compare outputs
# Run from: ~/dev/cosmo-class/

set -e

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

# Get script directory (~/dev/cosmo-class/)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo -e "${YELLOW}=== Running CLASS Comparison ===${NC}"
echo "Working directory: $SCRIPT_DIR"
echo ""

# Check if .ini file exists
if [ ! -f "src_test_LCDM.ini" ]; then
    echo -e "${RED}✗ Error: src_test_LCDM.ini not found in $SCRIPT_DIR${NC}"
    exit 1
fi

# Check if both CLASS directories exist
if [ ! -d "CLASS_clean" ]; then
    echo -e "${RED}✗ Error: CLASS_clean directory not found${NC}"
    echo "   Please ensure CLASS_clean/ exists"
    exit 1
fi

if [ ! -d "CLASS" ]; then
    echo -e "${RED}✗ Error: CLASS directory not found${NC}"
    echo "   Please ensure CLASS/ exists"
    exit 1
fi

# Check if both executables exist
if [ ! -f "CLASS_clean/class" ]; then
    echo -e "${RED}✗ Error: CLASS_clean/class executable not found${NC}"
    echo "   Please run ./verify_src_implementation.sh first"
    exit 1
fi

if [ ! -f "CLASS/class" ]; then
    echo -e "${RED}✗ Error: CLASS/class executable not found${NC}"
    echo "   Please run ./verify_src_implementation.sh first"
    exit 1
fi

# Clean any existing output files with the test root to prevent numbering
echo -e "${YELLOW}[1/4] Cleaning existing test output files...${NC}"
rm -f CLASS_clean/output/src_test_LCDM*.dat CLASS_clean/output/src_test_LCDM*.ini
rm -f CLASS/output/src_test_LCDM*.dat CLASS/output/src_test_LCDM*.ini
echo -e "${GREEN}✓ Existing test outputs cleaned${NC}"

# Copy .ini to both directories
echo -e "${YELLOW}[2/4] Copying test .ini file...${NC}"
cp src_test_LCDM.ini CLASS_clean/
cp src_test_LCDM.ini CLASS/
echo -e "${GREEN}✓ .ini files copied${NC}"

# Run clean CLASS
echo -e "${YELLOW}[3/4] Running clean CLASS...${NC}"
cd CLASS_clean
./class src_test_LCDM.ini > run_clean.log 2>&1
if [ $? -eq 0 ]; then
    echo -e "${GREEN}✓ Clean CLASS completed successfully${NC}"
else
    echo -e "${RED}✗ Clean CLASS failed (check CLASS_clean/run_clean.log)${NC}"
    exit 1
fi
cd ..

# Run modified CLASS
echo -e "${YELLOW}[4/4] Running modified CLASS...${NC}"
cd CLASS
./class src_test_LCDM.ini > run_mod.log 2>&1
if [ $? -eq 0 ]; then
    echo -e "${GREEN}✓ Modified CLASS completed successfully${NC}"
else
    echo -e "${RED}✗ Modified CLASS failed (check CLASS/run_mod.log)${NC}"
    exit 1
fi
cd ..

echo ""
echo -e "${GREEN}=== Runs Complete ===${NC}"
echo ""
echo "Output files generated:"
echo "  Clean:    CLASS_clean/output/src_test_LCDM_*.dat"
echo "  Modified: CLASS/output/src_test_LCDM_*.dat"
echo ""
echo "Next step: python3 compare_outputs.py"

