#!/usr/bin/env python3
"""
Compare CLASS output files between clean and modified versions.
This script verifies that modified CLASS with use_src=no produces
identical results to clean CLASS v3.3.3.
Run from: ~/dev/cosmo-class/
"""

import numpy as np
import sys
import os
from pathlib import Path

# Tolerance for numerical comparison
TOLERANCE = 1e-10

def read_class_file(filename):
    """
    Read a CLASS output file, skipping header lines.
    Returns data array and number of header lines.
    """
    data_lines = []
    header_lines = 0
    in_header = True
    
    with open(filename, 'r') as f:
        for line in f:
            # Skip empty lines
            if not line.strip():
                if in_header:
                    header_lines += 1
                continue
            
            # Check if line is a header (starts with # or contains non-numeric data)
            if in_header:
                if line.strip().startswith('#'):
                    header_lines += 1
                    continue
                # Try to parse as numbers - if it fails, it's still header
                try:
                    # Try to split and convert to float
                    parts = line.split()
                    [float(x) for x in parts]
                    in_header = False
                    data_lines.append(line)
                except (ValueError, IndexError):
                    header_lines += 1
                    continue
            else:
                data_lines.append(line)
    
    if not data_lines:
        return None, header_lines
    
    # Parse data
    try:
        data = np.loadtxt(data_lines, dtype=float)
        # Handle 1D case
        if data.ndim == 1:
            data = data.reshape(-1, 1)
        return data, header_lines
    except ValueError as e:
        print(f"Warning: Could not parse {filename} as numeric data")
        return None, header_lines

def compare_files(file_clean, file_mod, filename):
    """
    Compare two CLASS output files.
    Returns (match, message)
    """
    if not os.path.exists(file_clean):
        return False, f"Clean file not found: {file_clean}"
    
    if not os.path.exists(file_mod):
        return False, f"Modified file not found: {file_mod}"
    
    data_clean, hdr_clean = read_class_file(file_clean)
    data_mod, hdr_mod = read_class_file(file_mod)
    
    if data_clean is None:
        return False, f"Could not read clean file: {file_clean}"
    
    if data_mod is None:
        return False, f"Could not read modified file: {file_mod}"
    
    # Check shape
    if data_clean.shape != data_mod.shape:
        return False, (f"Shape mismatch: clean {data_clean.shape} vs "
                      f"modified {data_mod.shape}")
    
    # Compare numerical values
    diff = np.abs(data_clean - data_mod)
    max_diff = np.max(diff)
    max_diff_idx = np.unravel_index(np.argmax(diff), diff.shape)
    
    if max_diff > TOLERANCE:
        # Find the first row with a mismatch
        row_mismatch = None
        for i in range(len(data_clean)):
            row_diff = np.max(np.abs(data_clean[i] - data_mod[i]))
            if row_diff > TOLERANCE:
                row_mismatch = i
                break
        
        if row_mismatch is not None:
            return False, (f"Line {row_mismatch + 1 + hdr_clean}: "
                          f"max difference = {max_diff:.2e} at index {max_diff_idx}, "
                          f"clean = {data_clean[max_diff_idx]}, "
                          f"mod = {data_mod[max_diff_idx]}")
        else:
            return False, f"Max difference = {max_diff:.2e} at index {max_diff_idx}"
    
    return True, f"OK (max diff = {max_diff:.2e}, header lines: clean={hdr_clean}, mod={hdr_mod})"

def main():
    """Main comparison function"""
    
    # Get script directory (~/dev/cosmo-class/)
    script_dir = Path(__file__).parent.absolute()
    os.chdir(script_dir)
    
    # Define paths relative to script directory
    base_clean = Path("CLASS_clean/output")
    base_mod = Path("CLASS/output")
    
    if not base_clean.exists():
        print("Error: CLASS_clean/output directory not found")
        print(f"Current directory: {script_dir}")
        print("Please run ./run_comparison.sh first")
        sys.exit(1)
    
    if not base_mod.exists():
        print("Error: CLASS/output directory not found")
        print(f"Current directory: {script_dir}")
        print("Please run ./run_comparison.sh first")
        sys.exit(1)
    
    # Find all output files
    pattern = "src_test_LCDM_*.dat"
    files_clean = list(base_clean.glob(pattern))
    
    if not files_clean:
        print(f"Error: No output files found matching {pattern} in {base_clean}")
        print("Please run ./run_comparison.sh first")
        sys.exit(1)
    
    print("=" * 70)
    print("CLASS Output Comparison (use_src = no)")
    print("=" * 70)
    print(f"Working directory: {script_dir}")
    print(f"Tolerance: {TOLERANCE}")
    print()
    
    all_match = True
    results = []
    
    # Compare each file
    for file_clean in sorted(files_clean):
        filename = file_clean.name
        file_mod = base_mod / filename
        
        match, message = compare_files(str(file_clean), str(file_mod), filename)
        results.append((filename, match, message))
        
        if not match:
            all_match = False
    
    # Print results
    print("File Comparison Results:")
    print("-" * 70)
    for filename, match, message in results:
        status = "✓" if match else "✗"
        print(f"{status} {filename:40s} {message}")
    print("-" * 70)
    print()
    
    # Final verdict
    if all_match:
        print("=" * 70)
        print("OK: modified CLASS is identical to clean CLASS when use_src = no")
        print("=" * 70)
        return 0
    else:
        print("=" * 70)
        print("✗ FAILURE: Differences detected between clean and modified CLASS")
        print("=" * 70)
        return 1

if __name__ == "__main__":
    sys.exit(main())

