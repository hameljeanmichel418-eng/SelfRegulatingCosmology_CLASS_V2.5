#!/usr/bin/env python3
"""
Plotting script for Self-Regulating Cosmology (SRC) results.
Generates plots of H(a), ΔH/H, and other background quantities.

Run from: ~/dev/cosmo-class/
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys
import os

def read_class_background(filename):
    """
    Read CLASS background.dat file.
    Returns tuple: (data_dict, header_lines)
    
    CLASS format:
    - Header line 4 contains column names: "#    1:z  2:proper time ... 4:H [1/Mpc] ..."
    - Data starts after header
    """
    header_lines = []
    data_lines = []
    
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # Find the header line with column names (format: "#    1:z  2:proper time ...")
    header_idx = None
    for i, line in enumerate(lines):
        if line.strip().startswith('#') and ':' in line and ('z' in line.lower() or 'a' in line.lower()):
            header_idx = i
            break
    
    if header_idx is None:
        raise ValueError(f"Could not find column header in {filename}")
    
    # Extract column names from header
    header_line = lines[header_idx]
    # Format: "#    1:z  2:proper time [Gyr]  3:conf. time [Mpc]  4:H [1/Mpc] ..."
    # Extract names after colons
    col_names = []
    parts = header_line.strip('#').split()
    current_name = None
    for part in parts:
        if ':' in part:
            # Extract name after colon (e.g., "1:z" -> "z", "4:H" -> "H")
            if current_name:
                col_names.append(current_name.strip())
            name = part.split(':', 1)[1]  # Split only on first colon
            current_name = name
        elif current_name:  # Continuation of previous column name
            current_name += ' ' + part
    if current_name:
        col_names.append(current_name.strip())
    
    # Collect all header lines
    header_lines = lines[:header_idx+1]
    
    # Collect data lines (everything after header)
    data_lines = lines[header_idx+1:]
    
    # Parse data
    if not data_lines:
        raise ValueError(f"No data lines found in {filename}")
    
    data = np.loadtxt(data_lines, dtype=float)
    
    # Create dictionary with column names
    result = {}
    for i, col in enumerate(col_names):
        if i < data.shape[1]:
            # Clean column name (remove brackets, etc.)
            clean_name = col.split('[')[0].strip()
            result[clean_name] = data[:, i]
    
    # Also store by index for fallback
    for i in range(data.shape[1]):
        result[f'col_{i}'] = data[:, i]
    
    return result, header_lines

def read_class_thermodynamics(filename):
    """
    Read CLASS thermodynamics.dat file.
    Returns tuple: (data_dict, header_lines)
    
    Similar format to background.dat
    """
    return read_class_background(filename)  # Same format

def plot_H_comparison(file_lcdm, file_src, output_dir='plots'):
    """
    Plot H(a) and ΔH/H comparison between LCDM and SRC.
    
    Parameters:
    - file_lcdm: path to LCDM background.dat
    - file_src: path to SRC background.dat
    - output_dir: directory to save plots
    """
    # Read data
    try:
        data_lcdm, _ = read_class_background(file_lcdm)
        data_src, _ = read_class_background(file_src)
    except Exception as e:
        print(f"Error reading files: {e}")
        import traceback
        traceback.print_exc()
        return
    
    # Find redshift and H columns
    # CLASS format: column 1 is 'z', column 4 is 'H [1/Mpc]'
    z_key = None
    H_key = None
    
    for key in data_lcdm.keys():
        if key.lower().strip() == 'z':
            z_key = key
        if key.lower().strip() == 'h' or (key.upper().startswith('H') and '[' in key):
            H_key = key
    
    if z_key is None:
        z_key = 'col_0'  # Fallback: first column
    if H_key is None:
        H_key = 'col_3'  # Fallback: fourth column (0-indexed: col_3)
    
    print(f"Using columns: z='{z_key}', H='{H_key}'")
    
    z_lcdm = data_lcdm[z_key]
    H_lcdm = data_lcdm[H_key]
    
    z_src = data_src[z_key]
    H_src = data_src[H_key]
    
    # Interpolate to common z grid for comparison
    # Focus on low-z where SRC is active (z < 2)
    z_max = min(z_lcdm.max(), z_src.max(), 2.0)
    z_min = max(z_lcdm.min(), z_src.min(), 0.0)
    z_common = np.linspace(z_min, z_max, 1000)
    
    H_lcdm_interp = np.interp(z_common, z_lcdm, H_lcdm)
    H_src_interp = np.interp(z_common, z_src, H_src)
    
    # Calculate ΔH/H exactly as specified: (H_SRC - H_LCDM) / H_LCDM
    delta_H_over_H = (H_src_interp - H_lcdm_interp) / H_lcdm_interp * 100  # in percent
    
    # Convert z to a for reference
    a_common = 1.0 / (1.0 + z_common)
    
    # Create plots
    Path(output_dir).mkdir(exist_ok=True, parents=True)
    
    # Plot 1: H(z) comparison
    plt.figure(figsize=(10, 6))
    plt.plot(z_common, H_lcdm_interp, 'b-', label='ΛCDM', linewidth=2)
    plt.plot(z_common, H_src_interp, 'r--', label='SRC', linewidth=2)
    plt.xlabel('Redshift z', fontsize=12)
    plt.ylabel('H(z) [1/Mpc]', fontsize=12)
    plt.title('Hubble Parameter: SRC vs ΛCDM', fontsize=14)
    plt.legend(fontsize=11)
    plt.grid(True, alpha=0.3)
    plt.xlim(0, 2.0)  # Focus on low-z where SRC is active
    plt.tight_layout()
    plt.savefig(f'{output_dir}/H_z_comparison.png', dpi=150)
    plt.close()
    
    # Plot 2: ΔH/H in percent
    plt.figure(figsize=(10, 6))
    plt.plot(z_common, delta_H_over_H, 'g-', linewidth=2)
    plt.axhline(y=0, color='k', linestyle=':', alpha=0.5)
    plt.xlabel('Redshift z', fontsize=12)
    plt.ylabel('ΔH/H [%]', fontsize=12)
    plt.title('Relative Hubble Parameter Difference: (H_SRC - H_ΛCDM) / H_ΛCDM', fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.xlim(0, 2.0)
    # Auto-scale y-axis based on data range
    y_range = np.max(np.abs(delta_H_over_H))
    if y_range > 0:
        plt.ylim(-1.5 * y_range, 1.5 * y_range)
    plt.tight_layout()
    plt.savefig(f'{output_dir}/delta_H_over_H.png', dpi=150)
    plt.close()
    
    # Plot 3: Focus on z ≈ 0.25 region
    mask = (z_common >= 0.0) & (z_common <= 0.5)
    if np.any(mask):
        plt.figure(figsize=(10, 6))
        plt.plot(z_common[mask], delta_H_over_H[mask], 'g-', linewidth=2, marker='o', markersize=3)
        plt.axhline(y=0, color='k', linestyle=':', alpha=0.5)
        plt.axvline(x=0.25, color='r', linestyle='--', alpha=0.5, label='z ≈ 0.25 (SRC peak)')
        plt.xlabel('Redshift z', fontsize=12)
        plt.ylabel('ΔH/H [%]', fontsize=12)
        plt.title('SRC Signal: ΔH/H around z ≈ 0.25', fontsize=14)
        plt.legend(fontsize=11)
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(f'{output_dir}/delta_H_over_H_zoom.png', dpi=150)
        plt.close()
    
    print(f"✓ Plots saved to {output_dir}/")
    print(f"  - H_z_comparison.png")
    print(f"  - delta_H_over_H.png")
    print(f"  - delta_H_over_H_zoom.png")
    
    # Print summary statistics
    if np.any(mask):
        peak_idx = np.argmax(np.abs(delta_H_over_H[mask]))
        peak_z = z_common[mask][peak_idx]
        peak_delta = delta_H_over_H[mask][peak_idx]
        
        print(f"\nSummary:")
        print(f"  Peak |ΔH/H|: {abs(peak_delta):.4f}% at z = {peak_z:.3f}")
        print(f"  Expected: ~0.3-1.0% around z ≈ 0.25")
    else:
        print("\nWarning: No data in z ∈ [0, 0.5] range")

def main():
    """Main function"""
    
    # Get script directory
    script_dir = Path(__file__).parent.absolute()
    os.chdir(script_dir)
    
    # Default file paths
    file_lcdm = Path("output/lcdm/lcdm_background.dat")
    file_src = Path("output/src/src_background.dat")
    
    # Allow command-line arguments
    if len(sys.argv) >= 3:
        file_lcdm = Path(sys.argv[1])
        file_src = Path(sys.argv[2])
    
    if len(sys.argv) >= 4:
        output_dir = sys.argv[3]
    else:
        output_dir = "plots"
    
    if not file_lcdm.exists():
        print(f"Error: LCDM file not found: {file_lcdm}")
        print("Usage: python3 plot_SRC_results.py [lcdm_file] [src_file] [output_dir]")
        sys.exit(1)
    
    if not file_src.exists():
        print(f"Error: SRC file not found: {file_src}")
        print("Usage: python3 plot_SRC_results.py [lcdm_file] [src_file] [output_dir]")
        sys.exit(1)
    
    print("=" * 70)
    print("SRC Results Plotting")
    print("=" * 70)
    print(f"LCDM file: {file_lcdm}")
    print(f"SRC file:  {file_src}")
    print()
    
    plot_H_comparison(file_lcdm, file_src, output_dir)

if __name__ == "__main__":
    main()
