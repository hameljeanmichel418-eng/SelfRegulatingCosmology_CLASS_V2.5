#!/usr/bin/env python3
"""
Plot ΔH/H(z) for SRC gamma0 scan in the physical redshift range z in [0, 1].

Plots runs A2, A3, A4 (gamma0 = 0.020, 0.026, 0.030) to visualize
the late-time modulation of the Hubble parameter.
"""

import numpy as np
import matplotlib.pyplot as plt
import os


def main():
    """Main plotting function."""
    
    # Load the zoomed deltaH data
    print("Loading deltaH_gamma_scan_z0_1.dat...")
    data = np.loadtxt("output/src/deltaH_gamma_scan_z0_1.dat", comments="#")
    
    # Extract columns
    z = data[:, 0]
    dA1 = data[:, 1]
    dA2 = data[:, 2]
    dA3 = data[:, 3]
    dA4 = data[:, 4]
    dA5 = data[:, 5]
    
    print(f"Loaded {len(z)} data points in range z = [{z.min():.3f}, {z.max():.3f}]")
    
    # Create figure
    plt.figure(figsize=(8, 6))
    
    # Plot A2, A3, A4 (gamma0 = 0.020, 0.026, 0.030)
    plt.plot(z, dA2, 'b-', linewidth=1.5, label='γ₀ = 0.020')
    plt.plot(z, dA3, 'r--', linewidth=1.5, label='γ₀ = 0.026')
    plt.plot(z, dA4, 'g-.', linewidth=1.5, label='γ₀ = 0.030')
    
    # Label axes
    plt.xlabel('redshift z', fontsize=12)
    plt.ylabel('ΔH/H [%]', fontsize=12)
    
    # Add title
    plt.title('SRC: late-time modulation of H(z)', fontsize=14, fontweight='bold')
    
    # Add legend
    plt.legend(loc='best', fontsize=11)
    
    # Add grid
    plt.grid(True, alpha=0.3, linestyle='--')
    
    # Set x-axis limits to show full range
    plt.xlim(0.0, 1.0)
    
    # Adjust layout to prevent label cutoff
    plt.tight_layout()
    
    # Create figures directory if it doesn't exist
    output_dir = "figures"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created directory: {output_dir}")
    
    # Save figure
    output_file = "figures/deltaH_over_H_gamma_scan_z0_1.png"
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"Figure saved to: {output_file}")
    
    plt.close()


if __name__ == "__main__":
    main()

