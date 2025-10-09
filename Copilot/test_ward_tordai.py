"""
Simple test script for Ward-Tordai model implementation
"""

import numpy as np
import matplotlib.pyplot as plt
from Ward_Tordai import WardTordaiModel, create_example_surfactant

def simple_test():
    """Simple test of the Ward-Tordai model"""
    print("Testing Ward-Tordai Model Implementation")
    print("=" * 50)
    
    # Create a simple model without Langmuir isotherm
    D = 1e-9  # m²/s
    c0 = 0.1  # mol/m³
    model = WardTordaiModel(D=D, c0=c0)
    
    print(f"Model parameters:")
    print(f"  D = {D:.2e} m²/s")
    print(f"  c0 = {c0} mol/m³")
    print()
    
    # Test short time simulation
    t_max = 100  # seconds
    n_points = 50
    
    print("Running short simulation...")
    t, Gamma, c_surface = model.solve_ward_tordai(t_max, n_points)
    
    # Print first few results
    print("First few time points:")
    for i in range(0, min(10, len(t))):
        print(f"  t = {t[i]:6.1f} s: Γ = {Gamma[i]*1e6:8.3f} µmol/m², c_s = {c_surface[i]:8.4f} mol/m³")
    
    # Test with Langmuir isotherm
    print("\nTesting with Langmuir isotherm...")
    model_langmuir = create_example_surfactant()
    t2, gamma, Gamma2 = model_langmuir.simulate_dynamic_surface_tension(t_max, n_points)
    
    print("Dynamic surface tension results:")
    for i in [5, 10, 20, -1]:
        print(f"  t = {t2[i]:6.1f} s: γ = {gamma[i]:.5f} N/m, Γ = {Gamma2[i]*1e6:8.3f} µmol/m²")
    
    # Simple plot
    plt.figure(figsize=(12, 4))
    
    plt.subplot(1, 3, 1)
    plt.semilogx(t[1:], Gamma[1:] * 1e6, 'b-', linewidth=2)
    plt.xlabel('Time (s)')
    plt.ylabel('Surface Excess (µmol/m²)')
    plt.title('Without Langmuir')
    plt.grid(True, alpha=0.3)
    
    plt.subplot(1, 3, 2)
    plt.semilogx(t2[1:], Gamma2[1:] * 1e6, 'r-', linewidth=2)
    plt.xlabel('Time (s)')
    plt.ylabel('Surface Excess (µmol/m²)')
    plt.title('With Langmuir')
    plt.grid(True, alpha=0.3)
    
    plt.subplot(1, 3, 3)
    plt.semilogx(t2[1:], gamma[1:], 'g-', linewidth=2)
    plt.xlabel('Time (s)')
    plt.ylabel('Surface Tension (N/m)')
    plt.title('Dynamic Surface Tension')
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('ward_tordai_test.png', dpi=150, bbox_inches='tight')
    plt.show()
    
    print("\nTest completed successfully!")
    print("Plot saved as 'ward_tordai_test.png'")

if __name__ == "__main__":
    simple_test()