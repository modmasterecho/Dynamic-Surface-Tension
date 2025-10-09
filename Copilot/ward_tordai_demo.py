"""
Comprehensive demonstration of the Ward-Tordai model implementation
This script shows various features and use cases of the model.
"""

import numpy as np
import matplotlib.pyplot as plt
from Ward_Tordai import WardTordaiModel, create_example_surfactant

def demonstrate_different_scenarios():
    """Demonstrate Ward-Tordai model with different scenarios"""
    
    print("Ward-Tordai Model - Comprehensive Demonstration")
    print("=" * 60)
    
    # Scenario 1: Diffusion-limited adsorption (no isotherm)
    print("\n1. Diffusion-limited adsorption (no equilibrium isotherm)")
    print("-" * 50)
    
    D1 = 5e-10  # m²/s
    c0_1 = 0.05  # mol/m³
    model1 = WardTordaiModel(D=D1, c0=c0_1)
    
    t_max = 1000  # seconds
    t1, Gamma1, c_surf1 = model1.solve_ward_tordai(t_max, 200)
    
    print(f"Parameters: D = {D1:.2e} m²/s, c0 = {c0_1} mol/m³")
    print(f"Final surface excess: {Gamma1[-1]*1e6:.2f} µmol/m²")
    
    # Scenario 2: Langmuir isotherm with strong adsorption
    print("\n2. Langmuir isotherm with strong adsorption")
    print("-" * 50)
    
    D2 = 3e-10  # m²/s
    c0_2 = 0.1   # mol/m³
    Gamma_max2 = 5e-6  # mol/m²
    K2 = 5000   # m³/mol (strong adsorption)
    
    model2 = WardTordaiModel(D=D2, c0=c0_2, Gamma_max=Gamma_max2, K=K2)
    t2, gamma2, Gamma2 = model2.simulate_dynamic_surface_tension(t_max, 200)
    
    print(f"Parameters: D = {D2:.2e} m²/s, c0 = {c0_2} mol/m³")
    print(f"Γ_max = {Gamma_max2*1e6:.1f} µmol/m², K = {K2} m³/mol")
    print(f"Final: γ = {gamma2[-1]:.4f} N/m, Γ = {Gamma2[-1]*1e6:.2f} µmol/m²")
    
    # Scenario 3: Langmuir isotherm with weak adsorption
    print("\n3. Langmuir isotherm with weak adsorption")
    print("-" * 50)
    
    D3 = 8e-10  # m²/s
    c0_3 = 0.2   # mol/m³
    Gamma_max3 = 3e-6  # mol/m²
    K3 = 100    # m³/mol (weak adsorption)
    
    model3 = WardTordaiModel(D=D3, c0=c0_3, Gamma_max=Gamma_max3, K=K3)
    t3, gamma3, Gamma3 = model3.simulate_dynamic_surface_tension(t_max, 200)
    
    print(f"Parameters: D = {D3:.2e} m²/s, c0 = {c0_3} mol/m³")
    print(f"Γ_max = {Gamma_max3*1e6:.1f} µmol/m², K = {K3} m³/mol")
    print(f"Final: γ = {gamma3[-1]:.4f} N/m, Γ = {Gamma3[-1]*1e6:.2f} µmol/m²")
    
    # Create comprehensive plots
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    
    # Surface excess plots
    axes[0, 0].semilogx(t1[1:], Gamma1[1:] * 1e6, 'b-', linewidth=2, label='Diffusion-limited')
    axes[0, 0].set_xlabel('Time (s)')
    axes[0, 0].set_ylabel('Surface Excess (µmol/m²)')
    axes[0, 0].set_title('Scenario 1: Diffusion-limited')
    axes[0, 0].grid(True, alpha=0.3)
    axes[0, 0].legend()
    
    axes[0, 1].semilogx(t2[1:], Gamma2[1:] * 1e6, 'r-', linewidth=2, label='Strong adsorption')
    axes[0, 1].set_xlabel('Time (s)')
    axes[0, 1].set_ylabel('Surface Excess (µmol/m²)')
    axes[0, 1].set_title('Scenario 2: Strong Adsorption')
    axes[0, 1].grid(True, alpha=0.3)
    axes[0, 1].legend()
    
    axes[0, 2].semilogx(t3[1:], Gamma3[1:] * 1e6, 'g-', linewidth=2, label='Weak adsorption')
    axes[0, 2].set_xlabel('Time (s)')
    axes[0, 2].set_ylabel('Surface Excess (µmol/m²)')
    axes[0, 2].set_title('Scenario 3: Weak Adsorption')
    axes[0, 2].grid(True, alpha=0.3)
    axes[0, 2].legend()
    
    # Surface tension plots
    axes[1, 0].semilogx(t1[1:], np.full_like(t1[1:], 0.072), 'b--', linewidth=2, label='Pure water')
    axes[1, 0].set_xlabel('Time (s)')
    axes[1, 0].set_ylabel('Surface Tension (N/m)')
    axes[1, 0].set_title('Surface Tension (constant)')
    axes[1, 0].grid(True, alpha=0.3)
    axes[1, 0].legend()
    
    axes[1, 1].semilogx(t2[1:], gamma2[1:], 'r-', linewidth=2, label='Strong adsorption')
    axes[1, 1].set_xlabel('Time (s)')
    axes[1, 1].set_ylabel('Surface Tension (N/m)')
    axes[1, 1].set_title('Dynamic Surface Tension')
    axes[1, 1].grid(True, alpha=0.3)
    axes[1, 1].legend()
    
    axes[1, 2].semilogx(t3[1:], gamma3[1:], 'g-', linewidth=2, label='Weak adsorption')
    axes[1, 2].set_xlabel('Time (s)')
    axes[1, 2].set_ylabel('Surface Tension (N/m)')
    axes[1, 2].set_title('Dynamic Surface Tension')
    axes[1, 2].grid(True, alpha=0.3)
    axes[1, 2].legend()
    
    plt.tight_layout()
    plt.savefig('ward_tordai_comprehensive.png', dpi=150, bbox_inches='tight')
    plt.show()
    
    # Compare all scenarios on one plot
    plt.figure(figsize=(12, 5))
    
    plt.subplot(1, 2, 1)
    plt.semilogx(t1[1:], Gamma1[1:] * 1e6, 'b-', linewidth=2, label='Diffusion-limited')
    plt.semilogx(t2[1:], Gamma2[1:] * 1e6, 'r-', linewidth=2, label='Strong adsorption (K=5000)')
    plt.semilogx(t3[1:], Gamma3[1:] * 1e6, 'g-', linewidth=2, label='Weak adsorption (K=100)')
    plt.xlabel('Time (s)')
    plt.ylabel('Surface Excess (µmol/m²)')
    plt.title('Surface Excess Comparison')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.subplot(1, 2, 2)
    plt.semilogx(t2[1:], (0.072 - gamma2[1:]) * 1000, 'r-', linewidth=2, label='Strong adsorption')
    plt.semilogx(t3[1:], (0.072 - gamma3[1:]) * 1000, 'g-', linewidth=2, label='Weak adsorption')
    plt.xlabel('Time (s)')
    plt.ylabel('Surface Pressure (mN/m)')
    plt.title('Surface Pressure Comparison')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('ward_tordai_comparison.png', dpi=150, bbox_inches='tight')
    plt.show()
    
    print(f"\\nPlots saved as 'ward_tordai_comprehensive.png' and 'ward_tordai_comparison.png'")

def validate_short_time_behavior():
    """Validate against analytical short-time solution"""
    print("\\n4. Validation against analytical short-time solution")
    print("-" * 50)
    
    # Create model
    D = 1e-9  # m²/s
    c0 = 0.1  # mol/m³
    model = WardTordaiModel(D=D, c0=c0)
    
    # Short time simulation
    t_max = 10  # seconds
    t, Gamma, _ = model.solve_ward_tordai(t_max, 100)
    
    # Analytical solution for short times
    Gamma_analytical = np.sqrt(D / np.pi) * c0 * np.sqrt(t)
    
    # Compare
    plt.figure(figsize=(10, 4))
    
    plt.subplot(1, 2, 1)
    plt.plot(t, Gamma * 1e6, 'b-', linewidth=2, label='Numerical')
    plt.plot(t, Gamma_analytical * 1e6, 'r--', linewidth=2, label='Analytical')
    plt.xlabel('Time (s)')
    plt.ylabel('Surface Excess (µmol/m²)')
    plt.title('Linear Scale Comparison')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.subplot(1, 2, 2)
    plt.loglog(t[1:], Gamma[1:] * 1e6, 'b-', linewidth=2, label='Numerical')
    plt.loglog(t[1:], Gamma_analytical[1:] * 1e6, 'r--', linewidth=2, label='Analytical (√t)')
    plt.xlabel('Time (s)')
    plt.ylabel('Surface Excess (µmol/m²)')
    plt.title('Log-Log Scale (slope should be 0.5)')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('ward_tordai_validation.png', dpi=150, bbox_inches='tight')
    plt.show()
    
    # Calculate relative error
    rel_error = np.abs(Gamma - Gamma_analytical) / (Gamma_analytical + 1e-12)
    print(f"Maximum relative error: {np.max(rel_error[1:]) * 100:.2f}%")
    print(f"Average relative error: {np.mean(rel_error[1:]) * 100:.2f}%")

if __name__ == "__main__":
    demonstrate_different_scenarios()
    validate_short_time_behavior()
    print("\\nDemonstration completed successfully!")