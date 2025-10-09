"""
Realistic Ward-Tordai model demonstration with typical surfactant parameters
"""

import numpy as np
import matplotlib.pyplot as plt
from Ward_Tordai import WardTordaiModel

def create_realistic_examples():
    """Create Ward-Tordai models with realistic surfactant parameters"""
    
    print("Ward-Tordai Model - Realistic Surfactant Examples")
    print("=" * 55)
    
    # Example 1: C12E4 (nonionic surfactant)
    print("\\n1. C12E4 (polyethylene glycol dodecyl ether)")
    print("-" * 45)
    
    # Realistic parameters for C12E4
    D_C12E4 = 4e-10      # m²/s (typical for surfactants)
    c0_C12E4 = 0.01      # mol/m³ (0.01 M)
    Gamma_max_C12E4 = 2.5e-6  # mol/m² (typical maximum packing)
    K_C12E4 = 500        # m³/mol (moderate adsorption)
    
    model_C12E4 = WardTordaiModel(
        D=D_C12E4, 
        c0=c0_C12E4, 
        Gamma_max=Gamma_max_C12E4, 
        K=K_C12E4
    )
    
    # Simulate for 1 hour
    t_max = 3600  # 1 hour
    t1, gamma1, Gamma1 = model_C12E4.simulate_dynamic_surface_tension(t_max, 300)
    
    print(f"Parameters:")
    print(f"  D = {D_C12E4:.1e} m²/s")
    print(f"  c0 = {c0_C12E4} mol/m³ = {c0_C12E4*1000:.1f} mM")
    print(f"  Γ_max = {Gamma_max_C12E4*1e6:.1f} µmol/m²")
    print(f"  K = {K_C12E4} m³/mol")
    print(f"Results at t = 1 hour:")
    print(f"  γ = {gamma1[-1]:.4f} N/m")
    print(f"  Γ = {Gamma1[-1]*1e6:.2f} µmol/m²")
    print(f"  Surface pressure = {(0.072 - gamma1[-1])*1000:.1f} mN/m")
    
    # Example 2: SDS (anionic surfactant)
    print("\\n2. SDS (sodium dodecyl sulfate)")
    print("-" * 35)
    
    # Realistic parameters for SDS
    D_SDS = 6e-10        # m²/s
    c0_SDS = 0.005       # mol/m³ (5 mM, below CMC)
    Gamma_max_SDS = 3.2e-6  # mol/m²
    K_SDS = 1000         # m³/mol (strong adsorption)
    
    model_SDS = WardTordaiModel(
        D=D_SDS, 
        c0=c0_SDS, 
        Gamma_max=Gamma_max_SDS, 
        K=K_SDS
    )
    
    t2, gamma2, Gamma2 = model_SDS.simulate_dynamic_surface_tension(t_max, 300)
    
    print(f"Parameters:")
    print(f"  D = {D_SDS:.1e} m²/s")
    print(f"  c0 = {c0_SDS} mol/m³ = {c0_SDS*1000:.1f} mM")
    print(f"  Γ_max = {Gamma_max_SDS*1e6:.1f} µmol/m²")
    print(f"  K = {K_SDS} m³/mol")
    print(f"Results at t = 1 hour:")
    print(f"  γ = {gamma2[-1]:.4f} N/m")
    print(f"  Γ = {Gamma2[-1]*1e6:.2f} µmol/m²")
    print(f"  Surface pressure = {(0.072 - gamma2[-1])*1000:.1f} mN/m")
    
    # Example 3: Low concentration case
    print("\\n3. Low concentration example")
    print("-" * 30)
    
    D_low = 5e-10        # m²/s
    c0_low = 0.001       # mol/m³ (1 mM)
    Gamma_max_low = 2.0e-6  # mol/m²
    K_low = 200          # m³/mol (weak adsorption)
    
    model_low = WardTordaiModel(
        D=D_low, 
        c0=c0_low, 
        Gamma_max=Gamma_max_low, 
        K=K_low
    )
    
    t3, gamma3, Gamma3 = model_low.simulate_dynamic_surface_tension(t_max, 300)
    
    print(f"Parameters:")
    print(f"  D = {D_low:.1e} m²/s")
    print(f"  c0 = {c0_low} mol/m³ = {c0_low*1000:.1f} mM")
    print(f"  Γ_max = {Gamma_max_low*1e6:.1f} µmol/m²")
    print(f"  K = {K_low} m³/mol")
    print(f"Results at t = 1 hour:")
    print(f"  γ = {gamma3[-1]:.4f} N/m")
    print(f"  Γ = {Gamma3[-1]*1e6:.2f} µmol/m²")
    print(f"  Surface pressure = {(0.072 - gamma3[-1])*1000:.1f} mN/m")
    
    # Create publication-quality plots
    plt.style.use('default')
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Surface excess vs time
    axes[0, 0].semilogx(t1, Gamma1 * 1e6, 'b-', linewidth=2, label='C12E4 (10 mM)')
    axes[0, 0].semilogx(t2, Gamma2 * 1e6, 'r-', linewidth=2, label='SDS (5 mM)')
    axes[0, 0].semilogx(t3, Gamma3 * 1e6, 'g-', linewidth=2, label='Low conc. (1 mM)')
    axes[0, 0].set_xlabel('Time (s)')
    axes[0, 0].set_ylabel('Surface Excess (µmol/m²)')
    axes[0, 0].set_title('Surface Excess Evolution')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)
    
    # Surface tension vs time
    axes[0, 1].semilogx(t1, gamma1 * 1000, 'b-', linewidth=2, label='C12E4')
    axes[0, 1].semilogx(t2, gamma2 * 1000, 'r-', linewidth=2, label='SDS')
    axes[0, 1].semilogx(t3, gamma3 * 1000, 'g-', linewidth=2, label='Low conc.')
    axes[0, 1].axhline(y=72, color='k', linestyle='--', alpha=0.5, label='Pure water')
    axes[0, 1].set_xlabel('Time (s)')
    axes[0, 1].set_ylabel('Surface Tension (mN/m)')
    axes[0, 1].set_title('Dynamic Surface Tension')
    axes[0, 1].legend()
    axes[0, 1].grid(True, alpha=0.3)
    
    # Surface pressure vs time
    axes[1, 0].semilogx(t1, (0.072 - gamma1) * 1000, 'b-', linewidth=2, label='C12E4')
    axes[1, 0].semilogx(t2, (0.072 - gamma2) * 1000, 'r-', linewidth=2, label='SDS')
    axes[1, 0].semilogx(t3, (0.072 - gamma3) * 1000, 'g-', linewidth=2, label='Low conc.')
    axes[1, 0].set_xlabel('Time (s)')
    axes[1, 0].set_ylabel('Surface Pressure (mN/m)')
    axes[1, 0].set_title('Surface Pressure Evolution')
    axes[1, 0].legend()
    axes[1, 0].grid(True, alpha=0.3)
    
    # Short time behavior (first 100 seconds)
    mask = t1 <= 100
    axes[1, 1].loglog(t1[mask], Gamma1[mask] * 1e6, 'b-', linewidth=2, label='C12E4')
    axes[1, 1].loglog(t2[mask], Gamma2[mask] * 1e6, 'r-', linewidth=2, label='SDS')
    axes[1, 1].loglog(t3[mask], Gamma3[mask] * 1e6, 'g-', linewidth=2, label='Low conc.')
    
    # Add √t reference line
    t_ref = np.logspace(-1, 2, 50)
    Gamma_ref = 0.1 * np.sqrt(t_ref)  # Arbitrary scaling for reference
    axes[1, 1].loglog(t_ref, Gamma_ref, 'k--', alpha=0.5, linewidth=1, label='√t behavior')
    
    axes[1, 1].set_xlabel('Time (s)')
    axes[1, 1].set_ylabel('Surface Excess (µmol/m²)')
    axes[1, 1].set_title('Short-time Behavior (log-log)')
    axes[1, 1].legend()
    axes[1, 1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('realistic_ward_tordai_examples.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    print(f"\\nPlot saved as 'realistic_ward_tordai_examples.png'")
    
    # Summary table
    print(f"\\nSummary Table (at t = 1 hour):")
    print(f"{'Surfactant':<12} {'Conc (mM)':<10} {'γ (mN/m)':<10} {'Γ (µmol/m²)':<12} {'Π (mN/m)':<10}")
    print("-" * 60)
    print(f"{'C12E4':<12} {c0_C12E4*1000:<10.1f} {gamma1[-1]*1000:<10.1f} {Gamma1[-1]*1e6:<12.2f} {(0.072-gamma1[-1])*1000:<10.1f}")
    print(f"{'SDS':<12} {c0_SDS*1000:<10.1f} {gamma2[-1]*1000:<10.1f} {Gamma2[-1]*1e6:<12.2f} {(0.072-gamma2[-1])*1000:<10.1f}")
    print(f"{'Low conc.':<12} {c0_low*1000:<10.1f} {gamma3[-1]*1000:<10.1f} {Gamma3[-1]*1e6:<12.2f} {(0.072-gamma3[-1])*1000:<10.1f}")

if __name__ == "__main__":
    create_realistic_examples()
    print("\\nRealistic demonstration completed successfully!")