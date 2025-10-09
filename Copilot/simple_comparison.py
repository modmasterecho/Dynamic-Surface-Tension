"""
Simple example comparing planar and convex Ward-Tordai models
This script demonstrates the key differences between the two approaches.
"""

import numpy as np
import matplotlib.pyplot as plt
from Ward_Tordai import WardTordaiModel
from convex_ward_tordai import ConvexWardTordaiModel

def simple_comparison():
    """Simple comparison of planar vs convex models"""
    
    print("Ward-Tordai Models: Planar vs Convex Comparison")
    print("=" * 50)
    
    # Common parameters
    D = 4e-10      # m²/s
    c0 = 0.005     # mol/m³ (5 mM)
    Gamma_max = 2.5e-6  # mol/m²
    K = 600        # m³/mol
    
    print(f"Common parameters:")
    print(f"  D = {D:.1e} m²/s")
    print(f"  c0 = {c0*1000:.1f} mM")
    print(f"  Γ_max = {Gamma_max*1e6:.1f} µmol/m²")
    print(f"  K = {K} m³/mol")
    
    # Create models
    planar_model = WardTordaiModel(D=D, c0=c0, Gamma_max=Gamma_max, K=K)
    
    convex_model = ConvexWardTordaiModel(
        D=D, c0=c0, r0=0.6e-3, bubble_growth_rate=0.05e-3,
        Gamma_max=Gamma_max, K=K
    )
    
    # Simulation parameters
    t_max = 10.0  # 10 seconds
    n_points = 200
    
    print(f"\\nSimulation setup:")
    print(f"  Time: 0 - {t_max} s")
    print(f"  Points: {n_points}")
    
    # Run simulations
    print("\\nRunning simulations...")
    
    # Planar simulation
    t_planar, gamma_planar, Gamma_planar = planar_model.simulate_dynamic_surface_tension(t_max, n_points)
    
    # Convex simulation  
    t_convex, gamma_convex, Gamma_convex, radius, pressure = convex_model.simulate_bubble_tensiometry(t_max, n_points)
    
    # Results comparison
    print(f"\\nResults at t = {t_max} s:")
    print(f"  Planar model:")
    print(f"    Surface tension: {gamma_planar[-1]*1000:.2f} mN/m")
    print(f"    Surface excess:  {Gamma_planar[-1]*1e6:.3f} µmol/m²")
    print(f"  Convex model:")
    print(f"    Surface tension: {gamma_convex[-1]*1000:.2f} mN/m")
    print(f"    Surface excess:  {Gamma_convex[-1]*1e6:.3f} µmol/m²")
    print(f"    Final radius:    {radius[-1]*1000:.2f} mm")
    print(f"    Final pressure:  {pressure[-1]:.0f} Pa")
    
    # Create comparison plot
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    
    # Surface tension comparison
    axes[0, 0].plot(t_planar, gamma_planar*1000, 'b-', linewidth=2, label='Planar')
    axes[0, 0].plot(t_convex, gamma_convex*1000, 'r-', linewidth=2, label='Convex')
    axes[0, 0].set_xlabel('Time (s)')
    axes[0, 0].set_ylabel('Surface Tension (mN/m)')
    axes[0, 0].set_title('Surface Tension Evolution')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)
    
    # Surface excess comparison
    axes[0, 1].plot(t_planar, Gamma_planar*1e6, 'b-', linewidth=2, label='Planar')
    axes[0, 1].plot(t_convex, Gamma_convex*1e6, 'r-', linewidth=2, label='Convex')
    axes[0, 1].set_xlabel('Time (s)')
    axes[0, 1].set_ylabel('Surface Excess (µmol/m²)')
    axes[0, 1].set_title('Surface Excess Evolution')
    axes[0, 1].legend()
    axes[0, 1].grid(True, alpha=0.3)
    
    # Bubble radius (convex only)
    axes[1, 0].plot(t_convex, radius*1000, 'r-', linewidth=2)
    axes[1, 0].set_xlabel('Time (s)')
    axes[1, 0].set_ylabel('Bubble Radius (mm)')
    axes[1, 0].set_title('Bubble Growth (Convex Model)')
    axes[1, 0].grid(True, alpha=0.3)
    
    # Bubble pressure (convex only)
    axes[1, 1].plot(t_convex, pressure, 'r-', linewidth=2)
    axes[1, 1].set_xlabel('Time (s)')
    axes[1, 1].set_ylabel('Bubble Pressure (Pa)')
    axes[1, 1].set_title('Bubble Pressure (Convex Model)')
    axes[1, 1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('simple_planar_vs_convex_comparison.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    print(f"\\nPlot saved as 'simple_planar_vs_convex_comparison.png'")
    
    # Key differences summary
    print(f"\\nKey Differences:")
    print(f"  1. Interface geometry: Planar vs Spherical")
    print(f"  2. Area change: Constant vs Growing (4πR²)")
    print(f"  3. Pressure: Not measured vs Young-Laplace (2γ/R)")
    print(f"  4. Application: Drop/plate methods vs Bubble tensiometry")
    
    return {
        'planar': {'time': t_planar, 'gamma': gamma_planar, 'Gamma': Gamma_planar},
        'convex': {'time': t_convex, 'gamma': gamma_convex, 'Gamma': Gamma_convex, 
                  'radius': radius, 'pressure': pressure}
    }

if __name__ == "__main__":
    results = simple_comparison()
    print("\\nSimple comparison completed successfully!")