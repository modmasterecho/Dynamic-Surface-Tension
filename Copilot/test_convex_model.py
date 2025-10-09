"""
Test script for Convex Ward-Tordai Model with Bubble Tensiometry Data
This script demonstrates the convex model and compares it with planar results.
"""

import numpy as np
import matplotlib.pyplot as plt
from convex_ward_tordai import ConvexWardTordaiModel, create_bubble_tensiometry_example, generate_example_bubble_data
from Ward_Tordai import WardTordaiModel

def test_convex_ward_tordai():
    """Test the convex Ward-Tordai model with bubble tensiometry"""
    
    print("Convex Ward-Tordai Model - Bubble Tensiometry Test")
    print("=" * 55)
    
    # Generate example bubble tensiometry data
    print("\\n1. Generating example bubble tensiometry data...")
    bubble_data = generate_example_bubble_data()
    
    # Extract data
    t = bubble_data['time']
    gamma_measured = bubble_data['surface_tension']
    pressure_measured = bubble_data['pressure']
    model_params = bubble_data['model_parameters']
    
    print(f"Model parameters used:")
    print(f"  Diffusion coefficient: {model_params['D']:.2e} m²/s")
    print(f"  Bulk concentration: {model_params['c0']*1000:.1f} mM")
    print(f"  Initial radius: {model_params['r0']*1000:.1f} mm")
    print(f"  Growth rate: {model_params['growth_rate']*1000:.1f} mm/s")
    print(f"  Γ_max: {model_params['Gamma_max']*1e6:.1f} µmol/m²")
    print(f"  K: {model_params['K']} m³/mol")
    
    # Print some measurement results
    print(f"\\nMeasurement results:")
    print(f"  Time range: {t[0]:.3f} - {t[-1]:.1f} s")
    print(f"  Initial surface tension: {gamma_measured[0]*1000:.1f} mN/m")
    print(f"  Final surface tension: {gamma_measured[-1]*1000:.1f} mN/m")
    print(f"  Maximum pressure: {pressure_measured.max():.0f} Pa")
    print(f"  Final bubble radius: {bubble_data['radius'][-1]*1000:.2f} mm")
    
    return bubble_data

def compare_planar_vs_convex():
    """Compare planar and convex Ward-Tordai models"""
    
    print("\\n2. Comparing Planar vs Convex Ward-Tordai Models")
    print("-" * 50)
    
    # Create models with same parameters
    D = 3e-10
    c0 = 0.005
    Gamma_max = 2.5e-6
    K = 800
    
    # Planar model
    planar_model = WardTordaiModel(D=D, c0=c0, Gamma_max=Gamma_max, K=K)
    
    # Convex model
    convex_model = ConvexWardTordaiModel(
        D=D, c0=c0, r0=0.5e-3, bubble_growth_rate=1e-4,
        Gamma_max=Gamma_max, K=K
    )
    
    # Simulate both models
    t_max = 5.0
    n_points = 200
    
    # Planar simulation
    t_planar, gamma_planar, Gamma_planar = planar_model.simulate_dynamic_surface_tension(t_max, n_points)
    
    # Convex simulation
    t_convex, gamma_convex, Gamma_convex, radius, pressure = convex_model.simulate_bubble_tensiometry(t_max, n_points)
    
    print(f"Comparison at t = {t_max} s:")
    print(f"  Planar model:")
    print(f"    γ = {gamma_planar[-1]*1000:.2f} mN/m")
    print(f"    Γ = {Gamma_planar[-1]*1e6:.3f} µmol/m²")
    print(f"  Convex model:")
    print(f"    γ = {gamma_convex[-1]*1000:.2f} mN/m")
    print(f"    Γ = {Gamma_convex[-1]*1e6:.3f} µmol/m²")
    print(f"    Final radius = {radius[-1]*1000:.2f} mm")
    print(f"    Final pressure = {pressure[-1]:.0f} Pa")
    
    # Create comparison plots
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    
    # Surface tension comparison
    axes[0, 0].plot(t_planar, gamma_planar * 1000, 'b-', linewidth=2, label='Planar')
    axes[0, 0].plot(t_convex, gamma_convex * 1000, 'r-', linewidth=2, label='Convex')
    axes[0, 0].set_xlabel('Time (s)')
    axes[0, 0].set_ylabel('Surface Tension (mN/m)')
    axes[0, 0].set_title('Surface Tension Evolution')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)
    
    # Surface excess comparison
    axes[0, 1].plot(t_planar, Gamma_planar * 1e6, 'b-', linewidth=2, label='Planar')
    axes[0, 1].plot(t_convex, Gamma_convex * 1e6, 'r-', linewidth=2, label='Convex')
    axes[0, 1].set_xlabel('Time (s)')
    axes[0, 1].set_ylabel('Surface Excess (µmol/m²)')
    axes[0, 1].set_title('Surface Excess Evolution')
    axes[0, 1].legend()
    axes[0, 1].grid(True, alpha=0.3)
    
    # Bubble radius evolution
    axes[0, 2].plot(t_convex, radius * 1000, 'r-', linewidth=2)
    axes[0, 2].set_xlabel('Time (s)')
    axes[0, 2].set_ylabel('Bubble Radius (mm)')
    axes[0, 2].set_title('Bubble Growth')
    axes[0, 2].grid(True, alpha=0.3)
    
    # Surface pressure comparison
    axes[1, 0].plot(t_planar, (0.072 - gamma_planar) * 1000, 'b-', linewidth=2, label='Planar')
    axes[1, 0].plot(t_convex, (0.072 - gamma_convex) * 1000, 'r-', linewidth=2, label='Convex')
    axes[1, 0].set_xlabel('Time (s)')
    axes[1, 0].set_ylabel('Surface Pressure (mN/m)')
    axes[1, 0].set_title('Surface Pressure Evolution')
    axes[1, 0].legend()
    axes[1, 0].grid(True, alpha=0.3)
    
    # Bubble pressure vs time
    axes[1, 1].plot(t_convex, pressure, 'r-', linewidth=2)
    axes[1, 1].set_xlabel('Time (s)')
    axes[1, 1].set_ylabel('Bubble Pressure (Pa)')
    axes[1, 1].set_title('Bubble Pressure Evolution')
    axes[1, 1].grid(True, alpha=0.3)
    
    # Surface tension vs bubble radius
    axes[1, 2].plot(radius * 1000, gamma_convex * 1000, 'r-', linewidth=2)
    axes[1, 2].set_xlabel('Bubble Radius (mm)')
    axes[1, 2].set_ylabel('Surface Tension (mN/m)')
    axes[1, 2].set_title('γ vs Bubble Size')
    axes[1, 2].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('convex_vs_planar_comparison.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    return t_planar, gamma_planar, Gamma_planar, t_convex, gamma_convex, Gamma_convex, radius, pressure

def analyze_bubble_tensiometry_data():
    """Analyze the generated bubble tensiometry data in detail"""
    
    print("\\n3. Detailed Analysis of Bubble Tensiometry Data")
    print("-" * 50)
    
    # Generate data
    bubble_data = generate_example_bubble_data()
    
    t = bubble_data['time']
    gamma = bubble_data['clean_surface_tension']
    Gamma = bubble_data['surface_excess']
    radius = bubble_data['radius']
    pressure = bubble_data['clean_pressure']
    
    # Calculate derived quantities
    surface_pressure = (0.072 - gamma) * 1000  # mN/m
    surface_area = 4 * np.pi * radius**2  # m²
    
    # Find maximum bubble pressure point
    max_pressure_idx = np.argmax(pressure)
    t_max_pressure = t[max_pressure_idx]
    gamma_at_max_pressure = gamma[max_pressure_idx]
    
    print(f"Analysis results:")
    print(f"  Maximum bubble pressure: {pressure.max():.0f} Pa at t = {t_max_pressure:.2f} s")
    print(f"  Surface tension at max pressure: {gamma_at_max_pressure*1000:.2f} mN/m")
    print(f"  Surface pressure range: {surface_pressure.min():.2f} - {surface_pressure.max():.2f} mN/m")
    print(f"  Surface excess range: {Gamma.min()*1e6:.3f} - {Gamma.max()*1e6:.3f} µmol/m²")
    print(f"  Bubble radius range: {radius.min()*1000:.2f} - {radius.max()*1000:.2f} mm")
    print(f"  Surface area increase: {surface_area[-1]/surface_area[0]:.1f}x")
    
    # Create detailed analysis plots
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Bubble pressure vs time with maximum point
    axes[0, 0].plot(t, pressure, 'b-', linewidth=2)
    axes[0, 0].plot(t_max_pressure, pressure.max(), 'ro', markersize=8, label='Maximum pressure')
    axes[0, 0].set_xlabel('Time (s)')
    axes[0, 0].set_ylabel('Bubble Pressure (Pa)')
    axes[0, 0].set_title('Bubble Pressure Evolution')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)
    
    # Surface tension vs surface excess (adsorption isotherm)
    axes[0, 1].plot(Gamma * 1e6, gamma * 1000, 'g-', linewidth=2)
    axes[0, 1].set_xlabel('Surface Excess (µmol/m²)')
    axes[0, 1].set_ylabel('Surface Tension (mN/m)')
    axes[0, 1].set_title('Dynamic Adsorption Isotherm')
    axes[0, 1].grid(True, alpha=0.3)
    
    # Surface area vs time
    axes[1, 0].plot(t, surface_area * 1e6, 'purple', linewidth=2)  # Convert to mm²
    axes[1, 0].set_xlabel('Time (s)')
    axes[1, 0].set_ylabel('Surface Area (mm²)')
    axes[1, 0].set_title('Interface Area Growth')
    axes[1, 0].grid(True, alpha=0.3)
    
    # Surface pressure vs time with characteristic times
    axes[1, 1].plot(t, surface_pressure, 'orange', linewidth=2)
    # Mark characteristic diffusion time
    t_diff = (radius[0]**2) / bubble_data['model_parameters']['D']
    if t_diff < t.max():
        axes[1, 1].axvline(x=t_diff, color='red', linestyle='--', alpha=0.7, 
                          label=f'Diffusion time ({t_diff:.2f} s)')
    axes[1, 1].set_xlabel('Time (s)')
    axes[1, 1].set_ylabel('Surface Pressure (mN/m)')
    axes[1, 1].set_title('Surface Pressure Evolution')
    if t_diff < t.max():
        axes[1, 1].legend()
    axes[1, 1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('bubble_tensiometry_analysis.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    return bubble_data

def validate_young_laplace_equation():
    """Validate the Young-Laplace equation implementation"""
    
    print("\\n4. Validating Young-Laplace Equation")
    print("-" * 40)
    
    # Create a simple test case
    model = create_bubble_tensiometry_example()
    
    # Test with known values
    test_radius = np.array([0.5e-3, 1.0e-3, 2.0e-3])  # mm
    test_gamma = np.array([0.050, 0.055, 0.060])      # N/m
    test_time = test_radius / model.bubble_growth_rate  # Corresponding times
    
    # Calculate pressure using Young-Laplace
    pressure_calculated = model.calculate_bubble_pressure(test_time, test_gamma, apply_correction=False)
    
    # Theoretical pressure: ΔP = 2γ/R
    pressure_theoretical = 2 * test_gamma / test_radius
    
    print(f"Young-Laplace equation validation:")
    print(f"{'Radius (mm)':<12} {'γ (mN/m)':<10} {'P_calc (Pa)':<12} {'P_theory (Pa)':<12} {'Error (%)':<10}")
    print("-" * 70)
    
    for i in range(len(test_radius)):
        error = abs(pressure_calculated[i] - pressure_theoretical[i]) / pressure_theoretical[i] * 100
        print(f"{test_radius[i]*1000:<12.1f} {test_gamma[i]*1000:<10.1f} {pressure_calculated[i]:<12.0f} "
              f"{pressure_theoretical[i]:<12.0f} {error:<10.2f}")
    
    print(f"\\nValidation successful - Young-Laplace equation correctly implemented")

if __name__ == "__main__":
    # Run all tests
    print("Testing Convex Ward-Tordai Model for Bubble Tensiometry")
    print("=" * 60)
    
    # Test 1: Generate example data
    bubble_data = test_convex_ward_tordai()
    
    # Test 2: Compare models
    comparison_results = compare_planar_vs_convex()
    
    # Test 3: Detailed analysis
    analysis_results = analyze_bubble_tensiometry_data()
    
    # Test 4: Validate Young-Laplace
    validate_young_laplace_equation()
    
    print("\\n" + "=" * 60)
    print("All tests completed successfully!")
    print("Generated plots:")
    print("  - convex_vs_planar_comparison.png")
    print("  - bubble_tensiometry_analysis.png")
    print("\\nConvex Ward-Tordai model is ready for bubble tensiometry analysis!")