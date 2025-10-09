"""
Bubble Tensiometry Data Analysis with Convex Ward-Tordai Model
This script demonstrates how to analyze real bubble tensiometry measurements.
"""

import numpy as np
import matplotlib.pyplot as plt
from convex_ward_tordai import ConvexWardTordaiModel, create_bubble_tensiometry_example
from Ward_Tordai import WardTordaiModel

def create_realistic_bubble_measurement():
    """
    Create realistic bubble tensiometry data based on typical experimental conditions.
    """
    print("Creating Realistic Bubble Tensiometry Measurement")
    print("=" * 52)
    
    # Experimental parameters for SDS solution
    print("\\nExperimental Setup:")
    print("  Surfactant: Sodium Dodecyl Sulfate (SDS)")
    print("  Concentration: 3 mM (below CMC)")
    print("  Temperature: 25°C")
    print("  Capillary diameter: 1.0 mm")
    print("  Bubble formation time: ~8 seconds")
    
    # Model parameters based on SDS properties
    D_sds = 5.8e-10      # m²/s (literature value for SDS)
    c0_sds = 0.003       # mol/m³ (3 mM)
    r0_capillary = 0.5e-3  # m (0.5 mm radius)
    growth_rate = 0.08e-3  # m/s (slow bubble formation)
    Gamma_max_sds = 3.4e-6  # mol/m² (literature value)
    K_sds = 1200         # m³/mol (estimated from CMC)
    
    model = ConvexWardTordaiModel(
        D=D_sds, c0=c0_sds, r0=r0_capillary,
        bubble_growth_rate=growth_rate,
        Gamma_max=Gamma_max_sds, K=K_sds
    )
    
    # Simulate bubble formation cycle
    t_max = 8.0  # 8 seconds bubble formation
    n_points = 400
    
    t, gamma, Gamma, radius, pressure = model.simulate_bubble_tensiometry(t_max, n_points)
    
    # Add realistic experimental noise
    np.random.seed(123)
    
    # Pressure measurement noise (±5 Pa typical for good tensiometer)
    pressure_noise = np.random.normal(0, 5, len(pressure))
    pressure_measured = pressure + pressure_noise
    
    # Surface tension calculation uncertainty (±0.5 mN/m)
    gamma_uncertainty = np.random.normal(0, 0.0005, len(gamma))
    gamma_measured = gamma + gamma_uncertainty
    
    return {
        'time': t,
        'radius': radius,
        'pressure_clean': pressure,
        'pressure_measured': pressure_measured,
        'surface_tension_clean': gamma,
        'surface_tension_measured': gamma_measured,
        'surface_excess': Gamma,
        'model': model,
        'experimental_params': {
            'surfactant': 'SDS',
            'concentration_mM': c0_sds * 1000,
            'D': D_sds,
            'Gamma_max': Gamma_max_sds,
            'K': K_sds,
            'capillary_radius_mm': r0_capillary * 1000,
            'growth_rate_mm_s': growth_rate * 1000
        }
    }

def analyze_maximum_bubble_pressure():
    """
    Analyze the maximum bubble pressure method for surface tension determination.
    """
    print("\\n" + "="*60)
    print("Maximum Bubble Pressure Method Analysis")
    print("="*60)
    
    # Generate measurement data
    data = create_realistic_bubble_measurement()
    
    t = data['time']
    pressure = data['pressure_measured']
    gamma_true = data['surface_tension_clean']
    radius = data['radius']
    
    # Find maximum pressure point
    max_pressure_idx = np.argmax(pressure)
    t_max = t[max_pressure_idx]
    p_max = pressure[max_pressure_idx]
    r_max = radius[max_pressure_idx]
    gamma_at_max = gamma_true[max_pressure_idx]
    
    # Calculate surface tension from maximum pressure using Young-Laplace
    gamma_from_pressure = p_max * r_max / 2  # From ΔP = 2γ/R
    
    print(f"\\nMaximum Bubble Pressure Analysis:")
    print(f"  Time at maximum pressure: {t_max:.2f} s")
    print(f"  Maximum pressure: {p_max:.1f} Pa")
    print(f"  Bubble radius at max pressure: {r_max*1000:.2f} mm")
    print(f"  True surface tension: {gamma_at_max*1000:.2f} mN/m")
    print(f"  Surface tension from pressure: {gamma_from_pressure*1000:.2f} mN/m")
    print(f"  Error: {abs(gamma_from_pressure - gamma_at_max)*1000:.2f} mN/m")
    print(f"  Relative error: {abs(gamma_from_pressure - gamma_at_max)/gamma_at_max*100:.1f}%")
    
    # Create analysis plots
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Pressure vs time with maximum point marked
    axes[0, 0].plot(t, pressure, 'b-', linewidth=1.5, alpha=0.7, label='Measured')
    axes[0, 0].plot(t, data['pressure_clean'], 'r-', linewidth=2, label='True')
    axes[0, 0].plot(t_max, p_max, 'ro', markersize=10, label=f'Max P = {p_max:.0f} Pa')
    axes[0, 0].set_xlabel('Time (s)')
    axes[0, 0].set_ylabel('Bubble Pressure (Pa)')
    axes[0, 0].set_title('Bubble Pressure Evolution')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)
    
    # Surface tension vs time
    axes[0, 1].plot(t, gamma_true*1000, 'g-', linewidth=2, label='True γ')
    axes[0, 1].axhline(y=gamma_from_pressure*1000, color='red', linestyle='--', 
                       linewidth=2, label=f'From max P: {gamma_from_pressure*1000:.1f} mN/m')
    axes[0, 1].plot(t_max, gamma_at_max*1000, 'ro', markersize=8)
    axes[0, 1].set_xlabel('Time (s)')
    axes[0, 1].set_ylabel('Surface Tension (mN/m)')
    axes[0, 1].set_title('Dynamic Surface Tension')
    axes[0, 1].legend()
    axes[0, 1].grid(True, alpha=0.3)
    
    # Bubble radius vs time
    axes[1, 0].plot(t, radius*1000, 'purple', linewidth=2)
    axes[1, 0].plot(t_max, r_max*1000, 'ro', markersize=8, label=f'R at max P: {r_max*1000:.2f} mm')
    axes[1, 0].set_xlabel('Time (s)')
    axes[1, 0].set_ylabel('Bubble Radius (mm)')
    axes[1, 0].set_title('Bubble Growth')
    axes[1, 0].legend()
    axes[1, 0].grid(True, alpha=0.3)
    
    # Surface excess evolution
    axes[1, 1].plot(t, data['surface_excess']*1e6, 'orange', linewidth=2)
    axes[1, 1].plot(t_max, data['surface_excess'][max_pressure_idx]*1e6, 'ro', markersize=8,
                   label=f'Γ at max P: {data["surface_excess"][max_pressure_idx]*1e6:.2f} µmol/m²')
    axes[1, 1].set_xlabel('Time (s)')
    axes[1, 1].set_ylabel('Surface Excess (µmol/m²)')
    axes[1, 1].set_title('Surface Excess Evolution')
    axes[1, 1].legend()
    axes[1, 1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('bubble_tensiometry_max_pressure_analysis.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    return data

def compare_different_concentrations():
    """
    Compare bubble tensiometry results for different surfactant concentrations.
    """
    print("\\n" + "="*60)
    print("Concentration Effect Analysis")
    print("="*60)
    
    concentrations = [0.001, 0.003, 0.005, 0.008]  # mol/m³ (1, 3, 5, 8 mM)
    colors = ['blue', 'green', 'red', 'orange']
    
    # Fixed parameters
    D = 5.8e-10
    Gamma_max = 3.4e-6
    K = 1200
    r0 = 0.5e-3
    growth_rate = 0.08e-3
    
    results = []
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    print(f"\\nComparing concentrations:")
    print(f"{'Conc (mM)':<10} {'Max P (Pa)':<12} {'γ at max P (mN/m)':<18} {'Final Γ (µmol/m²)':<18}")
    print("-" * 65)
    
    for i, c0 in enumerate(concentrations):
        # Create model
        model = ConvexWardTordaiModel(
            D=D, c0=c0, r0=r0, bubble_growth_rate=growth_rate,
            Gamma_max=Gamma_max, K=K
        )
        
        # Simulate
        t, gamma, Gamma, radius, pressure = model.simulate_bubble_tensiometry(8.0, 400)
        
        # Find maximum pressure
        max_idx = np.argmax(pressure)
        max_pressure = pressure[max_idx]
        gamma_at_max = gamma[max_idx]
        final_Gamma = Gamma[-1]
        
        results.append({
            'concentration': c0,
            'time': t,
            'gamma': gamma,
            'Gamma': Gamma,
            'pressure': pressure,
            'max_pressure': max_pressure,
            'gamma_at_max': gamma_at_max
        })
        
        print(f"{c0*1000:<10.1f} {max_pressure:<12.0f} {gamma_at_max*1000:<18.2f} {final_Gamma*1e6:<18.3f}")
        
        # Plot results
        label = f'{c0*1000:.1f} mM'
        
        # Pressure evolution
        axes[0, 0].plot(t, pressure, color=colors[i], linewidth=2, label=label)
        
        # Surface tension evolution
        axes[0, 1].plot(t, gamma*1000, color=colors[i], linewidth=2, label=label)
        
        # Surface excess evolution
        axes[1, 0].plot(t, Gamma*1e6, color=colors[i], linewidth=2, label=label)
        
        # Surface pressure evolution
        axes[1, 1].plot(t, (0.072 - gamma)*1000, color=colors[i], linewidth=2, label=label)
    
    # Format plots
    axes[0, 0].set_xlabel('Time (s)')
    axes[0, 0].set_ylabel('Bubble Pressure (Pa)')
    axes[0, 0].set_title('Pressure vs Time')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)
    
    axes[0, 1].set_xlabel('Time (s)')
    axes[0, 1].set_ylabel('Surface Tension (mN/m)')
    axes[0, 1].set_title('Surface Tension vs Time')
    axes[0, 1].legend()
    axes[0, 1].grid(True, alpha=0.3)
    
    axes[1, 0].set_xlabel('Time (s)')
    axes[1, 0].set_ylabel('Surface Excess (µmol/m²)')
    axes[1, 0].set_title('Surface Excess vs Time')
    axes[1, 0].legend()
    axes[1, 0].grid(True, alpha=0.3)
    
    axes[1, 1].set_xlabel('Time (s)')
    axes[1, 1].set_ylabel('Surface Pressure (mN/m)')
    axes[1, 1].set_title('Surface Pressure vs Time')
    axes[1, 1].legend()
    axes[1, 1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('bubble_tensiometry_concentration_comparison.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    return results

if __name__ == "__main__":
    print("Bubble Tensiometry Analysis with Convex Ward-Tordai Model")
    print("=" * 65)
    
    # Analysis 1: Maximum bubble pressure method
    max_pressure_data = analyze_maximum_bubble_pressure()
    
    # Analysis 2: Concentration effects
    concentration_results = compare_different_concentrations()
    
    print("\\n" + "="*65)
    print("Analysis completed successfully!")
    print("Generated plots:")
    print("  - bubble_tensiometry_max_pressure_analysis.png")
    print("  - bubble_tensiometry_concentration_comparison.png")
    print("\\nConvex Ward-Tordai model provides realistic bubble tensiometry predictions!")