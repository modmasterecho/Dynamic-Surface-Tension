"""
Ward Tordai Model - Numerical Solution for Integral Equations

This module provides numerical solutions for integral equations related to the Ward-Tordai model,
which is used for calculating the interface concentration of surfactant solutions.
"""

import numpy as np
import scipy.integrate as integrate
from scipy.optimize import fsolve
from scipy.special import erfc
import matplotlib.pyplot as plt

class WardTordaiModel:
    """
    Ward-Tordai model for diffusion-controlled adsorption kinetics.
    
    This class implements the numerical solution of the Ward-Tordai integral equation
    for calculating surface excess concentration as a function of time.
    """
    
    def __init__(self, D, c0, Gamma_max=None, K=None, R=8.314, T=298.15):
        """
        Initialize Ward-Tordai model parameters.
        
        Parameters:
        -----------
        D : float
            Diffusion coefficient (m²/s)
        c0 : float
            Bulk concentration (mol/m³)
        Gamma_max : float, optional
            Maximum surface excess for Langmuir isotherm (mol/m²)
        K : float, optional
            Langmuir adsorption constant (m³/mol)
        R : float
            Gas constant (J/(mol·K))
        T : float
            Temperature (K)
        """
        self.D = D
        self.c0 = c0
        self.Gamma_max = Gamma_max
        self.K = K
        self.R = R
        self.T = T
    
    def kernel_function(self, t, tau, c_tau):
        """
        Kernel function for the Ward-Tordai integral equation.
        
        The integrand is c(0,τ) / √(t-τ) for the Ward-Tordai equation:
        Γ(t) = √(D/π) * [c0*√t - ∫[0,t] c(0,τ)/√(t-τ) dτ]
        
        Parameters:
        -----------
        t : float
            Current time
        tau : float
            Integration variable (time)
        c_tau : float
            Surface concentration at time tau
            
        Returns:
        --------
        float
            Kernel function value
        """
        if t <= tau or t - tau <= 0:
            return 0.0
        return c_tau / np.sqrt(t - tau)
    
    def langmuir_isotherm(self, c):
        """
        Langmuir adsorption isotherm.
        
        Parameters:
        -----------
        c : float or array
            Surface concentration
            
        Returns:
        --------
        float or array
            Surface excess concentration
        """
        if self.Gamma_max is None or self.K is None:
            raise ValueError("Gamma_max and K must be set for Langmuir isotherm calculation")
        return self.Gamma_max * self.K * c / (1 + self.K * c)
    
    def surface_tension_from_gibbs(self, Gamma, gamma0=0.072):
        """
        Calculate surface tension using Gibbs equation.
        
        Parameters:
        -----------
        Gamma : float or array
            Surface excess concentration (mol/m²)
        gamma0 : float
            Surface tension of pure solvent (N/m)
            
        Returns:
        --------
        float or array
            Surface tension (N/m)
        """
        # Ensure Gamma is positive
        Gamma = np.maximum(Gamma, 0.0)
        
        # Gibbs equation: dγ = -Γ * R * T * d(ln c)
        if self.Gamma_max is None or self.K is None:
            # Simple approximation: γ = γ0 - R*T*Γ
            return gamma0 - self.R * self.T * Gamma
        else:
            # From Langmuir isotherm: Γ = Γ_max * K * c / (1 + K * c)
            # Solving for c: c = Γ / (K * (Γ_max - Γ))
            # Ensure we don't divide by zero or get negative concentrations
            Gamma_safe = np.minimum(Gamma, 0.99 * self.Gamma_max)
            c_surface = np.maximum(Gamma_safe / (self.K * (self.Gamma_max - Gamma_safe)), 1e-12)
            
            # Surface tension: γ = γ0 - R*T*Γ_max*ln(1 + K*c)
            surface_tension = gamma0 - self.R * self.T * self.Gamma_max * np.log(1 + self.K * c_surface)
            return np.maximum(surface_tension, 0.001)  # Ensure positive surface tension
    
    def solve_ward_tordai(self, t_max, n_points=1000, method='trapz'):
        """
        Solve the Ward-Tordai integral equation numerically.
        
        The Ward-Tordai equation: 
        Γ(t) = √(D/π) * [c0*√t - ∫[0,t] c(0,τ)/√(t-τ) dτ]
        Combined with adsorption isotherm: Γ = f(c(0,t))
        
        Parameters:
        -----------
        t_max : float
            Maximum time for simulation (s)
        n_points : int
            Number of time points
        method : str
            Integration method ('trapz' or 'simpson')
            
        Returns:
        --------
        tuple
            (time_array, Gamma_array, c_surface_array)
        """
        # Time array (avoid t=0 to prevent division by zero)
        t = np.linspace(1e-6, t_max, n_points)
        dt = t[1] - t[0]
        
        # Initialize arrays
        Gamma = np.zeros(n_points)
        c_surface = np.zeros(n_points)
        
        # Initial conditions
        c_surface[0] = self.c0
        Gamma[0] = 0.0
        
        # Numerical solution loop
        for i in range(1, n_points):
            current_time = t[i]
            
            # Define the equation to solve at current time
            def ward_tordai_equation(c_s):
                """
                Equation to solve: Ward-Tordai equation combined with isotherm
                """
                # Calculate the integral term
                if i > 1:
                    tau_array = t[:i]  # Previous time points
                    c_tau_array = c_surface[:i]  # Surface concentrations at previous times
                    
                    # Integrand: c(0,τ) / √(t-τ)
                    integrand = c_tau_array / np.sqrt(current_time - tau_array)
                    
                    # Numerical integration
                    if method == 'simpson' and len(tau_array) > 2:
                        integral_value = integrate.simpson(integrand, tau_array)
                    else:
                        integral_value = np.trapz(integrand, tau_array)
                else:
                    integral_value = 0.0
                
                # Ward-Tordai equation: Γ = √(D/π) * [c0*√t - integral]
                ward_tordai_gamma = np.sqrt(self.D / np.pi) * (
                    self.c0 * np.sqrt(current_time) - integral_value
                )
                
                # Adsorption isotherm: Γ = f(c_s)
                if self.Gamma_max is not None and self.K is not None:
                    # Langmuir isotherm
                    isotherm_gamma = self.langmuir_isotherm(c_s)
                else:
                    # Simple linear isotherm assumption: Γ ∝ c_s
                    # This is just for cases without specified isotherm
                    isotherm_gamma = ward_tordai_gamma  # Direct from Ward-Tordai
                    return 0.0  # No equation to solve
                
                # Return the difference (should be zero at solution)
                return ward_tordai_gamma - isotherm_gamma
            
            # Solve for surface concentration
            if self.Gamma_max is not None and self.K is not None:
                # Use root finding to solve the coupled equations
                try:
                    # Initial guess based on previous value
                    initial_guess = max(c_surface[i-1], 1e-10)
                    
                    # Solve the equation
                    solution = fsolve(ward_tordai_equation, initial_guess, full_output=True)
                    c_s_solution = solution[0][0]
                    
                    # Check if solution is valid
                    if solution[2] == 1 and c_s_solution > 0:
                        c_surface[i] = c_s_solution
                        Gamma[i] = self.langmuir_isotherm(c_s_solution)
                    else:
                        # Fallback: use diffusion-limited approximation
                        c_surface[i] = c_surface[i-1]
                        Gamma[i] = np.sqrt(self.D / np.pi) * self.c0 * np.sqrt(current_time)
                        
                except:
                    # Fallback for numerical issues
                    c_surface[i] = c_surface[i-1]
                    Gamma[i] = np.sqrt(self.D / np.pi) * self.c0 * np.sqrt(current_time)
            else:
                # Direct calculation without isotherm (diffusion-limited case)
                c_surface[i] = self.c0
                
                # Calculate integral
                if i > 1:
                    tau_array = t[:i]
                    c_tau_array = c_surface[:i]
                    integrand = c_tau_array / np.sqrt(current_time - tau_array)
                    
                    if method == 'simpson' and len(tau_array) > 2:
                        integral_value = integrate.simpson(integrand, tau_array)
                    else:
                        integral_value = np.trapz(integrand, tau_array)
                else:
                    integral_value = 0.0
                
                # Ward-Tordai equation
                Gamma[i] = max(0.0, np.sqrt(self.D / np.pi) * (
                    self.c0 * np.sqrt(current_time) - integral_value
                ))
        
        return t, Gamma, c_surface
    
    def _interpolate_concentration(self, tau, t_array, c_array):
        """Helper function to interpolate concentration at given time tau."""
        if tau <= 0 or len(t_array) == 0:
            return self.c0
        return np.interp(tau, t_array, c_array)
    
    def simulate_dynamic_surface_tension(self, t_max, n_points=1000, gamma0=0.072):
        """
        Simulate dynamic surface tension evolution.
        
        Parameters:
        -----------
        t_max : float
            Maximum time for simulation (s)
        n_points : int
            Number of time points
        gamma0 : float
            Surface tension of pure solvent (N/m)
            
        Returns:
        --------
        tuple
            (time_array, surface_tension_array, Gamma_array)
        """
        t, Gamma, c_surface = self.solve_ward_tordai(t_max, n_points)
        gamma = self.surface_tension_from_gibbs(Gamma, gamma0)
        
        return t, gamma, Gamma


def create_example_surfactant():
    """
    Create an example surfactant model with typical parameters.
    
    Returns:
    --------
    WardTordaiModel
        Configured model instance
    """
    # Typical parameters for a nonionic surfactant like C12E4
    D = 5e-10  # Diffusion coefficient (m²/s)
    c0 = 0.1   # Bulk concentration (mol/m³)
    Gamma_max = 3e-6  # Maximum surface excess (mol/m²)
    K = 1000   # Langmuir constant (m³/mol)
    
    return WardTordaiModel(D=D, c0=c0, Gamma_max=Gamma_max, K=K)


def plot_ward_tordai_results(t, gamma, Gamma, title="Ward-Tordai Model Results"):
    """
    Plot the results of Ward-Tordai simulation.
    
    Parameters:
    -----------
    t : array
        Time array (s)
    gamma : array
        Surface tension array (N/m)
    Gamma : array
        Surface excess array (mol/m²)
    title : str
        Plot title
    """
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
    
    # Surface tension vs time
    ax1.semilogx(t[1:], gamma[1:], 'b-', linewidth=2)
    ax1.set_xlabel('Time (s)')
    ax1.set_ylabel('Surface Tension (N/m)')
    ax1.set_title(f'{title} - Dynamic Surface Tension')
    ax1.grid(True, alpha=0.3)
    
    # Surface excess vs time
    ax2.semilogx(t[1:], Gamma[1:] * 1e6, 'r-', linewidth=2)  # Convert to µmol/m²
    ax2.set_xlabel('Time (s)')
    ax2.set_ylabel('Surface Excess (µmol/m²)')
    ax2.set_title(f'{title} - Surface Excess')
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.show()


def compare_with_analytical_solution(model, t_max=3600, n_points=1000):
    """
    Compare numerical solution with analytical approximations.
    
    Parameters:
    -----------
    model : WardTordaiModel
        Model instance
    t_max : float
        Maximum time (s)
    n_points : int
        Number of points
    """
    # Numerical solution
    t, gamma_num, Gamma_num = model.simulate_dynamic_surface_tension(t_max, n_points)
    
    # Analytical approximation for short times (diffusion-limited)
    # Γ ≈ √(D/π) * c0 * √t for short times
    Gamma_analytical = np.sqrt(model.D / np.pi) * model.c0 * np.sqrt(t)
    
    # Plot comparison
    plt.figure(figsize=(12, 5))
    
    plt.subplot(1, 2, 1)
    plt.loglog(t[1:], Gamma_num[1:] * 1e6, 'b-', label='Numerical', linewidth=2)
    plt.loglog(t[1:], Gamma_analytical[1:] * 1e6, 'r--', label='Analytical (short time)', linewidth=2)
    plt.xlabel('Time (s)')
    plt.ylabel('Surface Excess (µmol/m²)')
    plt.title('Surface Excess Comparison')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.subplot(1, 2, 2)
    plt.semilogx(t[1:], gamma_num[1:], 'b-', label='Numerical', linewidth=2)
    plt.xlabel('Time (s)')
    plt.ylabel('Surface Tension (N/m)')
    plt.title('Dynamic Surface Tension')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.show()


def demonstrate_ward_tordai_model():
    """
    Demonstrate the Ward-Tordai model with example calculations.
    """
    print("Ward-Tordai Model Demonstration")
    print("=" * 40)
    
    # Create model
    model = create_example_surfactant()
    
    print(f"Model parameters:")
    print(f"  Diffusion coefficient: {model.D:.2e} m²/s")
    print(f"  Bulk concentration: {model.c0} mol/m³")
    print(f"  Max surface excess: {model.Gamma_max:.2e} mol/m²")
    print(f"  Langmuir constant: {model.K} m³/mol")
    print()
    
    # Run simulation
    print("Running simulation...")
    t_max = 3600  # 1 hour
    n_points = 500
    
    t, gamma, Gamma = model.simulate_dynamic_surface_tension(t_max, n_points)
    
    # Print some results
    print(f"Results at selected times:")
    time_indices = [10, 50, 100, 200, -1]
    for i in time_indices:
        print(f"  t = {t[i]:8.1f} s: γ = {gamma[i]:.4f} N/m, Γ = {Gamma[i]*1e6:.2f} µmol/m²")
    
    # Plot results
    plot_ward_tordai_results(t, gamma, Gamma, "Example Surfactant")
    
    # Compare with analytical solution
    compare_with_analytical_solution(model, t_max, n_points)


if __name__ == "__main__":
    # Run demonstration
    demonstrate_ward_tordai_model()