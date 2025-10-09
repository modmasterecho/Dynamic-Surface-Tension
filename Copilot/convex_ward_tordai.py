"""
Convex Ward-Tordai Model for Bubble Tensiometry

This module provides the convex (spherical) version of the Ward-Tordai model,
specifically designed for bubble tensiometric measurements where the interface
area changes with time due to bubble growth.
"""

import numpy as np
import scipy.integrate as integrate
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from Ward_Tordai import WardTordaiModel

class ConvexWardTordaiModel(WardTordaiModel):
    """
    Convex Ward-Tordai model for bubble tensiometry.
    
    This class extends the planar Ward-Tordai model to account for the convex
    (spherical) geometry of growing bubbles in tensiometric measurements.
    """
    
    def __init__(self, D, c0, r0=1e-4, bubble_growth_rate=1e-5, 
                 Gamma_max=None, K=None, R=8.314, T=298.15):
        """
        Initialize convex Ward-Tordai model parameters.
        
        Parameters:
        -----------
        D : float
            Diffusion coefficient (m²/s)
        c0 : float
            Bulk concentration (mol/m³)
        r0 : float
            Initial bubble radius (m)
        bubble_growth_rate : float
            Bubble growth rate (m/s) - constant growth approximation
        Gamma_max : float, optional
            Maximum surface excess for Langmuir isotherm (mol/m²)
        K : float, optional
            Langmuir adsorption constant (m³/mol)
        R : float
            Gas constant (J/(mol·K))
        T : float
            Temperature (K)
        """
        super().__init__(D, c0, Gamma_max, K, R, T)
        self.r0 = r0
        self.bubble_growth_rate = bubble_growth_rate
    
    def bubble_radius(self, t):
        """
        Calculate bubble radius at time t.
        
        Parameters:
        -----------
        t : float or array
            Time (s)
            
        Returns:
        --------
        float or array
            Bubble radius (m)
        """
        return self.r0 + self.bubble_growth_rate * t
    
    def bubble_surface_area(self, t):
        """
        Calculate bubble surface area at time t.
        
        Parameters:
        -----------
        t : float or array
            Time (s)
            
        Returns:
        --------
        float or array
            Surface area (m²)
        """
        r = self.bubble_radius(t)
        return 4 * np.pi * r**2
    
    def convex_kernel_function(self, t, tau, c_tau):
        """
        Kernel function for the convex Ward-Tordai integral equation.
        
        For spherical geometry, the kernel is modified to account for
        the curved interface and changing surface area.
        
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
            Convex kernel function value
        """
        if t <= tau or t - tau <= 0:
            return 0.0
        
        r_t = self.bubble_radius(t)
        r_tau = self.bubble_radius(tau)
        
        # Area ratio correction for expanding interface
        area_ratio = (r_tau / r_t)**2
        
        # Modified kernel for spherical geometry
        kernel = area_ratio * c_tau / np.sqrt(t - tau)
        
        # Additional correction for curvature effects
        curvature_correction = 1 + self.D / (r_t * np.sqrt(np.pi * self.D * (t - tau)))
        
        return kernel * curvature_correction
    
    def solve_convex_ward_tordai(self, t_max, n_points=1000, method='trapz'):
        """
        Solve the convex Ward-Tordai integral equation numerically.
        
        The convex Ward-Tordai equation accounts for spherical diffusion
        and changing interface area:
        
        Γ(t) = √(D/π) * [c0*√t - ∫[0,t] K_convex(t,τ,c(τ)) dτ]
        
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
            (time_array, Gamma_array, c_surface_array, radius_array)
        """
        # Time array (avoid t=0 to prevent division by zero)
        t = np.linspace(1e-6, t_max, n_points)
        
        # Initialize arrays
        Gamma = np.zeros(n_points)
        c_surface = np.zeros(n_points)
        radius = self.bubble_radius(t)
        
        # Initial conditions
        c_surface[0] = self.c0
        Gamma[0] = 0.0
        
        # Numerical solution loop
        for i in range(1, n_points):
            current_time = t[i]
            current_radius = radius[i]
            
            # Define the equation to solve at current time
            def convex_ward_tordai_equation(c_s):
                """
                Convex Ward-Tordai equation combined with isotherm
                """
                # Calculate the integral term
                if i > 1:
                    tau_array = t[:i]  # Previous time points
                    c_tau_array = c_surface[:i]  # Surface concentrations
                    
                    # Convex integrand
                    integrand = np.array([
                        self.convex_kernel_function(current_time, tau, c_tau)
                        for tau, c_tau in zip(tau_array, c_tau_array)
                    ])
                    
                    # Numerical integration
                    if method == 'simpson' and len(tau_array) > 2:
                        integral_value = integrate.simpson(integrand, tau_array)
                    else:
                        integral_value = np.trapz(integrand, tau_array)
                else:
                    integral_value = 0.0
                
                # Convex Ward-Tordai equation with curvature correction
                curvature_factor = 1 + np.sqrt(self.D / (np.pi * current_radius**2))
                
                ward_tordai_gamma = np.sqrt(self.D / np.pi) * curvature_factor * (
                    self.c0 * np.sqrt(current_time) - integral_value
                )
                
                # Adsorption isotherm
                if self.Gamma_max is not None and self.K is not None:
                    isotherm_gamma = self.langmuir_isotherm(c_s)
                else:
                    isotherm_gamma = ward_tordai_gamma
                    return 0.0
                
                return ward_tordai_gamma - isotherm_gamma
            
            # Solve for surface concentration
            if self.Gamma_max is not None and self.K is not None:
                try:
                    initial_guess = max(c_surface[i-1], 1e-10)
                    solution = fsolve(convex_ward_tordai_equation, initial_guess, full_output=True)
                    c_s_solution = solution[0][0]
                    
                    if solution[2] == 1 and c_s_solution > 0:
                        c_surface[i] = c_s_solution
                        Gamma[i] = self.langmuir_isotherm(c_s_solution)
                    else:
                        # Fallback
                        c_surface[i] = c_surface[i-1]
                        curvature_factor = 1 + np.sqrt(self.D / (np.pi * current_radius**2))
                        Gamma[i] = np.sqrt(self.D / np.pi) * curvature_factor * self.c0 * np.sqrt(current_time)
                        
                except:
                    c_surface[i] = c_surface[i-1]
                    curvature_factor = 1 + np.sqrt(self.D / (np.pi * current_radius**2))
                    Gamma[i] = np.sqrt(self.D / np.pi) * curvature_factor * self.c0 * np.sqrt(current_time)
            else:
                # Direct calculation without isotherm
                c_surface[i] = self.c0
                
                # Calculate integral
                if i > 1:
                    tau_array = t[:i]
                    c_tau_array = c_surface[:i]
                    integrand = np.array([
                        self.convex_kernel_function(current_time, tau, c_tau)
                        for tau, c_tau in zip(tau_array, c_tau_array)
                    ])
                    
                    if method == 'simpson' and len(tau_array) > 2:
                        integral_value = integrate.simpson(integrand, tau_array)
                    else:
                        integral_value = np.trapz(integrand, tau_array)
                else:
                    integral_value = 0.0
                
                # Convex Ward-Tordai equation
                curvature_factor = 1 + np.sqrt(self.D / (np.pi * current_radius**2))
                Gamma[i] = max(0.0, np.sqrt(self.D / np.pi) * curvature_factor * (
                    self.c0 * np.sqrt(current_time) - integral_value
                ))
        
        return t, Gamma, c_surface, radius
    
    def calculate_bubble_pressure(self, t, gamma, apply_correction=True):
        """
        Calculate bubble pressure using Young-Laplace equation.
        
        For a spherical bubble: ΔP = 2γ/R
        
        Parameters:
        -----------
        t : array
            Time array (s)
        gamma : array
            Surface tension array (N/m)
        apply_correction : bool
            Apply correction for dynamic effects
            
        Returns:
        --------
        array
            Pressure difference (Pa)
        """
        radius = self.bubble_radius(t)
        pressure = 2 * gamma / radius
        
        if apply_correction:
            # Dynamic pressure correction for growing bubble
            dr_dt = self.bubble_growth_rate
            dynamic_correction = 4 * np.pi * radius * dr_dt * gamma / self.bubble_surface_area(t)
            pressure += dynamic_correction
        
        return pressure
    
    def simulate_bubble_tensiometry(self, t_max, n_points=1000, gamma0=0.072):
        """
        Simulate complete bubble tensiometry measurement.
        
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
            (time, surface_tension, surface_excess, radius, pressure)
        """
        t, Gamma, c_surface, radius = self.solve_convex_ward_tordai(t_max, n_points)
        gamma = self.surface_tension_from_gibbs(Gamma, gamma0)
        pressure = self.calculate_bubble_pressure(t, gamma)
        
        return t, gamma, Gamma, radius, pressure


def create_bubble_tensiometry_example():
    """
    Create example bubble tensiometry model with realistic parameters.
    
    Returns:
    --------
    ConvexWardTordaiModel
        Configured model instance for bubble measurements
    """
    # Typical parameters for bubble tensiometry
    D = 3e-10          # Diffusion coefficient (m²/s)
    c0 = 0.005         # Bulk concentration (mol/m³) - 5 mM
    r0 = 0.5e-3        # Initial radius (0.5 mm capillary)
    growth_rate = 1e-4 # Growth rate (0.1 mm/s)
    Gamma_max = 2.5e-6 # Maximum surface excess (mol/m²)
    K = 800            # Langmuir constant (m³/mol)
    
    return ConvexWardTordaiModel(
        D=D, c0=c0, r0=r0, 
        bubble_growth_rate=growth_rate,
        Gamma_max=Gamma_max, K=K
    )


def generate_example_bubble_data():
    """
    Generate example bubble tensiometry data that mimics real measurements.
    
    Returns:
    --------
    dict
        Dictionary with time, pressure, and other measurement data
    """
    # Simulate a typical bubble tensiometry measurement
    model = create_bubble_tensiometry_example()
    
    # Typical measurement time for one bubble cycle
    t_max = 5.0  # 5 seconds
    n_points = 200
    
    t, gamma, Gamma, radius, pressure = model.simulate_bubble_tensiometry(t_max, n_points)
    
    # Add some realistic noise to the data
    np.random.seed(42)  # For reproducible results
    pressure_noise = np.random.normal(0, pressure.max() * 0.02, len(pressure))  # 2% noise
    gamma_noise = np.random.normal(0, 0.001, len(gamma))  # 1 mN/m noise
    
    return {
        'time': t,
        'surface_tension': gamma + gamma_noise,
        'surface_excess': Gamma,
        'radius': radius,
        'pressure': pressure + pressure_noise,
        'clean_pressure': pressure,
        'clean_surface_tension': gamma,
        'model_parameters': {
            'D': model.D,
            'c0': model.c0,
            'r0': model.r0,
            'growth_rate': model.bubble_growth_rate,
            'Gamma_max': model.Gamma_max,
            'K': model.K
        }
    }