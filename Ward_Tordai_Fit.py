from __future__ import annotations
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit, minimize
import math
import sys
from typing import Callable, Tuple, Optional
from dataclasses import dataclass


PI = np.pi
R = 8.314  # Gas constant (J/(mol*K))


@dataclass
class WardTordaiParams:
    """
    Data class to hold parameters for the Ward-Tordai model with various isotherms.
    """
    # Global parameters
    Cb: float = 1.0e-2        # Bulk concentration (mol/m^3)
    D: float = 8.8e-10        # Diffusion coefficient (m^2/s)
    rb: float = 1.12e-3       # Radius of the bubble (m)
    Tmpr: float = 298.15      # Temperature (K)
    gamma_0: float = 72.2e-3  # Solvent surface tension (N/m)

    # Isotherm constants
    gamma_m: float = 3.91e-6  # Maximum surface concentration (mol/m^2)
    kh: float = 2146.0        # Henry constant (m)
    kl: float = 2146.0        # Langmuir constant (m^3/mol)
    kf: float = 0.0           # Frumkin constant (m^3/mol)
    A: float = 0.0            # Frumkin interaction parameter
    kfl: float = 0.0          # Freundlich constant (mol^x m^y)
    knl: float = 1.0          # Freundlich consistency index
    kv: float = 0.0           # Volmer constant (m^3/mol)
    nn: float = 1.0           # Value of n (-)

    # Time stepping controls
    T: float = 3000.0
    h: float = 1.0
    err: float = 1e-30
    M: int = 100

    # Selections
    isotherm: int = 1   # 1=Henry, 2=Langmuir, 3=Frumkin, 4=Freundlich, 5=Volmer
    geometry: int = 0   # 0=spherical, 1=planar


# Isotherm helper functions

def henry(gamma: float, p: WardTordaiParams) -> float:
    """Henry isotherm: gamma/kh"""
    if p.kh == 0:
        return float('inf')
    return gamma / p.kh


def langmuir(gamma: float, p: WardTordaiParams) -> float:
    """Langmuir isotherm: gamma / (kl * (gamma_m - gamma))"""
    denom = (p.gamma_m - gamma)
    if p.kl == 0 or denom == 0:
        return float('inf')
    return gamma / (p.kl * denom)


def frumkin(gamma: float, p: WardTordaiParams) -> float:
    """Frumkin isotherm"""
    denom = (p.gamma_m - gamma)
    if p.kf == 0 or denom == 0 or p.gamma_m == 0:
        return float('inf')
    return (1.0 / p.kf) * (gamma / denom) * math.exp(-p.A * (gamma / p.gamma_m))


def freundlich(gamma: float, p: WardTordaiParams) -> float:
    """Freundlich isotherm: (gamma/kfl)^knl"""
    if p.kfl == 0:
        return float('inf')
    if gamma < 0:
        return float('nan')
    return (gamma / p.kfl) ** p.knl


def volmer(gamma: float, p: WardTordaiParams) -> float:
    """Volmer isotherm"""
    denom = (p.gamma_m - gamma)
    if p.kv == 0 or denom == 0:
        return float('inf')
    frac = gamma / denom
    return p.kv * frac * math.exp(frac)


def geom(geometry: int, t: float, p: WardTordaiParams) -> float:
    """Geometry term for Ward-Tordai equation"""
    if t < 0:
        return 0.0
    
    if geometry == 0:  # spherical (convex)
        rt = max(t, 0.0)
        core = 2.0 * math.sqrt(rt * p.D / PI) + (p.D / p.rb) * rt if p.rb != 0 else float('inf')
        return core * p.Cb
    elif geometry == 1:  # planar
        return math.sqrt(p.D / PI) * 2.0 * p.Cb * math.sqrt(max(t, 0.0))
    else:
        raise ValueError("Geometry must be 0 or 1")


def stn(isotherm: int, nn: float, gamma: float, p: WardTordaiParams) -> float:
    """Surface excess to surface tension change via isotherm"""
    if isotherm == 1:  # Henry
        return -nn * R * p.Tmpr * gamma
    elif isotherm == 2:  # Langmuir
        if p.gamma_m == 0:
            return float('inf')
        return nn * R * p.Tmpr * p.gamma_m * math.log(1.0 - gamma / p.gamma_m)
    elif isotherm == 3:  # Frumkin
        if p.gamma_m == 0:
            return float('nan')
        x = gamma / p.gamma_m
        return nn * R * p.Tmpr * p.gamma_m * math.log(1.0 - x) + R * nn * p.Tmpr * p.A / 2.0 * p.gamma_m * x * x
    elif isotherm == 4:  # Freundlich
        return nn * p.knl * R * p.Tmpr * gamma
    elif isotherm == 5:  # Volmer
        if p.gamma_m == 0:
            return float('nan')
        return nn * p.gamma_m * p.gamma_m / (p.gamma_m - gamma) * R * p.Tmpr
    else:
        raise ValueError("Isotherm must be 1..5")


def K(t: float, tau: float, gamma: float, p: WardTordaiParams) -> float:
    """Kernel function for Ward-Tordai integral"""
    if t <= tau:
        return 0.0
    
    sqrt_term = math.sqrt(p.D / PI) / math.sqrt(max(t - tau, 0.0))
    add_spherical = p.D / p.rb if (p.geometry == 0 and p.rb != 0) else 0.0
    
    if p.geometry == 0:  # spherical
        factor = -(sqrt_term + add_spherical)
    elif p.geometry == 1:  # planar
        factor = -sqrt_term
    else:
        raise ValueError("Geometry must be 0 or 1")
    
    # Select isotherm function
    if p.isotherm == 1:
        return factor * henry(gamma, p)
    elif p.isotherm == 2:
        return factor * langmuir(gamma, p)
    elif p.isotherm == 3:
        return factor * frumkin(gamma, p)
    elif p.isotherm == 4:
        return factor * freundlich(gamma, p)
    elif p.isotherm == 5:
        return factor * volmer(gamma, p)
    else:
        raise ValueError("Isotherm must be 1..5")


def rtbis(tn: float, x1: float, x2: float, p: WardTordaiParams) -> Tuple[float, bool]:
    """Root Bisection Method to solve the Ward-Tordai equation for gamma at time tn"""
    def f(x: float) -> float:
        return x - p.h / 2.0 * K(tn, 0.9999 * tn, x, p) - x2
    
    dx = 0.0
    f1 = f(x1)
    
    if f1 < 0.0:
        dx = x2 - x1
        rtb = x1
    else:
        dx = x1 - x2
        rtb = x2
    
    acc = False
    for _ in range(p.M):
        dx *= 0.5
        xmid = rtb + dx
        fmid = f(xmid)
        if fmid < 0.0:
            rtb = xmid
        if abs(dx) < p.err or fmid == 0.0:
            acc = True
            return rtb, acc
    
    return rtb, acc


def ward_tordai(p: WardTordaiParams) -> Tuple[list[float], list[float], list[float], bool]:
    """Run Ward-Tordai simulation"""
    if p.h <= 0:
        raise ValueError("Time step h must be > 0")
    
    N = int(math.ceil(p.T / p.h))
    t = [0.0] * (N + 1)
    gamma = [0.0] * (N + 1)
    st = [p.gamma_0] * (N + 1)
    
    for i in range(N + 1):
        t[i] = p.h * i
    
    acc = True
    
    for n in range(1, N + 1):
        ssum = 0.0
        for j in range(1, n):
            ssum += K(t[n], t[j], gamma[j], p)
        
        x = geom(p.geometry, t[n], p) + p.h * ssum
        root, ok = rtbis(t[n], gamma[n - 1], x, p)
        gamma[n] = root
        st[n] = p.gamma_0 + stn(p.isotherm, p.nn, gamma[n], p)
        acc = acc and ok
    
    return t, gamma, st, acc


# ==================== Fitting Functions ====================

def fit_ward_tordai(t_data: np.ndarray, st_data: np.ndarray, p: WardTordaiParams, 
                    fit_params: list[str]) -> Tuple[WardTordaiParams, dict]:
    """
    Fit Ward-Tordai model to experimental data.
    
    Parameters:
    -----------
    t_data : array of time points (s)
    st_data : array of surface tension values (N/m)
    p : initial parameter guess
    fit_params : list of parameter names to fit, e.g., ['D', 'kl', 'gamma_m']
    
    Returns:
    --------
    fitted_params : WardTordaiParams with optimized values
    fit_info : dict with fitting statistics
    """
    
    # Create parameter vector from fit_params
    param_names = fit_params
    x0 = [getattr(p, name) for name in param_names]
    
    def objective(x_vec):
        # Update parameters
        p_temp = WardTordaiParams(**vars(p))
        for name, val in zip(param_names, x_vec):
            setattr(p_temp, name, val)
        
        # Set simulation time range to match data
        p_temp.T = float(np.max(t_data))
        p_temp.h = float(np.min(np.diff(t_data))) if len(t_data) > 1 else 1.0
        
        # Run simulation
        try:
            t_sim, _, st_sim, _ = ward_tordai(p_temp)
            
            # Interpolate to data points
            st_interp = np.interp(t_data, t_sim, st_sim)
            
            # Calculate residuals
            residuals = st_data - st_interp
            return np.sum(residuals**2)
        except:
            return 1e10  # Large penalty for failed simulation
    
    # Optimize
    result = minimize(objective, x0, method='Nelder-Mead')
    
    # Update parameters with fitted values
    p_fitted = WardTordaiParams(**vars(p))
    for name, val in zip(param_names, result.x):
        setattr(p_fitted, name, val)
    
    # Calculate final fit statistics
    p_fitted.T = float(np.max(t_data))
    p_fitted.h = float(np.min(np.diff(t_data))) if len(t_data) > 1 else 1.0
    t_sim, _, st_sim, _ = ward_tordai(p_fitted)
    st_interp = np.interp(t_data, t_sim, st_sim)
    
    residuals = st_data - st_interp
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((st_data - np.mean(st_data))**2)
    r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
    
    fit_info = {
        'r_squared': r_squared,
        'residual_sum_squares': ss_res,
        'fitted_params': {name: getattr(p_fitted, name) for name in param_names},
        'success': result.success
    }
    
    return p_fitted, fit_info


# ==================== Interactive Functions ====================

def safe_int(prompt: str, valid: Tuple[int, ...] | None = None) -> int:
    while True:
        try:
            v = int(input(prompt).strip())
            if valid and v not in valid:
                print(f"Please enter one of: {valid}")
                continue
            return v
        except ValueError:
            print("Invalid integer. Try again.")


def safe_float(prompt: str) -> float:
    while True:
        try:
            return float(input(prompt).strip())
        except ValueError:
            print("Invalid number. Use e.g. 1e-6 for scientific notation.")


def banner():
    print("         Diffusion-controlled adsorption rate and DST calculator")
    print("-------------------------------------------------------------------------------")
    print("This code is a Python translation of the program by X.B. Li and P. Stevenson,\n"
          "University of Newcastle, Australia. If used for publication, please cite their\n"
          "contribution. Any questions: paul.stevenson@newcastle.edu.au")
    print("-------------------------------------------------------------------------------")
    print("DISCLAIMER")
    print("This program is provided as is for educational and indicative purposes only.\n"
          "No guarantee is provided for accuracy or applicability. Use at your own risk.")
    print("-------------------------------------------------------------------------------\n")


def explain_models():
    print("-------------------------------------------------------------------------------")
    print("Planar adsorption uses Ward & Tordai (J Chem Phys 14, 453-461, 1946)")
    print("Spherical/convex uses Lin et al. (AIChE J 36, 1785-1795, 1990)")
    print("Isotherms per Eastoe & Dalton (Adv Colloid Interface Sci 84, 103-144, 2000)")
    print("-------------------------------------------------------------------------------\n")


def collect_parameters(p: WardTordaiParams) -> None:
    # Advanced mode
    option = input("Advanced mode to change tolerance and max iterations? (Y/N): ").strip().lower()
    if option == 'y':
        p.err = safe_float("Tolerance (default 1e-30): ")
        p.M = safe_int("Maximum iterations (default 100): ")

    print("\n-------------------------------------------------------------------------------")
    print("Enter model parameters. Use suggested units and scientific notation (e.g., 1e-6).")

    p.Tmpr = safe_float("Thermodynamic temperature, T (K): ")
    p.gamma_0 = safe_float("Solvent surface tension, sigma0 (N/m): ")
    p.Cb = safe_float("Bulk concentration, cb (mol/m^3): ")
    p.D = safe_float("Molecular diffusivity, D (m^2/s): ")
    p.nn = safe_float("Value of n (-): ")
    p.T = safe_float("Total simulation time, T (s): ")
    p.h = safe_float("Time increment, h (s): ")

    print("\n-------------------------------------------------------------------------------")
    p.geometry = safe_int("Interface shape: 0 = spherical (convex), 1 = planar: ", valid=(0, 1))

    print("\n-------------------------------------------------------------------------------")
    print("Adsorption isotherm:")
    print(" 1 = Henry's Law")
    print(" 2 = Langmuir")
    print(" 3 = Frumkin")
    print(" 4 = Freundlich")
    print(" 5 = Volmer")
    p.isotherm = safe_int("Choose 1..5: ", valid=(1, 2, 3, 4, 5))

    if p.geometry == 0:
        p.rb = safe_float("Radius of curvature, rb (m): ")

    if p.isotherm == 1:  # Henry
        p.kh = safe_float("Henry's Law constant, kh (m): ")
        p.gamma_m = safe_float("Saturation surface excess, gamma_m (mol/m^2): ")
    elif p.isotherm == 2:  # Langmuir
        p.kl = safe_float("Langmuir constant, kl (m^3/mol): ")
        p.gamma_m = safe_float("Saturation surface excess, gamma_m (mol/m^2): ")
    elif p.isotherm == 3:  # Frumkin
        p.kf = safe_float("Frumkin constant, kf (m^3/mol): ")
        p.A = safe_float("Surface interaction parameter, A (-): ")
        p.gamma_m = safe_float("Saturation surface excess, gamma_m (mol/m^2): ")
    elif p.isotherm == 4:  # Freundlich
        p.kfl = safe_float("Freundlich constant, kfl (mol^x m^y): ")
        p.knl = safe_float("Consistency index, knl (-): ")
    elif p.isotherm == 5:  # Volmer
        p.kv = safe_float("Volmer constant, kv (m^3/mol): ")
        p.gamma_m = safe_float("Saturation surface excess, gamma_m (mol/m^2): ")

    print("\nPlease check the values you entered above.\n")


def write_results(path: str, t: list[float], gamma: list[float], st: list[float]) -> None:
    with open(path, 'w', encoding='utf-8') as f:
        f.write("time_s\tsurface_excess_mol_m2\tdynamic_surface_tension_N_m\n")
        for ti, gi, si in zip(t, gamma, st):
            f.write(f"{ti}\t{gi}\t{si}\n")


def fit_mode():
    """Interactive fitting mode"""
    print("\n===============================================================================")
    print("                            DATA FITTING MODE")
    print("===============================================================================\n")
    
    # Load data
    file_path = input("Enter path to data file (CSV or TSV with columns: time, surface_tension): ").strip()
    try:
        if file_path.endswith('.csv'):
            df = pd.read_csv(file_path)
        else:
            df = pd.read_csv(file_path, sep='\t')
        
        # Assume first column is time, second is surface tension
        t_data = df.iloc[:, 0].values
        st_data = df.iloc[:, 1].values
        print(f"Loaded {len(t_data)} data points from {file_path}")
    except Exception as e:
        print(f"Error loading file: {e}")
        return
    
    # Setup initial parameters
    p = WardTordaiParams()
    collect_parameters(p)
    
    # Select parameters to fit
    print("\n-------------------------------------------------------------------------------")
    print("Select parameters to fit (comma-separated):")
    print("Available: D, kh, kl, kf, A, gamma_m, kfl, knl, kv, Cb, nn")
    fit_params_str = input("Enter parameters: ").strip()
    fit_params = [s.strip() for s in fit_params_str.split(',')]
    
    print(f"\nFitting parameters: {fit_params}")
    print("This may take a while...")
    
    # Perform fit
    p_fitted, fit_info = fit_ward_tordai(t_data, st_data, p, fit_params)
    
    # Display results
    print("\n===============================================================================")
    print("                            FITTING RESULTS")
    print("===============================================================================")
    print(f"R-squared: {fit_info['r_squared']:.6f}")
    print(f"Residual sum of squares: {fit_info['residual_sum_squares']:.6e}")
    print(f"Success: {fit_info['success']}")
    print("\nFitted parameters:")
    for name, value in fit_info['fitted_params'].items():
        print(f"  {name}: {value:.6e}")
    
    # Plot comparison
    p_fitted.T = float(np.max(t_data))
    p_fitted.h = float(np.min(np.diff(t_data))) if len(t_data) > 1 else 1.0
    t_sim, _, st_sim, _ = ward_tordai(p_fitted)
    
    plt.figure(figsize=(10, 6))
    plt.plot(t_data, st_data, 'o', label='Experimental data', markersize=5)
    plt.plot(t_sim, st_sim, '-', label='Fitted model', linewidth=2)
    plt.xlabel('Time (s)')
    plt.ylabel('Surface Tension (N/m)')
    plt.legend()
    plt.grid(True)
    plt.title(f'Ward-Tordai Fit (R² = {fit_info["r_squared"]:.4f})')
    plt.savefig('fit_result.png', dpi=150, bbox_inches='tight')
    plt.show()
    print("\nPlot saved as 'fit_result.png'")
    
    # Save fitted results
    write_results('Fitted_Simulation.tsv', t_sim, _, st_sim)
    print("Fitted simulation saved as 'Fitted_Simulation.tsv'")


def simulation_mode():
    """Interactive simulation mode"""
    while True:
        p = WardTordaiParams()
        collect_parameters(p)

        t, gamma, st, acc = ward_tordai(p)
        out_path = "Simulation.tsv"
        write_results(out_path, t, gamma, st)

        if not acc:
            print("-------------------------------------------------------------------------------")
            print(f"WARNING: Required accuracy {p.err} was not reached in {p.M} iterations!")
            print("Try advanced mode to adjust numerical settings if results look unsatisfactory.")

        print("-------------------------------------------------------------------------------")
        print("A file 'Simulation.tsv' containing time, surface excess, and DST has been written\n"
              "to the current directory. Save or rename it before running another simulation\n"
              "or it will be overwritten.")
        
        again = input("Continue simulations? (Y/N): ").strip().lower()
        if again != 'y':
            break


def main(argv: list[str]) -> int:
    banner()
    explain_models()
    
    print("===============================================================================")
    print("Select mode:")
    print("  1 = Simulation mode (generate DST curves)")
    print("  2 = Fitting mode (fit model to experimental data)")
    mode = safe_int("Choose 1 or 2: ", valid=(1, 2))
    
    if mode == 1:
        simulation_mode()
    elif mode == 2:
        fit_mode()
    
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
