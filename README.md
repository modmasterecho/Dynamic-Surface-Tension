# Dynamic-Surface-Tension

This Program aims at finding possible extensions for the Ward Tordai model for calculating the interface concentration of a solution.

## Overview

The Ward-Tordai equation is a fundamental model in colloid and interface science that describes the adsorption of surfactants at liquid interfaces. This project provides tools for:
- Numerical solution of the Ward-Tordai integral equation
- Fitting experimental surface tension data to theoretical models
- Analyzing dynamic surface tension measurements

## Project Structure

```
Dynamic-Surface-Tension/
├── Ward_Tordai.py           # Numerical solution module for integral equations
├── plotting_and_fitting.py  # Data visualization and model fitting
├── testing.ipynb            # Jupyter notebook for interactive testing
├── devlog.md                # Development log
├── suggestions.md           # Measurement device improvement ideas
└── README.md                # This file
```

## Variables and Parameters

### Surface Tension Variables

- **γ (gamma)**: Surface tension [mN/m or mJ/m²]
  - `gamma_0`: Initial surface tension of pure solvent (typically water: ~72 mN/m at 25°C)
  - `gamma_eq`: Equilibrium surface tension at long times
  - `gamma(t)`: Time-dependent surface tension

- **t**: Time [s]
  - Measurement time from the moment of fresh interface creation

- **Γ (Gamma)**: Surface concentration (surface excess) [mol/m²]
  - Amount of surfactant adsorbed per unit area at the interface

### Diffusion and Kinetic Parameters

- **D**: Diffusion coefficient [m²/s]
  - Describes the rate at which surfactant molecules diffuse to the interface

- **c₀ (c_0)**: Bulk concentration [mol/m³ or mol/L]
  - Initial concentration of surfactant in the bulk solution

- **c(x,t)**: Concentration profile as a function of distance x from interface and time t

### Model-Specific Parameters

- **k**: Rate constant [various units depending on model]
  - Empirical parameter in simplified models

- **R**: Universal gas constant [8.314 J/(mol·K)]

- **T**: Absolute temperature [K]

## Code Snippets

### Basic Usage Example

```python
import numpy as np
import Ward_Tordai
import plotting_and_fitting as pf

# Generate or load experimental data
time = np.array([1, 5, 10, 20, 50, 100])  # seconds
surface_tension = np.array([68, 60, 55, 50, 45, 42])  # mN/m

# Plot the data
pf.plot_surface_tension_vs_time(time, surface_tension)
```

### Fitting a Model

```python
# Define a model function (example: exponential decay)
def model_function(t, gamma_eq, gamma_0, k):
    """
    Simple exponential model for dynamic surface tension.
    
    Parameters:
    -----------
    t : float or array
        Time values
    gamma_eq : float
        Equilibrium surface tension
    gamma_0 : float
        Initial surface tension
    k : float
        Decay rate constant
    """
    return gamma_eq + (gamma_0 - gamma_eq) * np.exp(-k * np.sqrt(t))

# Initial parameter guesses
initial_params = [40.0, 72.0, 0.1]

# Fit the model
popt, pcov = pf.fit_model(time, surface_tension, model_function, initial_params)

# Generate fitted curve
time_fine = np.linspace(time.min(), time.max(), 100)
fitted_values = model_function(time_fine, *popt)

# Visualize
pf.plot_data_with_fit(time, surface_tension, model_function(time, *popt))
```

### Ward-Tordai Equation

The Ward-Tordai equation relates the surface concentration to the bulk concentration and diffusion:

```python
# Mathematical form (for reference):
# Γ(t) = 2 * c₀ * sqrt(D*t/π) - 2 * sqrt(D/π) * ∫₀ᵗ c(0,τ) / sqrt(t-τ) dτ
```

This equation must be solved numerically as it is an integral equation where the unknown function appears inside the integral.

## Installation

### Requirements

```bash
pip install numpy scipy matplotlib pandas jupyter
```

### Running the Code

1. Clone this repository
2. Install required packages
3. Open `testing.ipynb` in Jupyter Notebook to explore examples
4. Import modules in your Python scripts as needed

## Usage

### Interactive Testing

Launch Jupyter Notebook:
```bash
jupyter notebook testing.ipynb
```

### Python Scripts

```python
# Import the modules
import Ward_Tordai as wt
import plotting_and_fitting as pf

# Your analysis code here
```

## Key Equations

### Gibbs Adsorption Equation

Relates surface tension change to surface concentration:

```
dγ = -Γ RT d(ln c)
```

Where:
- dγ: Change in surface tension
- Γ: Surface excess concentration
- R: Gas constant (8.314 J/(mol·K))
- T: Absolute temperature
- c: Bulk concentration

### Diffusion-Controlled Adsorption

For purely diffusion-controlled adsorption:

```
Γ(t) = 2c₀√(Dt/π)
```

This is valid for short times when the subsurface concentration hasn't been depleted.

## Contributing

See `suggestions.md` for ideas on improving the measurement device and analysis methods.

## Development Log

Check `devlog.md` for development history and planned features.

## License

[Specify license here]

## References

1. Ward, A.F.H., Tordai, L. (1946). "Time‐Dependence of Boundary Tensions of Solutions I. The Role of Diffusion in Time‐Effects." The Journal of Chemical Physics, 14(7), 453-461.

## Contact

[Add contact information]
