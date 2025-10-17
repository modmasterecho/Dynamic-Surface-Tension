# Ward-Tordai Model Wrapper for scipy.optimize.curve_fit

This module provides a simplified interface to fit the Ward-Tordai dynamic surface tension model to experimental data using `scipy.optimize.curve_fit` or other optimization routines.

## Overview

The `ward_tordai_wrapper.py` module makes it easy to:
- Fit Ward-Tordai model parameters to experimental data
- Use the model with `scipy.optimize.curve_fit` 
- Evaluate the model at specific time points
- Choose which parameters to fit and which to hold constant

## Quick Start

### Basic Usage with curve_fit

```python
import numpy as np
from scipy.optimize import curve_fit
from ward_tordai_wrapper import create_fit_function

# Load your experimental data
t_data = np.array([...])   # Time in seconds
st_data = np.array([...])  # Surface tension in N/m

# Create a fit function with fixed parameters
fit_func = create_fit_function(
    fit_params=['D', 'kl'],      # Parameters to fit
    Cb=1.0e-2,                   # Bulk concentration (mol/m^3) - FIXED
    Tmpr=298.15,                 # Temperature (K) - FIXED
    gamma_0=72.2e-3,             # Solvent surface tension (N/m) - FIXED
    rb=1.12e-3,                  # Bubble radius (m) - FIXED
    nn=1.0,                      # Value of n - FIXED
    isotherm=2,                  # Langmuir isotherm
    geometry=0,                  # Spherical geometry
    gamma_m=3.91e-6              # Max surface concentration (mol/m^2) - FIXED
)

# Initial guesses for D and kl
p0 = [8.8e-10, 2146.0]

# Fit the model
popt, pcov = curve_fit(fit_func, t_data, st_data, p0=p0)

# Extract fitted parameters
D_fit, kl_fit = popt
errors = np.sqrt(np.diag(pcov))

print(f"Fitted D: {D_fit:.3e} ± {errors[0]:.3e} m^2/s")
print(f"Fitted kl: {kl_fit:.3e} ± {errors[1]:.3e} m^3/mol")
```

### Using Convenience Functions

For common isotherms, you can use pre-defined convenience functions:

```python
from scipy.optimize import curve_fit
from ward_tordai_wrapper import langmuir_fit_function

# Fit using the Langmuir convenience function
popt, pcov = curve_fit(
    lambda t, D, kl: langmuir_fit_function(t, D, kl, Cb=1.0e-2),
    t_data, 
    st_data, 
    p0=[8.8e-10, 2146.0]
)

D_fit, kl_fit = popt
```

### Direct Model Evaluation

You can also evaluate the model directly without fitting:

```python
from ward_tordai_wrapper import ward_tordai_model
import numpy as np

t = np.linspace(0, 1000, 100)
st = ward_tordai_model(
    t,
    D=8.8e-10,
    kl=2146.0,
    Cb=1.0e-2,
    gamma_m=3.91e-6,
    isotherm=2,
    geometry=0
)
```

## Available Functions

### Main Functions

1. **`ward_tordai_model(t, **params)`**
   - Direct evaluation of the Ward-Tordai model
   - Returns surface tension at time(s) t
   - All parameters can be specified as keyword arguments

2. **`create_fit_function(fit_params, **fixed_params)`**
   - Creates a fitting function for use with `curve_fit`
   - Specify which parameters to fit and which to hold constant
   - Returns a function with signature `func(t, *fit_params)`

3. **`create_simple_fit_function(isotherm_type, geometry_type, **params)`**
   - Simplified version using string names for isotherms
   - isotherm_type: 'henry', 'langmuir', 'frumkin', 'freundlich', 'volmer'
   - geometry_type: 'spherical', 'planar'

### Convenience Functions

Pre-configured functions for common isotherms:
- `langmuir_fit_function(t, D, kl, ...)`
- `henry_fit_function(t, D, kh, ...)`
- `frumkin_fit_function(t, D, kf, A, ...)`

## Parameters

### Required Parameters (Isotherm-Specific)

Depending on the isotherm type, you need to provide specific constants:

| Isotherm | Code | Required Constants |
|----------|------|-------------------|
| Henry | 1 | `kh` (Henry constant, m) |
| Langmuir | 2 | `kl` (Langmuir constant, m³/mol) |
| Frumkin | 3 | `kf` (Frumkin constant, m³/mol), `A` (interaction parameter) |
| Freundlich | 4 | `kfl` (Freundlich constant, mol^x m^y), `knl` (consistency index) |
| Volmer | 5 | `kv` (Volmer constant, m³/mol) |

### Common Parameters

| Parameter | Description | Default | Units |
|-----------|-------------|---------|-------|
| `D` | Diffusion coefficient | 8.8e-10 | m²/s |
| `Cb` | Bulk concentration | 1.0e-2 | mol/m³ |
| `rb` | Bubble radius | 1.12e-3 | m |
| `Tmpr` | Temperature | 298.15 | K |
| `gamma_0` | Solvent surface tension | 72.2e-3 | N/m |
| `gamma_m` | Max surface concentration | 3.91e-6 | mol/m² |
| `nn` | Value of n | 1.0 | - |
| `isotherm` | Isotherm type (1-5) | 2 | - |
| `geometry` | Geometry (0=sphere, 1=plane) | 0 | - |
| `h` | Time step (auto if None) | None | s |
| `err` | Numerical tolerance | 1e-30 | - |
| `M` | Max iterations | 100 | - |

## Examples

### Example 1: Fit Diffusion Coefficient Only

```python
from scipy.optimize import curve_fit
from ward_tordai_wrapper import create_fit_function

# Fix all parameters except D
fit_func = create_fit_function(
    fit_params=['D'],
    kl=2146.0,
    Cb=1.0e-2,
    gamma_m=3.91e-6,
    isotherm=2
)

popt, pcov = curve_fit(fit_func, t_data, st_data, p0=[8.8e-10])
D_fit = popt[0]
```

### Example 2: Fit Multiple Parameters

```python
# Fit D, kl, and Cb simultaneously
fit_func = create_fit_function(
    fit_params=['D', 'kl', 'Cb'],
    gamma_m=3.91e-6,
    isotherm=2,
    geometry=0
)

p0 = [8.8e-10, 2146.0, 1.0e-2]
bounds = ([1e-11, 100, 1e-4], [1e-8, 10000, 1e-1])

popt, pcov = curve_fit(fit_func, t_data, st_data, p0=p0, bounds=bounds)
D_fit, kl_fit, Cb_fit = popt
```

### Example 3: Using Different Isotherms

```python
# Henry isotherm
from ward_tordai_wrapper import henry_fit_function

popt, pcov = curve_fit(
    lambda t, D, kh: henry_fit_function(t, D, kh, Cb=1.0e-2),
    t_data, st_data, 
    p0=[8.8e-10, 2146.0]
)

# Frumkin isotherm  
from ward_tordai_wrapper import frumkin_fit_function

popt, pcov = curve_fit(
    lambda t, D, kf, A: frumkin_fit_function(t, D, kf, A, Cb=1.0e-2),
    t_data, st_data,
    p0=[8.8e-10, 2146.0, 0.5]
)
```

### Example 4: Planar Geometry

```python
# For planar interfaces instead of spherical bubbles
fit_func = create_fit_function(
    fit_params=['D', 'kl'],
    Cb=1.0e-2,
    gamma_m=3.91e-6,
    isotherm=2,
    geometry=1  # Planar geometry
)

popt, pcov = curve_fit(fit_func, t_data, st_data, p0=[8.8e-10, 2146.0])
```

## Running the Examples

A complete example script is provided in `example_curve_fit.py`:

```bash
python example_curve_fit.py
```

This will:
1. Generate synthetic data and fit it
2. Demonstrate different fitting approaches
3. Try to load and fit real experimental data (if available)
4. Fit multiple parameters simultaneously
5. Generate plots showing fit quality

## Tips for Successful Fitting

### 1. Good Initial Guesses
Provide reasonable initial guesses (`p0`) based on literature values or physical intuition:
- D typically ranges from 1e-11 to 1e-9 m²/s
- Isotherm constants vary widely depending on surfactant

### 2. Use Bounds
Constrain parameters to physically meaningful ranges:
```python
bounds = (
    [1e-11, 100, 1e-4],      # Lower bounds for [D, kl, Cb]
    [1e-8, 10000, 1e-1]      # Upper bounds
)
popt, pcov = curve_fit(fit_func, t_data, st_data, p0=p0, bounds=bounds)
```

### 3. Start Simple
Begin by fitting fewer parameters:
1. First fit only D with other parameters fixed
2. Then fit D and isotherm constant together
3. Finally, add more parameters if needed

### 4. Data Quality
- Ensure data covers adequate time range
- More data points generally improve fits
- Remove outliers or noisy data points

### 5. Check Fit Quality
Always evaluate the fit quality:
```python
st_fit = fit_func(t_data, *popt)
residuals = st_data - st_fit
ss_res = np.sum(residuals**2)
ss_tot = np.sum((st_data - np.mean(st_data))**2)
r_squared = 1 - (ss_res / ss_tot)
print(f"R-squared: {r_squared:.4f}")
```

### 6. Increase maxfev if Needed
For complex fits, increase the maximum function evaluations:
```python
popt, pcov = curve_fit(fit_func, t_data, st_data, p0=p0, maxfev=10000)
```

## Troubleshooting

### "Required accuracy not reached" Warning
- Increase `M` (max iterations): `M=200`
- Adjust `err` (tolerance): `err=1e-25`
- Check if parameters are physically reasonable

### Fitting Fails or Returns Unrealistic Values
- Check initial guesses `p0`
- Add bounds to constrain parameters
- Reduce number of fitted parameters
- Verify data quality and units

### Slow Fitting
- Ward-Tordai model requires numerical integration, so fitting is computational
- Use coarser time steps for initial exploration
- Consider fitting on subset of data first
- Multi-parameter fits (3+) can be slow

### Import Errors
Ensure `Ward_Tordai_Fit.py` is in the same directory or Python path:
```python
import sys
sys.path.append('/path/to/directory')
```

## Integration with Other Programs

The wrapper can be used with any optimization routine that accepts a callable:

### Using scipy.optimize.minimize
```python
from scipy.optimize import minimize

def objective(params):
    D, kl = params
    st_pred = ward_tordai_model(t_data, D=D, kl=kl, Cb=1.0e-2, isotherm=2)
    return np.sum((st_data - st_pred)**2)

result = minimize(objective, x0=[8.8e-10, 2146.0], method='Nelder-Mead')
D_fit, kl_fit = result.x
```

### Using lmfit
```python
from lmfit import Model

# Create lmfit Model
model = Model(lambda t, D, kl: ward_tordai_model(
    t, D=D, kl=kl, Cb=1.0e-2, gamma_m=3.91e-6, isotherm=2
))

params = model.make_params(D=8.8e-10, kl=2146.0)
result = model.fit(st_data, params, t=t_data)
print(result.fit_report())
```

## Citation

If you use this wrapper in your research, please cite the original Ward-Tordai implementation:

> This code is based on a Python translation of the program by X.B. Li and P. Stevenson,
> University of Newcastle, Australia. For questions: paul.stevenson@newcastle.edu.au

## License

This wrapper follows the same license and disclaimer as the original `Ward_Tordai_Fit.py` program.

**DISCLAIMER**: This program is provided as-is for educational and indicative purposes only. No guarantee is provided for accuracy or applicability. Use at your own risk.

## Contact

For issues specific to the wrapper, please check the examples and documentation above. For questions about the underlying Ward-Tordai model, contact the original authors.
