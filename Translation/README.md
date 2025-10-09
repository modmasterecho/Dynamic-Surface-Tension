# Stevenson Code — Python Translation

This folder contains a faithful Python translation of the C++ program found in `Literature/Stevenson Code.doc`.

What it does:
- Computes diffusion-controlled adsorption and dynamic surface tension (DST)
- Supports spherical (convex) and planar geometries
- Includes five isotherms: Henry, Langmuir, Frumkin, Freundlich, Volmer
- Writes results to a tab-separated file `Simulation.tsv`

## How to run

- Requires Python 3.9+ (standard library only)
- Run the script and follow the interactive prompts:

```powershell
python .\Translation\stevenson_code_translation.py
```

Results will be saved to `Simulation.tsv` in the current working directory (overwritten each run).

## Notes
- The original C++ comments list diffusivity units as `m/s^2`; here we assume standard `m^2/s` for diffusion.
- The numerical method mirrors the original simple bisection routine and fixed time step.
- If spherical geometry is chosen, you will be asked for the radius `rb` as in the original.

## Provenance / Disclaimer
This translation is intended for educational and reproducibility purposes, based on the code by X.B. Li and P. Stevenson (University of Newcastle, Australia). If used for publication, please cite their contribution. No guarantees are made about accuracy or applicability; use at your own risk.