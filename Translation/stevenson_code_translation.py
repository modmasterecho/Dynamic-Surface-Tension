#!/usr/bin/env python3
"""
Python translation of the "Stevenson Code" C++ program included in Literature/Stevenson Code.doc

Original authors: Xueliang Bruce Li and Paul Stevenson (University of Newcastle, Australia)
Source summary: Diffusion-controlled adsorption rate and Dynamic Surface Tension (DST) calculator

Notes on translation:
- Preserves structure and behavior of the original interactive console program
- Uses double precision floats by default (Python float)
- Writes output to a TSV file named Simulation.tsv (tab-separated), instead of .xls
- The numerical routine mirrors the original simple fixed-step + bisection solver
- Minor guardrails added (input validation, safe sqrt, divide-by-zero guards)

This translation is intended for educational and reproducibility purposes only, matching
the legacy program behavior as closely as practical in Python.
"""

from __future__ import annotations
import numpy as np
import math
import sys
from dataclasses import dataclass
from typing import Callable, Tuple


PI = np.pi
R = 8.314  # Gas constant (J/mol/K)


@dataclass
class ModelParams:
    # Global parameters (initialized with defaults from the C++ code)
    cb: float = 0.0           # bulk concentration (mol/m^3)
    dif: float = 0.0          # molecular diffusivity (m^2/s) [C++ comment says m/s^2 but WT uses m^2/s]
    rb: float = 0.0           # radius of curvature for spherical geometry (m)
    Tmpr: float = 298.0       # temperature (K)
    st0: float = 0.072        # solvent surface tension (N/m)

    # Isotherm constants (varies by choice)
    gamma_m: float = 0.0
    kh: float = 0.0
    kl: float = 0.0
    kf: float = 0.0
    A: float = 0.0
    nn: float = 1.0
    kfl: float = 0.0
    knl: float = 0.0
    kv: float = 0.0

    # Time stepping controls
    T: float = 3000.0
    h: float = 1.0
    err: float = 1e-30
    M: int = 100

    # Selections
    isotherm: int = 1   # 1..5
    geometry: int = 0   # 0 spherical, 1 planar


# Isotherm helper functions (mirroring the inline functions in C++)

def hr(gamma: float, p: ModelParams) -> float:
    # Henry's law: gamma/kh
    if p.kh == 0.0:
        return float("inf")
    return gamma / p.kh


def lm(gamma: float, p: ModelParams) -> float:
    # Langmuir: 1/kl * gamma/(gamma_m - gamma)
    denom = (p.gamma_m - gamma)
    if p.kl == 0.0 or denom == 0.0:
        return float("inf")
    return (1.0 / p.kl) * (gamma / denom)


def fk(gamma: float, p: ModelParams) -> float:
    # Frumkin: (1/kf)*gamma/(gamma_m - gamma) * exp(-A*(gamma/gamma_m))
    denom = (p.gamma_m - gamma)
    if p.kf == 0.0 or denom == 0.0 or p.gamma_m == 0.0:
        return float("inf")
    return (1.0 / p.kf) * (gamma / denom) * math.exp(-p.A * (gamma / p.gamma_m))


def fl_(gamma: float, p: ModelParams) -> float:
    # Freundlich: (gamma/kfl)^knl
    if p.kfl == 0.0:
        return float("inf")
    if gamma < 0:
        # avoid complex for negative gamma with fractional powers
        return float("nan")
    return (gamma / p.kfl) ** p.knl


def vm(gamma: float, p: ModelParams) -> float:
    # Volmer: kv*(gamma/(gamma_m-gamma))*exp(gamma/(gamma_m-gamma))
    denom = (p.gamma_m - gamma)
    if p.kv == 0.0 or denom == 0.0:
        return float("inf")
    frac = gamma / denom
    return p.kv * frac * math.exp(frac)


# Geometry term g(geometry, t)

def g_term(geometry: int, t: float, p: ModelParams) -> float:
    if t < 0:
        return float("nan")
    if geometry == 0:  # spherical (convex)
        # (2*sqrt(t*dif/pi)+ dif/rb*t)*cb
        rt = max(t, 0.0)
        core = 2.0 * math.sqrt(rt * p.dif / PI) + (p.dif / p.rb) * rt if p.rb != 0.0 else float("inf")
        return core * p.cb
    elif geometry == 1:  # planar
        # sqrt(dif/pi)*2*cb*sqrt(t)
        return math.sqrt(p.dif / PI) * 2.0 * p.cb * math.sqrt(max(t, 0.0))
    else:
        raise ValueError("geometry must be 0 (spherical) or 1 (planar)")


# Surface tension change stn(isotherm, nn, gamma)

def stn(isotherm: int, nn: float, gamma: float, p: ModelParams) -> float:
    if isotherm == 1:  # Henry
        return -nn * R * p.Tmpr * gamma
    elif isotherm == 2:  # Langmuir
        if p.gamma_m == 0.0:
            return float("nan")
        return nn * R * p.Tmpr * p.gamma_m * math.log(1.0 - gamma / p.gamma_m)
    elif isotherm == 3:  # Frumkin
        if p.gamma_m == 0.0:
            return float("nan")
        x = gamma / p.gamma_m
        return nn * R * p.Tmpr * p.gamma_m * math.log(1.0 - x) + R * nn * p.Tmpr * p.A / 2.0 * p.gamma_m * x * x
    elif isotherm == 4:  # Freundlich
        return nn * p.knl * R * p.Tmpr * gamma
    elif isotherm == 5:  # Volmer
        if p.gamma_m == 0.0:
            return float("nan")
        return nn * p.gamma_m * p.gamma_m / (p.gamma_m - gamma) * R * p.Tmpr
    else:
        raise ValueError("isotherm must be 1..5")


# Kernel K(t, tau, gamma, isotherm, geometry)

def K(t: float, tau: float, gamma: float, p: ModelParams) -> float:
    if t <= tau:
        return 0.0

    sqrt_term = math.sqrt(p.dif / PI) / math.sqrt(max(t - tau, 0.0))
    add_spherical = p.dif / p.rb if (p.geometry == 0 and p.rb != 0.0) else 0.0

    if p.geometry == 0:  # spherical
        factor = -(sqrt_term + add_spherical)
    elif p.geometry == 1:  # planar
        factor = -sqrt_term
    else:
        raise ValueError("geometry must be 0 or 1")

    if p.isotherm == 1:
        return factor * hr(gamma, p)
    elif p.isotherm == 2:
        return factor * lm(gamma, p)
    elif p.isotherm == 3:
        return factor * fk(gamma, p)
    elif p.isotherm == 4:
        return factor * fl_(gamma, p)
    elif p.isotherm == 5:
        return factor * vm(gamma, p)
    else:
        raise ValueError("isotherm must be 1..5")


# Bisection root finder equivalent to rtbis in the C++ code

"""Alternativ mit scipy.optimize.bisect?"""

def rtbis(tn: float, x1: float, x2: float, p: ModelParams) -> Tuple[float, bool]:
    def f(x: float) -> float:
        return x - p.h / 2.0 * K(tn, 0.9999 * tn, x, p) - x2

    dx = 0.0
    f1 = f(x1)
    # Mimic original logic choosing the side based on sign of f(x1)
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


def collect_parameters(p: ModelParams) -> None:
    # Advanced mode
    option = input("Advanced mode to change tolerance and max iterations? (Y/N): ").strip().lower()
    if option == 'y':
        p.err = safe_float("Tolerance (default 1e-30): ")
        p.M = safe_int("Maximum iterations (default 100): ")

    print("\n-------------------------------------------------------------------------------")
    print("Enter model parameters. Use suggested units and scientific notation (e.g., 1e-6).")

    p.Tmpr = safe_float("Thermodynamic temperature, T (K): ")
    p.st0 = safe_float("Solvent surface tension, sigma0 (N/m): ")
    p.cb = safe_float("Bulk concentration, cb (mol/m^3): ")
    p.dif = safe_float("Molecular diffusivity, D (m^2/s): ")
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
        # spherical needs rb for all; others also need their constants
        p.rb = safe_float("Radius of curvature, rb (m): ")

    # per-isotherm constants
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


def run_simulation(p: ModelParams) -> Tuple[list[float], list[float], list[float], bool]:
    # Arrays sized by N = ceil(T/h)
    if p.h <= 0.0:
        raise ValueError("Time step h must be positive")

    N = int(math.ceil(p.T / p.h))
    t = [0.0] * (N + 1)
    gamma = [0.0] * (N + 1)
    st = [p.st0] * (N + 1)

    for i in range(N + 1):
        t[i] = p.h * i

    acc = True  # Will be set False if any bisection fails to meet accuracy

    for n in range(1, N + 1):
        # Sum over previous steps
        ssum = 0.0
        for j in range(1, n):
            ssum += K(t[n], t[j], gamma[j], p)

        x = g_term(p.geometry, t[n], p) + p.h * ssum
        root, ok = rtbis(t[n], gamma[n - 1], x, p)
        gamma[n] = root
        st[n] = p.st0 + stn(p.isotherm, p.nn, gamma[n], p)
        acc = acc and ok

    return t, gamma, st, acc


def write_results(path: str, t: list[float], gamma: list[float], st: list[float]) -> None:
    with open(path, 'w', encoding='utf-8') as f:
        # Use TSV with header for clarity
        f.write("time_s\tsurface_excess_mol_m2\tdynamic_surface_tension_N_m\n")
        for ti, gi, si in zip(t, gamma, st):
            f.write(f"{ti}\t{gi}\t{si}\n")


def main(argv: list[str]) -> int:
    banner()
    explain_models()

    while True:
        p = ModelParams()
        collect_parameters(p)

        t, gamma, st, acc = run_simulation(p)
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
        again = input("Continue? (Y/N): ").strip().lower()
        if again != 'y':
            break

    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
