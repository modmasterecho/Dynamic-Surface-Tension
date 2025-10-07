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
