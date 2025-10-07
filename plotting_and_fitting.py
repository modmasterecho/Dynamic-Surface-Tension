"""
Plotting and Model Fitting Module

This module provides functions for plotting experimental data and fitting
the Ward-Tordai model to measured surface tension data.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import linregress
import pandas as pd


def plot_surface_tension_vs_time(time_data, surface_tension_data, title="Dynamic Surface Tension"):
    """
    Plot surface tension as a function of time.
    
    Parameters:
    -----------
    time_data : array-like
        Time values (typically in seconds)
    surface_tension_data : array-like
        Surface tension values (typically in mN/m)
    title : str
        Plot title
    """
    plt.figure(figsize=(10, 6))
    plt.plot(time_data, surface_tension_data, 'o-', label='Experimental Data')
    plt.xlabel('Time (s)')
    plt.ylabel('Surface Tension (mN/m)')
    plt.title(title)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.show()


def fit_model(time_data, surface_tension_data, model_function, initial_params):
    """
    Fit a model function to experimental data.
    
    Parameters:
    -----------
    time_data : array-like
        Time values
    surface_tension_data : array-like
        Surface tension values
    model_function : callable
        Model function to fit
    initial_params : array-like
        Initial parameter guesses
        
    Returns:
    --------
    popt : array
        Optimal parameters
    pcov : array
        Covariance matrix
    """
    popt, pcov = curve_fit(model_function, time_data, surface_tension_data, p0=initial_params)
    return popt, pcov


def plot_data_with_fit(time_data, surface_tension_data, fitted_data, title="Data with Model Fit"):
    """
    Plot experimental data along with fitted model.
    
    Parameters:
    -----------
    time_data : array-like
        Time values
    surface_tension_data : array-like
        Experimental surface tension values
    fitted_data : array-like
        Fitted surface tension values
    title : str
        Plot title
    """
    plt.figure(figsize=(10, 6))
    plt.plot(time_data, surface_tension_data, 'o', label='Experimental Data', markersize=8)
    plt.plot(time_data, fitted_data, '-', label='Model Fit', linewidth=2)
    plt.xlabel('Time (s)')
    plt.ylabel('Surface Tension (mN/m)')
    plt.title(title)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.show()


def residual_plot(time_data, surface_tension_data, fitted_data):
    """
    Plot residuals between experimental data and model fit.
    
    Parameters:
    -----------
    time_data : array-like
        Time values
    surface_tension_data : array-like
        Experimental surface tension values
    fitted_data : array-like
        Fitted surface tension values
    """
    residuals = surface_tension_data - fitted_data
    
    plt.figure(figsize=(10, 6))
    plt.plot(time_data, residuals, 'o-')
    plt.axhline(y=0, color='r', linestyle='--', alpha=0.5)
    plt.xlabel('Time (s)')
    plt.ylabel('Residuals (mN/m)')
    plt.title('Residual Plot')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()
