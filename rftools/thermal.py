"""Thermal continuum models."""

import numpy as np 
import scipy.constants as sc 

WIEN = sc.physical_constants['Wien frequency displacement law constant'][0]


def blackbody(variable, temp, model='planck', variable_type='frequency'):
    """Black body spectral radiance.

    Args:
        variable: frequency in [Hz] or wavelength in [m]
        temp: temperature in [K]
        model: black body model, either 'planck' or 'rj'
        variable_type: variable type, either 'frequency' or 'wavelength'

    Returns:
        spectral radiance

    """

    # Planck vs frequency
    if model.lower() == 'planck' and variable_type.lower() == 'frequency':
        term1 = 2 * sc.h * variable ** 3 / sc.c ** 2
        term2 = np.exp(sc.h * variable / (sc.k * temp)) - 1
        return term1 / term2

    # Planck vs wavelength
    if model.lower() == 'planck' and variable_type.lower() == 'wavelength':
        term1 = 2 * sc.h * sc.c ** 2 / variable ** 5
        term2 = np.exp(sc.h * sc.c / (variable * sc.k * temp)) - 1
        return term1 / term2

    # Rayleigh-Jeans vs frequency
    if model.lower() == 'rj' and variable_type.lower() == 'frequency':
        return 2 * variable ** 2 * sc.k * temp / sc.c ** 2

    # Rayleigh-Jeans vs wavelength
    if model.lower() == 'rj' and variable_type.lower() == 'wavelength':
        return 2 * sc.c * sc.k * temp / variable ** 4


def wien(temp):
    """Wien's displacement law.

    Args:
        temp: temperature in [K]

    Returns:
        frequency at which the blackbody curve peaks 

    """

    return temp * WIEN


# def temp_cw(freq, tphys):

#     freq = float(freq)
#     tphys = float(tphys)

#     return sc.h * freq / 2 / sc.k / np.tanh(sc.h * freq / 2 / sc.k / tphys)

