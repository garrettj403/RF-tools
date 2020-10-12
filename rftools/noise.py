"""Functions to calculate noise properties."""

import scipy.constants as sc
import numpy as np 


def calculate_tn(yfac, thot=293, tcold=78):
    """Calculate the noise temperature from the Y-factor.

    Args:
        yfac (ndarray): Y-factor

    Keyword Args:
        thot (float): hot load RJ temperature in [K], default is 293 K
        tcold (float): cold load RJ temperature in [K], default is 78 K

    Returns:
        ndarray: noise temperature

    """

    return (thot - yfac * tcold) / (yfac - 1)


def temp_cw(freq, tphys):
    """Callen-Welton temperature. 

    Uses Planck distribution with half photon.

    Args:
        freq (ndarray): frequency in [Hz]
        tphys (float): physical temperature in [K]

    Returns:
        ndarray: equivalent temperature

    """

    return sc.h * freq / 2 / sc.k / np.tanh(sc.h * freq / 2 / sc.k / tphys)


def temp_rj(freq, tphys):
    """Rayleigh-Jeans (RJ) equivalent temperature. 

    Args:
        freq (ndarray): frequency in [Hz]
        tphys (float): physical temperature in [K]

    Returns:
        ndarray: RJ equivalent temperature

    """

    return (sc.h * freq / sc.k) / (np.exp(sc.h * freq / (sc.k * tphys)) - 1)
