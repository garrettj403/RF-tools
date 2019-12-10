"""Calculate the electrical properties of various transmission lines.

Ref:

   D.M. Pozar, "Microwave Engineering", 3rd edition, 2005.

"""

import cmath
import numpy as np
import scipy.constants as sc
from rftools.util import pvalf, header

Z0 = sc.physical_constants['characteristic impedance of vacuum'][0]


# Waveguides -----------------------------------------------------------------

class RectangularWaveguide:
    """Class for a rectangular waveguide.
    
    Args:
        a (float): dimension a
        b (float): dimension b

    """

    def __init__(self, a, b, **kwargs):

        verbose = kwargs.pop('verbose', True)
        comment = kwargs.pop('comment', '')

        self.a = a
        self.b = b

        if verbose:
            header("Rectangular Waveguide: {0}".format(comment))
            pvalf('a', a / sc.milli, 'mm')
            pvalf('b', b / sc.milli, 'mm')
            print("")

    def wavelength(self, frequency, m, n):
        """Calculate guided wavelength.

        Args:
            frequency (float): frequency in units Hz
            m (int): mode m
            n (int): mode n

        Returns:
            float: guided wavelength

        """

        fs_wavelength = sc.c / frequency
        k = 2 * np.pi / fs_wavelength
        kc = np.sqrt((m * np.pi / self.a)**2 + (n * np.pi / self.b)**2)
        beta = cmath.sqrt(k**2 - kc**2)
        wavelength_g = 2 * np.pi / beta

        return wavelength_g.real

    def impedance(self, frequency, m, n):
        """Calculate characteristic impedance.

        Args:
            frequency (float): frequency in units Hz
            m (int): mode m
            n (int): mode n

        Returns:
            float: characteristic impedance

        """

        fs_wavelength = sc.c / frequency
        k = 2 * np.pi / fs_wavelength
        kc = np.sqrt((m * np.pi / self.a)**2 + (n * np.pi / self.b)**2)
        beta = cmath.sqrt(k**2 - kc**2)
        z_te = k * Z0 / beta

        return z_te

    def cutoff(self, m, n):
        """Calculate cutoff frequency for mode (m,n).

        Args:
            frequency (float): frequency in units Hz
            m (int): mode m
            n (int): mode n

        Returns:
            float: cutoff frequency

        """

        kc = np.sqrt((m * np.pi / self.a)**2 + (n * np.pi / self.b)**2)
        fc = sc.c / (2 * np.pi) * kc

        return fc 
