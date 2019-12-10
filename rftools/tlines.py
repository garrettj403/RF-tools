"""Calculate the electrical properties of various transmission lines.

Ref:

   D.M. Pozar, "Microwave Engineering", 3rd edition, 2005.

"""

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

    def wavelength(self, frequency, mode):
        """Calculate guided wavelength.

        Args:
            frequency (float): frequency in units Hz
            mode (str): waveguide mode

        Returns:
            float: guided wavelength

        """

        assert mode[0:2] == 'TE' or mode[0:2] == 'TM', \
            "mode must be either TE or TM"
        m, n = int(mode[2]), int(mode[3])

        fs_wavelength = sc.c / frequency
        k = 2 * np.pi / fs_wavelength
        kc = np.sqrt((m * np.pi / self.a)**2 + (n * np.pi / self.b)**2)
        beta = np.sqrt(k**2 - kc**2)
        wavelength_g = 2 * np.pi / beta

        return wavelength_g.real

    def impedance(self, frequency, mode):
        """Calculate characteristic impedance.

        Args:
            frequency (float): frequency in units Hz
            mode (str): waveguide mode

        Returns:
            float: characteristic impedance

        """

        assert mode[0:2] == 'TE' or mode[0:2] == 'TM', \
            "mode must be either TE or TM"
        m, n = int(mode[2]), int(mode[3])

        fs_wavelength = sc.c / frequency
        k = 2 * np.pi / fs_wavelength
        kc = np.sqrt((m * np.pi / self.a)**2 + (n * np.pi / self.b)**2)
        beta = np.sqrt(k**2 - kc**2)
        z_te = k * Z0 / beta

        return z_te

    def cutoff(self, mode):
        """Calculate cutoff frequency for mode (m,n).

        Args:
            frequency (float): frequency in units Hz
            mode (str): waveguide mode

        Returns:
            float: cutoff frequency

        """

        assert mode[0:2] == 'TE' or mode[0:2] == 'TM', \
            "mode must be either TE or TM"
        m, n = int(mode[2]), int(mode[3])

        kc = np.sqrt((m * np.pi / self.a)**2 + (n * np.pi / self.b)**2)
        fc = sc.c / (2 * np.pi) * kc

        return fc 


class CircularWaveguide:
    """Class for a circular waveguide.

    Args:
        a (float): inner radius a

    """

    # TODO: add impedance 

    def __init__(self, a, **kwargs):

        verbose = kwargs.pop('verbose', True)
        comment = kwargs.pop('comment', '')

        self.a = a

        if verbose:
            header("Circular Waveguide: {0}".format(comment))
            pvalf('a', a / sc.milli, 'mm')
            print("")

        # Constants
        # For TE modes (table 3.3 in pozar)
        self._pp = np.array([
            [0, 3.832, 7.016, 10.174],
            [0, 1.841, 5.331, 8.536],
            [0, 3.054, 6.706, 9.970]])
        # For TM modes (table 3.4 in pozar)
        self._p = np.array([
            [0, 2.405, 5.520, 8.654],
            [0, 3.832, 7.016, 10.174],
            [0, 5.135, 8.417, 11.620]])

    def wavelength(self, frequency, mode='TE11'):
        """Calculate guided wavelength.

        Args:
            frequency (float): frequency
            mode (str): waveguide mode, e.g., 'TE11'

        Returns:
            float: guided wavelength

        """

        assert mode[0:2] == 'TE' or mode[0:2] == 'TM', \
            "mode must be either TE or TM"

        if mode[0:2] == 'TE':
            p_temp = self._pp 
        else:
            p_temp = self._p

        p_coeff = p_temp[int(mode[2:3]), int(mode[3:4])]

        k = 2 * np.pi * frequency * np.sqrt(sc.mu_0 * sc.epsilon_0)
        kc = p_coeff / self.a
        beta = np.sqrt(k**2 - kc**2)
        wavelength_guided = 2 * np.pi / beta

        return wavelength_guided

    def cutoff(self, mode='TE11'):
        """Calculate cutoff frequency for given mode.

        Args:
            mode (str): waveguide mode

        Returns: 
            float: cutoff frequency

        """

        assert mode[0:2] == 'TE' or mode[0:2] == 'TM', \
            "mode must be either TE or TM"

        if mode[0:2] == 'TE':
            p_temp = self._pp 
        else:
            p_temp = self._p

        p_coeff = p_temp[int(mode[2:3]), int(mode[3:4])]

        return p_coeff * sc.c / (2 * np.pi * self.a)
