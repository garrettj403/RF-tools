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

    Keyword Args:
        er (float): relative permittivity
        ur (float): relative permeability
        verbose (bool): verbosity
        comment (str): comment to describe waveguide

    """

    def __init__(self, a, b, er=1, ur=1, **kwargs):

        verbose = kwargs.pop('verbose', True)
        comment = kwargs.pop('comment', '')

        self.a = a
        self.b = b

        self.er = er 
        self.ur = ur

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

        assert mode[0:2].lower() == 'te' or mode[0:2].lower() == 'tm', \
            "mode must be either TE or TM"
        m, n = int(mode[2]), int(mode[3])

        fs_wavelength = sc.c / np.sqrt(self.er * self.ur) / frequency
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

        assert mode[0:2].lower() == 'te' or mode[0:2].lower() == 'tm', \
            "mode must be either TE or TM"
        m, n = int(mode[2]), int(mode[3])

        fs_wavelength = sc.c / np.sqrt(self.er * self.ur) / frequency
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

        assert mode[0:2].lower() == 'te' or mode[0:2].lower() == 'tm', \
            "mode must be either TE or TM"
        m, n = int(mode[2]), int(mode[3])

        kc = np.sqrt((m * np.pi / self.a)**2 + (n * np.pi / self.b)**2)
        fc = sc.c / np.sqrt(self.er * self.ur) / (2 * np.pi) * kc

        return fc 


class CircularWaveguide:
    """Class for a circular waveguide.

    Args:
        a (float): inner radius a

    Keyword Args:
        er (float): relative permittivity
        ur (float): relative permeability
        verbose (bool): verbosity
        comment (str): comment to describe waveguide

    """

    # TODO: add impedance 

    def __init__(self, a, er=1, ur=1, **kwargs):

        verbose = kwargs.pop('verbose', True)
        comment = kwargs.pop('comment', '')

        self.a = a

        self.er = er 
        self.ur = ur

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

        if mode[0:2].lower() == 'te':
            p_temp = self._pp 
        elif mode[0:2].lower() == 'tm':
            p_temp = self._p
        else:
            print("Error: mode must be either TE or TM")
            raise

        p_coeff = p_temp[int(mode[2:3]), int(mode[3:4])]

        k = 2 * np.pi * sc.c * frequency * np.sqrt(self.er * self.ur)
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

        if mode[0:2].lower() == 'te':
            p_temp = self._pp 
        elif mode[0:2].lower() == 'tm':
            p_temp = self._p
        else:
            print("Error: mode must be either TE or TM")
            raise

        p_coeff = p_temp[int(mode[2:3]), int(mode[3:4])]

        return p_coeff * sc.c / np.sqrt(self.er * self.ur) / (2 * np.pi * self.a)


# class ParallelPlateWaveguide:
#     """
#     Class for a parallel plate waveguide.
#     """

#     def __init__(self, er, d, w, **kwargs):
#         """
#         Initialize microstrip.

#         :param er: relative permittivity of the dielectric
#         :param d: thickness of dielectric slab
#         :param w: width of microstrip
#         """

#         verbose = kwargs.pop('verbose', True)
#         comment = kwargs.pop('comment', '')

#         self.er = er
#         self.d = d
#         self.w = w 

#         self.z0 = d / w * np.sqrt(sc.mu_0 / sc.epsilon_0 / er)
#         # self.ee = _ee(er, d, w)
#         # self.vp = sc.c / np.sqrt(self.ee)
#         # self.z0 = find_microstrip_z0(er, d, w)

#         if verbose:
#             parallel_plate_z0 = d / w * Z0 / np.sqrt(self.ee)
#             fringing_factor = self.z0 / parallel_plate_z0
            
#             header("Parallel plate {0}".format(comment))
#             pvalf('w', w / sc.micro, 'um')
#             pvalf('d', d / sc.micro, 'um')
#             # pvalf('w / d', w / d)
#             # pvalf('e_r', er)
#             # pvalf('e_eff', self.ee)
#             # pvalf('v_p', self.vp / sc.c, 'c')
#             pvalf('Z0', self.z0, 'ohms')
#             print("")
#             # pvalf('par. plate Z0', parallel_plate_z0, 'ohms')
#             # pvalf('finging fac.', (fringing_factor - 1.) * 100.)
#             # print ""

#     def impedance(self):
#         """
#         Find characteristic impedance of microstrip.

#         :return: characteristic impedance
#         """

#         return self.z0

#     # def wavelength(self, frequency):
#     #     """
#     #     Find wavelength of microstrip.

#     #     :param float frequency: frequency (in Hz)
#     #     :return: wavelength (in m)
#     #     """

#     #     return sc.c / frequency / np.sqrt(self.ee)


# Microstrip -----------------------------------------------------------------

class Microstrip:
    """Class for a microstrip.

    Args:
        er (float): relative permittivity of the dielectric
        d (float): thickness of the dielectric
        w (float): width of the microstrip

    """

    def __init__(self, er, d, w, **kwargs):

        verbose = kwargs.pop('verbose', True)
        comment = kwargs.pop('comment', '')

        self.er = er
        self.d = d
        self.w = w 

        self.ee = _ee(er, d, w)
        self.vp = sc.c / np.sqrt(self.ee)
        self.z0 = find_microstrip_z0(er, d, w)

        if verbose:
            parallel_plate_z0 = d / w * Z0 / np.sqrt(self.ee)
            fringing_factor = self.z0 / parallel_plate_z0
            
            header("\nMicrostrip {0}".format(comment))
            pvalf('w', w / sc.micro, 'um')
            pvalf('d', d / sc.micro, 'um')
            pvalf('w / d', w / d)
            print("")

    def impedance(self):
        """Calculate characteristic impedance.

        Returns:
            float: characteristic impedance

        """

        return self.z0

    def wavelength(self, frequency):
        """Calculate wavelength.

        Args:
            frequency (float): frequency

        Returns:
            float: wavelength

        """

        return sc.c / frequency / np.sqrt(self.ee)


def find_microstrip_width(er, d, z0=50.):
    """Calculate microstrip width required for a given characteristic 
    impedance.

    Args:
        er (float): relative permittivity
        d (float): thickness of the dielectric
        z0 (float): desired characterisitic impedance

    Returns:
        float: width of microstrip

    """

    a = z0 / 60 * np.sqrt((er + 1) / 2) + (er - 1) / (er + 1) * (0.23 + 0.11 / er)
    b = 377 * sc.pi / (2 * z0 * np.sqrt(er))

    wd1 = 8 * np.exp(a) / (np.exp(2 * a) - 2)
    wd2 = 2 / sc.pi * (b - 1 - np.log(2 * b - 1) + (er - 1) / (2 * er) * (np.log(b - 1) + 0.39 - 0.61 / er))

    if wd1 < 2:
        return wd1 * d
    elif wd2 > 2:
        return wd2 * d
    else:
        raise ValueError


def find_microstrip_z0(er, d, w):
    """Calculate characterisitic impedance.

    Args:
        er (float): relative permittivity of the dielectric
        d (float): thickness of the dielectric
        w (float): width of the microstrip

    Returns:
        float: characteristic impedance

    """

    ee = _ee(er, d, w)

    if w / d <= 1:
        return 60 / np.sqrt(ee) * np.log(8 * d / w + w / 4 / d)
    else:
        return 120 * sc.pi / (np.sqrt(ee) * \
               (w / d + 1.393 + 0.667 * np.log(w / d + 1.444)))


def _ee(er, d, w):
    """Effective dielectric permittivity."""

    return (er + 1) / 2 + (er - 1) / 2 / np.sqrt(1 + 12 * d / w)
