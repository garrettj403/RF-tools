"""Calculate the electrical properties of various transmission lines.

Ref:

   D.M. Pozar, "Microwave Engineering", 3rd edition, 2005.

"""

import numpy as np
import scipy.constants as sc
from rftools.util import pvalf, header
from rftools.conduction import surface_resistance

from numpy import sqrt, exp
from numpy import log as ln
from scipy.constants import pi 

Z0 = sc.physical_constants['characteristic impedance of vacuum'][0]


# Waveguides -----------------------------------------------------------------

class RectangularWaveguide:
    """Class for a rectangular waveguide.
    
    Args:
        a (float): dimension a (broad wall width) in [m]
        b (float): dimension b (narrow wall width) in [m]

    Keyword Args:
        er (float): relative permittivity
        ur (float): relative permeability
        verbose (bool): verbosity
        comment (str): comment to describe waveguide

    """

    def __init__(self, a, b, er=1, ur=1, **kwargs):

        verbose = kwargs.pop('verbose', True)
        comment = kwargs.pop('comment', '')

        # Waveguide dimensions
        self.a = a  # broad
        self.b = b  # narrow

        # Dielectric fill
        self.er = er 
        self.ur = ur
        self.eta = sqrt(self.ur * sc.mu_0 / self.er / sc.epsilon_0)

        # Cutoff frequency
        self.fc = sc.c / sqrt(self.er * self.ur) / 2 / self.a

        # Frequency range
        self.f1 = 1.25 * self.fc
        self.f2 = 1.89 * self.fc
        self.fmid = (self.f1 + self.f2) / 2

        if verbose:
            header("Rectangular Waveguide: {0}".format(comment))
            print("\n\tDimensions (metric):")
            if a / sc.milli > 1:
                pvalf('a', a / sc.milli, 'mm')
                pvalf('b', b / sc.milli, 'mm')
            else:
                pvalf('a', a / sc.micro, 'um')
                pvalf('b', b / sc.micro, 'um')
            print("\n\tDimensions (imperial):")
            pvalf('a', a / sc.mil, 'mil')
            pvalf('b', b / sc.mil, 'mil')
            print("")
            print("\tStandard frequency range:")
            pvalf('low', self.f1 / sc.giga, 'GHz')
            pvalf('mid', self.fmid / sc.giga, 'GHz')
            pvalf('high', self.f2 / sc.giga, 'GHz')
            print("")

    def k(self, frequency):
        """Calculate wavenumber.

        Args:
            frequency (float): frequency in [Hz]

        Returns:
            float: wavenumber in [m-1]

        """

        k = 2 * pi * frequency * sqrt(self.er * self.ur) / sc.c

        return k

    def kc(self, mode):
        """Calculate cutoff wavenumber.

        Args:
            mode (str): waveguide mode, e.g., 'TE10'

        Returns:
            float: cutoff wavenumber in [m-1]

        """

        assert mode[0:2].lower() == 'te' or mode[0:2].lower() == 'tm', \
            "mode must be either TE or TM"
        m, n = int(mode[2]), int(mode[3])

        return sqrt((m * pi / self.a)**2 + (n * pi / self.b)**2)

    def beta(self, frequency, mode):
        """Calculate propagation constant.

        Args:
            frequency (float): frequency in [Hz]
            mode (str): waveguide mode, e.g., 'TE10'

        Returns:
            float: propagation constant in [m-1]

        """

        k = self.k(frequency)
        kc = self.kc(mode)

        return sqrt(k**2 - kc**2)

    def wavelength(self, frequency, mode='TE10'):
        """Calculate guided wavelength.

        Args:
            frequency (float): frequency in [Hz]
            mode (str): waveguide mode, e.g., 'TE10'

        Returns:
            float: guided wavelength in [m]

        """

        beta = self.beta(frequency, mode)
        wavelength_g = 2 * pi / beta

        return wavelength_g.real

    def impedance(self, frequency, mode='TE10'):
        """Calculate characteristic impedance.

        Args:
            frequency (float): frequency in [Hz]
            mode (str): waveguide mode, e.g., 'TE10'

        Returns:
            float: characteristic impedance in [ohm]

        """

        k = self.k(frequency)
        beta = self.beta(frequency, mode)

        if mode[0:2].lower() == 'te':
            return k * self.eta / beta
        elif mode[0:2].lower() == 'tm':
            return beta * self.eta / k
        else:
            print("mode must be either TE or TM")
            raise

    def cutoff(self, mode='TE10'):
        """Calculate cutoff frequency for mode (m,n).

        Args:
            frequency (float): frequency, in [Hz]
            mode (str): waveguide mode, e.g., 'TE10'

        Returns:
            float: cutoff frequency, in [Hz]

        """

        kc = self.kc(mode)
        fc = sc.c / (2 * pi) / sqrt(self.er * self.ur) * kc

        return fc 

    def attenuation(self, freq, cond):
        """Calculate waveguide attenuation in Np/m.

        Using equation from:

            E. Maxwell, "Conductivity of Metallic Surfaces at Microwave 
            Frequencies," Journal of Applied Physics, vol. 18, no. 7, 
            Jul. 1947.

        Args:
            frequency (float): frequency, in [Hz]
            cond (float): conductivity at desired frequency, in [S/m]

        Returns:
            float: waveguide attenuation, in [Np/m]

        """

        # # cutoff wavelength of TE10
        # lambda_c = sc.c / self.fc
        
        # # freespace wavelength
        # lambda_0 = sc.c / freq

        # # Eqn from pg 630 of Maxwell 1947
        # att = ((1 / (2 * self.b)) *
        #        (1 / (1 - (lambda_0 / lambda_c) ** 2) ** 0.5) *
        #        ((4 * pi / (lambda_0 * sc.mu_0 * sc.c * cond)) ** 0.5) *
        #        (1 + 2 * self.b / self.a * (lambda_0 / lambda_c) ** 2))

        k = self.k(freq)
        beta = self.beta(freq, 'TE10')
        rs = surface_resistance(freq, cond)

        att = rs / self.a**3 / self.b / beta / k / self.eta * (2 * self.b * pi**2 + self.a**3 * k**2)

        return att


class CircularWaveguide:
    """Class for circular waveguides.

    Args:
        a (float): inner radius 'a'

    Keyword Args:
        er (float): relative permittivity of dielectric fill (er=1 for air)
        ur (float): relative permeability of dielectric fill (ur=1 for air)
        verbose (bool): verbosity
        comment (str): comment to describe waveguide

    """

    def __init__(self, a, er=1, ur=1, **kwargs):

        verbose = kwargs.pop('verbose', True)
        comment = kwargs.pop('comment', '')

        # Radius
        self.a = a

        # Dielectric fill
        self.er = er 
        self.ur = ur
        self.eta = sqrt(self.ur * sc.mu_0 / self.er / sc.epsilon_0)

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

        # Cutoff frequency
        self.fc = self.cutoff()

        # Frequency range (approximate)
        # same equations as rectangular waveguides
        self.f1 = 1.25 * self.fc
        self.f2 = 1.89 * self.fc
        self.fmid = (self.f1 + self.f2) / 2

        if verbose:
            header("Circular Waveguide: {0}".format(comment))
            print("\n\tDimensions:")
            pvalf('radius a', a / sc.milli, 'mm')
            print("")
            # TODO: find better rule of thumb!
            # pvalf('low freq.', self.f1 / sc.giga, 'GHz')
            # pvalf('mid freq.', self.fmid / sc.giga, 'GHz')
            # pvalf('high freq.', self.f2 / sc.giga, 'GHz')
            # print("")

    def k(self, frequency):
        """Calculate wavenumber.

        Args:
            frequency (float): frequency in [Hz]

        Returns:
            float: wavenumber in [m-1]

        """

        k = 2 * pi * frequency * sqrt(self.er * self.ur) / sc.c

        return k

    def kc(self, mode):
        """Calculate cutoff wavenumber.

        Args:
            mode (str): waveguide mode, e.g., 'TE11'

        Returns:
            float: cutoff wavenumber in [m-1]

        """

        if mode[0:2].lower() == 'te':
            p_temp = self._pp 
        elif mode[0:2].lower() == 'tm':
            p_temp = self._p
        else:
            print("Error: mode must be either TE or TM")
            raise

        p_coeff = p_temp[int(mode[2]), int(mode[3])]

        return p_coeff / self.a

    def beta(self, frequency, mode='TE11'):
        """Calculate propagation constant.

        Args:
            frequency (float): frequency in [Hz]
            mode (str): waveguide mode, e.g., 'TE11'

        Returns:
            float: propagation constant in [m-1]

        """

        k = self.k(frequency)
        kc = self.kc(mode)

        return sqrt(k**2 - kc**2)

    def wavelength(self, frequency, mode='TE11'):
        """Calculate guided wavelength.

        Args:
            frequency (float): frequency
            mode (str): waveguide mode, e.g., 'TE11'

        Returns:
            float: guided wavelength

        """

        beta = self.beta(frequency, mode)

        return 2 * pi / beta

    def cutoff(self, mode='TE11'):
        """Calculate cutoff frequency for given mode.

        Args:
            mode (str): waveguide mode

        Returns: 
            float: cutoff frequency

        """

        return self.kc(mode) * sc.c / sqrt(self.er * self.ur) / (2 * pi)

    def impedance(self, frequency, mode='TE11'):
        """Calculate characteristic impedance.

        Args:
            frequency (float): frequency in [Hz]
            mode (str): waveguide mode, e.g., 'TE11'

        Returns:
            float: characteristic impedance in [ohms]

        """

        k = self.k(frequency)
        beta = self.beta(frequency, mode)

        if mode[0:2].lower() == 'te':
            return self.eta * k / beta
        elif mode[0:2].lower() == 'tm':
            return self.eta * beta / k


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

#         self.z0 = d / w * sqrt(sc.mu_0 / sc.epsilon_0 / er)
#         # self.ee = _ee(er, d, w)
#         # self.vp = sc.c / sqrt(self.ee)
#         # self.z0 = find_microstrip_z0(er, d, w)

#         if verbose:
#             parallel_plate_z0 = d / w * Z0 / sqrt(self.ee)
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

#     #     return sc.c / frequency / sqrt(self.ee)


# Microstrip -----------------------------------------------------------------

class Microstrip:
    """Class for a microstrip.

    Args:
        er (float): relative permittivity of the dielectric
        d (float): thickness of the dielectric in [m]
        w (float): width of the microstrip in [m]

    """

    def __init__(self, er, d, w, **kwargs):

        verbose = kwargs.pop('verbose', True)
        comment = kwargs.pop('comment', '')

        self.er = er
        self.d = d
        self.w = w 

        self.ee = _ee(er, d, w)
        self.vp = sc.c / sqrt(self.ee)
        self.z0 = find_microstrip_z0(er, d, w)

        if verbose:
            parallel_plate_z0 = d / w * Z0 / sqrt(self.ee)
            fringing_factor = self.z0 / parallel_plate_z0
            
            header("\nMicrostrip {0}".format(comment))
            pvalf('w', w / sc.micro, 'um')
            pvalf('d', d / sc.micro, 'um')
            pvalf('w / d', w / d)
            print("")

    def impedance(self):
        """Calculate characteristic impedance.

        Returns:
            float: characteristic impedance in [ohm]

        """

        return self.z0

    def wavelength(self, frequency):
        """Calculate wavelength.

        Args:
            frequency (float): frequency in [Hz]

        Returns:
            float: wavelength in [m]

        """

        return sc.c / frequency / sqrt(self.ee)

    def dielectric_loss(self, skin_depth):
        """Calculate the attenuation due to dielectric loss in Np/m.

        Equation 3.198 in Pozar.

        Args:
            skin_depth (float): skin depth

        Returns:
            attenuation due to dielectric loss in Np/m

        """

        alpha_d = (k0 * self.er * (self.ee - 1) * tan(skin_depth) / 
                   (2 * sqrt(self.ee) * (self.er - 1)))

        return alpha_d

    def conductor_loss(self, freq, cond):
        """Calculate the attenuation due to conductor loss in Np/m.

        Equation 3.199 in Pozar.

        Args:
            freq (float): frequency
            cond (float): conductivity of conductor

        Returns:
            attenuation due to conductor loss

        """

        # Surface resistance of conductor
        rs = sqrt(2 * pi * freq * sc.mu_0 / 2 / cond)

        return rs / self.z0 / self.w

    def attenuation(self, rs, skin_depth):
        """Calculate attenuation in Np/m (including dielectric and conductor 
        loss).

        Args:
            rs (float): surface resistance of the conductor
            skin_depth (float): skin depth of the dielectric

        Returns:
            attenuation in Np/m

        """

        return self.conductor_loss(rs) + self.dielectric_loss(skin_depth)


def find_microstrip_width(er, d, z0=50.):
    """Calculate microstrip width required for a given characteristic 
    impedance.

    Args:
        er (float): relative permittivity
        d (float): thickness of the dielectric in [m]
        z0 (float): desired characterisitic impedance in [ohm]

    Returns:
        float: width of microstrip in [m]

    """

    a = z0 / 60 * sqrt((er + 1) / 2) + (er - 1) / (er + 1) * (0.23 + 0.11 / er)
    b = 60 * pi ** 2 / (z0 * sqrt(er))

    # Eqn. 3.197 in Pozar
    wd1 = 8 * exp(a) / (exp(2 * a) - 2)
    wd2 = 2 / pi * (b - 1 - ln(2 * b - 1) + (er - 1) / (2 * er) * (ln(b - 1) + 0.39 - 0.61 / er))

    if isinstance(wd1, float):
        if a > 1.52:
            return wd1 * d
        elif wd2 >= 2:
            return wd2 * d
    elif isinstance(wd1, np.ndarray):
        wd = np.empty_like(wd1)
        mask = a > 1.52
        wd[mask] = wd1[mask]
        wd[~mask] = wd2[~mask]
        return wd * d 


def find_microstrip_z0(er, d, w):
    """Calculate characterisitic impedance of microstrip.

    Closed-form expressions from Wheeler and Schneider:
    - Eqn. 3.196 in Pozar
    - Eqn. 2.116 in Gupta

    Args:
        er (float): relative permittivity of the dielectric
        d (float): thickness of the dielectric in [m]
        w (float): width of the microstrip in [m]

    Returns:
        float: characteristic impedance in [ohm]

    """

    # Effective dielectric constant
    # Note: if w / d is an array, ee will be an array as well
    ee = _ee(er, d, w)

    # Width to thickness ratio
    wd = w / d

    if isinstance(wd, float):
        if wd <= 1:
            return 60 / sqrt(ee) * ln(8 / wd + wd / 4)
        else:
            return 120 * pi / (sqrt(ee) * (wd + 1.393 + 0.667 * ln(wd + 1.444)))
    elif isinstance(wd, np.ndarray):
        z0 = np.empty_like(wd)
        for i, wdi in np.ndenumerate(wd):
            if wdi <= 1:
                z0[i] = 60 / sqrt(ee[i]) * ln(8 / wdi + wdi / 4)
            else:
                z0[i] = 120 * pi / (sqrt(ee[i]) * (wdi + 1.393 + 0.667 * ln(wdi + 1.444)))
        return z0


def _ee(er, d, w):
    """Effective dielectric permittivity for a microstrip.

    - Eqn. 3.195 in Pozar
    - Eqn. 2.117 in Gupta

    Args:
        er (float): relative permittivity of the dielectric material
        d (float): thickness of the dielectric in [m]
        w (float): width of the microstrip in [m]

    Returns:
        float: effective permittivity

    """

    # Width to thickness ratio
    wd = w / d

    if isinstance(wd, float):
        f = 1 / sqrt(1 + 12 / wd)
        if wd <= 1:
            f += 0.04 * (1 - wd) ** 2
        return (er + 1) / 2 + (er - 1) / 2 * f
    elif isinstance(wd, np.ndarray):
        mask = wd <= 1
        f = 1 / sqrt(1 + 12 / wd)
        f[mask] += 0.04 * (1 - wd[mask]) ** 2
        return (er + 1) / 2 + (er - 1) / 2 * f


# TODO: Move to util
def db10(value):
    """Calculate dB from power value.

    Args:
        value (float): value to convert (units of power)

    Returns:
        value in dB

    """

    return 10 * np.log10(np.abs(value))


def db20(value):
    """Calculate dB from voltage value.

    Args:
        value (float): value to convert (units of voltage)

    Returns:
        value in dB

    """

    return 20 * np.log10(np.abs(value))


# TODO: implement
# def np2db(value_in_np):
#     return 

# def db2np(value_in_db):
#     return


# ----------------------------------------------------------------------------

if __name__ == "__main__":

    import matplotlib.pyplot as plt 

    ### Microstrip ###

    # Frequency
    f = 320 * sc.giga

    # Relative permittivity
    er = 4.4 

    # Dielectric thickness
    d = 250 * sc.nano

    # Width (sweep)
    w = np.linspace(0.001, 20, 1001) * sc.micro

    # Calculate characteristic impedance for each width
    z0 = find_microstrip_z0(er, d, w)

    # Recover width from characteristic impedance
    w2 = find_microstrip_width(er, d, z0)

    # Effective dielectric constant
    ee = _ee(er, d, w)

    # Wavelength
    wavel = sc.c / f / sqrt(ee)

    # Print to terminal
    print("")
    header("Microstrip")
    print("")
    pvalf("Frequency", f/sc.giga, "GHz")
    print("")
    pvalf("Rel. permittivity", er, "")
    pvalf("Diel. thickness", d / sc.nano, "nm")
    print("")

    # Plot characteristic impedance vs. width
    plt.figure()
    plt.plot(w / sc.micro, z0)
    plt.xlabel("Width (um)")
    plt.ylabel("Characteristic impedance (ohms)")
    plt.xlim([0, w.max() / sc.micro])
    plt.ylim(ymin=0)

    # Effective permittivity vs. width
    plt.figure()
    plt.semilogx(w / sc.micro, ee)
    plt.axhline(er, c='k', ls='--')
    plt.axvline(d / sc.micro, c='k', lw=0.5)
    plt.xlabel("Width (um)")
    plt.ylabel("Relative Permittivity")
    plt.xlim([w.min() / sc.micro, w.max() / sc.micro])

    # # Error on recovered width (from going back and forth)
    # plt.figure()
    # plt.plot(w / sc.micro, (w - w2) / w * 100)
    # plt.xlabel("Width (um)")
    # plt.ylabel("Width Error (%)")
    # plt.xlim([0, w.max() / sc.micro])

    # Wavelength vs. width
    plt.figure()
    plt.plot(w / sc.micro, wavel / sc.micro)
    plt.xlabel("Width (um)")
    plt.ylabel("Wavelength (um)")
    plt.xlim([0, w.max() / sc.micro])

    plt.show()
