"""Functions for calculating the transmission/reflection
of uniform plane waves (UPW) through planar slabs of homogeneous dielectric
material. All of these equations are for **lossless** materials. Use Scikit-rf
for lossy materials (although it doesn't do angles).

This includes transmission/reflection through the Mylar beamsplitters, through
the vacuum window, through the Zitex RF filter, etc.

"""

import numpy as np
import scipy.constants as sc 

z0 = sc.physical_constants['characteristic impedance of vacuum'][0]

# Dielectric properties (approximate values from literature)

HDPE_N = 1.5247
MYLAR_N = 1.75
ZITEX_N = 1.22
TEFLON_N = 1.44


def lossless_slab(slab_n, slab_d, freq):
    """Lossless slab of dielectric material.

    See equation 5.4.3 in:

       http://www.ece.rutgers.edu/~orfanidi/ewa/ch05.pdf

    Args:
        slab_n (float): index of refraction
        slab_d (float): thickness, in [m]
        freq (ndarray): frequency, in [Hz]

    Returns:
        ndarray: reflected power
        ndarray: transmitted power

    """

    # Wave in freespace (fs)
    fs_c = sc.c  # velocity
    fs_z = z0  # impedance

    # Wave in dielectric slab
    slab_c = fs_c / slab_n
    slab_z = fs_z / slab_n

    # Reflection and transmission (in power, NOT E/H)
    beta = freq * 2 * np.pi / slab_c
    e_power = -1j * 2 * beta * slab_d
    r_1 = (slab_z - fs_z) / (slab_z + fs_z)
    r_2 = (fs_z - slab_z) / (slab_z + fs_z)

    reflected_power = np.abs(((r_1 + r_2 * np.exp(e_power)) /
                              (1 + r_1 * r_2 * np.exp(e_power)))) ** 2

    transmitted_power = 1 - np.abs(reflected_power)

    return reflected_power, transmitted_power


def lossless_slab_at_angle(slab_n, slab_d, freq, phii=sc.pi/4, pol='perpendicular'):
    """Lossless slab of dielectric material at an angle.

    Transmission/reflection of a uniform plane wave through a slab of
    lossless dielectric at some incident angle (phii).

    See below for equations used:

       https://www.eecis.udel.edu/~mirotzni/ELEG648/ELEG648_planewavesII.ppt

    Args:
        slab_n (float): index of refraction
        slab_d (float): thickness, in [m]
        freq (ndarray): frequency, in [Hz]
        phii (float): incident angle, in [radians]
        pol (str): type of polarization, either 'perpendicular' or 'parallel',
            default is 'perpendicular'

    Returns:
        ndarray: reflected power
        ndarray: transmitted power

    """

    assert pol == 'parallel' or pol == 'perpendicular', \
        "pol must equal either 'parallel' or 'perpendicular'"

    # Wave in freespace (fs)
    fs_c = sc.c  # velocity

    # Wave in dielectric slab
    slab_c = fs_c / slab_n

    # Reflection and transmission (in power, NOT E/H)
    beta = freq * 2 * np.pi / slab_c
    phi2 = np.arcsin(np.sin(phii) / slab_n)
    phit = np.arcsin(np.sin(phi2) * slab_n)
    if pol == 'parallel':
        r_12 = (np.cos(phii) - slab_n * np.cos(phi2)) / \
               (np.cos(phii) + slab_n * np.cos(phi2))
        r_23 = (slab_n * np.cos(phi2) - np.cos(phit)) / \
               (slab_n * np.cos(phi2) + np.cos(phit))
    else:
        r_12 = (np.cos(phii) * slab_n - np.cos(phi2)) / \
               (np.cos(phii) * slab_n + np.cos(phi2))
        r_23 = (np.cos(phi2) - slab_n * np.cos(phit)) / \
               (np.cos(phi2) + slab_n * np.cos(phit))

    e_power = -1j * 2 * beta * slab_d
    reflected_power = np.abs((r_12 + r_23 * np.exp(e_power)) /
                             (1 + r_12 * r_23 * np.exp(e_power))) ** 2

    transmitted_power = 1 - reflected_power

    return reflected_power, transmitted_power
