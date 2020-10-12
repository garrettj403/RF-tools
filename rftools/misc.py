"""Miscellaneous functions."""

import numpy as np 
import scipy.constants as sc 

from numpy import pi, sqrt, arctan
from scipy.constants import mu_0, m_e, e

def skin_depth(freq, cond, mu_r=1):
    """Calculate skin depth.

    Args:
        freq (float): frequency in [Hz]
        cond (float): conductivity in [S/m]
        mu_r (float): relative permeability

    Returns:
        float: skin depth in [m]

    """

    return 1 / sqrt(pi * mu_r * mu_0 * cond * freq)


def mean_free_path(cond, fermi_speed, e_density):
    """Calculate mean free path between collisions.

    Args:
        cond (float): conductivity in [S/m]
        fermi_speed (float): Fermi speed in [m/s]
        e_density (float): electron density [m-3]

    Returns:
        float: mean free path in [m]

    """

    return cond * m_e * fermi_speed / (e_density * e ** 2)


def conductivity_4k(freq, fermi_speed, e_density, mu_r=1):
    """Calculate the effective conductivity at 4K in the regime of the 
    anomalous skin effect.

    Args:
        freq (float): frequency in [Hz]
        fermi_speed (float): Fermi speed in [m/s]
        e_density (float): electron density [m-3]
        mu_r (float): relative permeability

    Returns:
        float: effective conductivity

    """

    return (e_density ** 2 * e ** 4 / 
        (pi * mu_r * mu_0 * m_e ** 2 * fermi_speed ** 2 * freq)) ** (1 / 3)


def conductivity_rough(freq, cond, h):
    """Calculate the effective conductivity of a rough metal.

    Using the Hammerstad-Bekkadal (HB) formula.

    Args:
        freq (float): frequency in [Hz]
        cond (float): conductivity in [S/m]
        h (float): rms surface roughness in [m]

    Returns:
        float: effective conductivity

    """
    
    d = skin_depth(freq, cond)
    
    return cond * (1 + 2 / pi * arctan(1.4 * (h / d) ** 2)) ** (-2)


if __name__ == "__main__":
    """Print out the properties of copper and gold at 300K and 4K."""

    import sys

    # Get frequency (if provided)
    if len(sys.argv) > 1:
        f = float(sys.argv[1]) * 1e9
    else:
        f = 345e9

    # Copper properties (see Finger2008)
    cond = 5.9e7  # conductivity
    rrr = 106      # residual resistance ratio
    vf = 1.57e6    # fermi speed
    ne = 8.47e28   # free electron density

    print(("\n\t" + "-" * 52))
    print("\tCopper at {} GHz...".format(f / 1e9))
    print("")
    print("\tConductivity at 300 K:\t\t{:.2e}\tS/m".format(cond))
    print("\tConductivity at 4.2 K:\t\t{:.2e}\tS/m".format(cond * rrr))
    print("")
    print("\tSkin depth at 300 K:\t\t{:5.1f}\t\tnm".format(skin_depth(f, cond) / sc.nano))
    print("\tSkin depth at 4.2 K:\t\t{:5.1f}\t\tnm".format(skin_depth(f, cond * rrr) / sc.nano))
    print("")
    print("\tMean free path at 300 K:\t{:5.1f}\t\tnm".format(mean_free_path(cond, vf, ne) / sc.nano))
    print("\tMean free path at 4.2 K:\t{:5.1f}\t\tnm".format(mean_free_path(cond * rrr, vf, ne) / sc.nano))
    print("")
    print("\tLength ratio at 300 K:\t\t{:5.3f}".format(mean_free_path(cond, vf, ne) / skin_depth(f, cond)))
    print("\tLength ratio at 4.2 K:\t\t{:5.3f}".format(mean_free_path(cond * rrr, vf, ne) / skin_depth(f, cond * rrr)))
    print("")
    print("\tEffective conductivity at 4K:\t{:.2e}\tS/m".format(effective_conductivity(f, vf, ne)))
    print("\tEffective RRR:\t\t\t{:5.1f}\t\tS/m".format(effective_conductivity(f, vf, ne)/cond))
    # print("\tEffective skin depth:\t\t{:5.1f}\t\tnm".format(skin_depth(f, effective_conductivity(f, vf, ne)) / sc.nano)) 
    print("")

    # Gold properties (see Lamb1996)
    cond = 2.1e7  # conductivity
    vf = 1.40e6    # fermi speed
    ne = 5.90e28   # free electron density

    print(("\n\t" + "-" * 52))
    print("\tGold at {} GHz...".format(f / 1e9))
    print("")
    print("\tConductivity at 300 K:\t\t{:.2e}\tS/m".format(cond))
    print("")
    print("\tSkin depth at 300 K:\t\t{:5.1f}\t\tnm".format(skin_depth(f, cond) / sc.nano))
    print("")
    print("\tMean free path at 300 K:\t{:5.1f}\t\tnm".format(mean_free_path(cond, vf, ne) / sc.nano))
    print("")
    print("\tLength ratio at 300 K:\t\t{:.2f}".format(mean_free_path(cond, vf, ne) / skin_depth(f, cond)))
    print("")
    print("\tEffective conductivity at 4K:\t{:.2e}\tS/m".format(effective_conductivity(f, vf, ne)))
    print("\tEffective RRR:\t\t\t{:.1f}\t\tS/m".format(effective_conductivity(f, vf, ne)/cond))
    print("")

    # Aluminum properties (see Lamb1996)
    cond = 1.6e7  # conductivity
    rrr = 5000
    vf = 2.03e6    # fermi speed
    ne = 18.10e28   # free electron density

    print(("\n\t" + "-" * 52))
    print("\tAluminum at {} GHz...".format(f / 1e9))
    print("")
    print("\tConductivity at 300 K:\t\t{:.2e}\tS/m".format(cond))
    print("\tConductivity at 4.2 K:\t\t{:.2e}\tS/m".format(cond * rrr))
    print("")
    print("\tSkin depth at 300 K:\t\t{:5.1f}\t\tnm".format(skin_depth(f, cond) / sc.nano))
    print("\tSkin depth at 4.2 K:\t\t{:5.1f}\t\tnm".format(skin_depth(f, cond * rrr) / sc.nano))
    print("")
    print("\tMean free path at 300 K:\t{:5.1f}\t\tnm".format(mean_free_path(cond, vf, ne) / sc.nano))
    print("\tMean free path at 4.2 K:\t{:5.1f}\t\tnm".format(mean_free_path(cond * rrr, vf, ne) / sc.nano))
    print("")
    print("\tLength ratio at 300 K:\t\t{:5.3f}".format(mean_free_path(cond, vf, ne) / skin_depth(f, cond)))
    print("\tLength ratio at 4.2 K:\t\t{:5.3f}".format(mean_free_path(cond * rrr, vf, ne) / skin_depth(f, cond * rrr)))
    print("")
    print("\tEffective conductivity at 4K:\t{:.2e}\tS/m".format(effective_conductivity(f, vf, ne)))
    print("\tEffective RRR:\t\t\t{:5.1f}\t\tS/m".format(effective_conductivity(f, vf, ne)/cond))
    # print("\tEffective skin depth:\t\t{:5.1f}\t\tnm".format(skin_depth(f, effective_conductivity(f, vf, ne)) / sc.nano)) 
    print("")
