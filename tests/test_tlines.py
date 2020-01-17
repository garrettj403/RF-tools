"""Test the transmission line module (``rf.tlines``).

This module includes classes for rectangular and circular waveguides, and
microstrips.

"""

import numpy as np
import pytest
import scipy.constants as sc

import rftools as rf 


def test_teflon_filled_waveguide():
    """Example 3.1 on page 112 of Pozar."""

    a = 1.07 * sc.centi
    b = 0.43 * sc.centi
    er = 2.08  # approx. teflon
    tand = 0.0004  # loss tangent
    wg = rf.RectangularWaveguide(a, b, er=er, verbose=False)

    fc_te10 = wg.cutoff('TE10') / sc.giga
    fc_te20 = wg.cutoff('TE20') / sc.giga
    fc_te01 = wg.cutoff('TE01') / sc.giga
    fc_te11 = wg.cutoff('TE11') / sc.giga
    fc_te21 = wg.cutoff('TE21') / sc.giga

    # Debug
    # print(round(fc_te10, 2))
    # print(round(fc_te20, 2))
    # print(round(fc_te01, 2))
    # print(round(fc_te11, 2))
    # print(round(fc_te21, 2))

    assert abs(fc_te10 -  9.72) < 0.025
    assert abs(fc_te20 - 19.44) < 0.025
    assert abs(fc_te01 - 24.19) < 0.025
    assert abs(fc_te11 - 26.07) < 0.025
    assert abs(fc_te21 - 31.03) < 0.025


def test_teflon_filled_circular_waveguide():
    """Example 3.2 on page 123 of Pozar."""

    # Build waveguide
    a = 0.5 * sc.centi
    er = 2.08  # teflon
    wg = rf.CircularWaveguide(a, er=er, verbose=False)

    # Cutoff frequencies
    fc_te11 = wg.cutoff('TE11') / sc.giga
    fc_tm01 = wg.cutoff('TM01') / sc.giga

    # Debug
    # print("Cutoff TE11: ", round(fc_te11, 2), " GHz")
    # print("Cutoff TM01: ", round(fc_tm01, 2), " GHz")

    assert abs(fc_te11 - 12.19) < 0.025
    assert abs(fc_tm01 - 15.92) < 0.025

    # Properties at 14 GHz
    freq = 14e9
    k    = wg.k(freq)
    beta = wg.beta(freq)

    # Debug
    # print("Wavenumber:   {:.2f} m-1".format(k))
    # print("beta:         {:.2f} m-1".format(beta))

    assert abs(k    - 422.9) < 0.3
    assert abs(beta - 208.0) < 0.6


def test_microstrip_parameters():
    """Example 3.7 in Pozar."""

    # Microstrip parameters
    z0 = 50.
    d = 0.127 * sc.centi
    er = 2.20

    # Get width
    w = rf.tlines.find_microstrip_width(er, d, z0)
    assert round(w/d, 3) == 3.081

    # Get effective dielectric constant
    ee = rf.tlines._ee(er, d, w)
    assert round(ee, 2) == 1.87


def test_microstrip_impedance():
    """Example 3.8 in Pozar."""

    er = 2.55

    d = 0.1 * sc.centi

    # print(round(rf.tlines.find_microstrip_z0(er, d, 0.05 * sc.centi), 3))
    # print(round(rf.tlines.find_microstrip_z0(er, d, 0.10 * sc.centi), 3))
    # print(round(rf.tlines.find_microstrip_z0(er, d, 0.20 * sc.centi), 3))
    # print(round(rf.tlines.find_microstrip_z0(er, d, 0.40 * sc.centi), 3))
    # print(round(rf.tlines.find_microstrip_z0(er, d, 0.70 * sc.centi), 3))
    # print(round(rf.tlines.find_microstrip_z0(er, d, 1.00 * sc.centi), 3))

    assert abs(round(rf.tlines.find_microstrip_z0(er, d, 0.05 * sc.centi), 3) - 119.8) < 0.3
    assert abs(round(rf.tlines.find_microstrip_z0(er, d, 0.10 * sc.centi), 3) -  89.8) < 0.3
    assert abs(round(rf.tlines.find_microstrip_z0(er, d, 0.20 * sc.centi), 3) -  62.2) < 0.3
    assert abs(round(rf.tlines.find_microstrip_z0(er, d, 0.40 * sc.centi), 3) -  39.3) < 0.3
    assert abs(round(rf.tlines.find_microstrip_z0(er, d, 0.70 * sc.centi), 3) -  25.6) < 0.3
    assert abs(round(rf.tlines.find_microstrip_z0(er, d, 1.00 * sc.centi), 3) -  19.1) < 0.3


def test_microstrip():
    """Build microstrip and then analyze."""

    er = 3.0
    d = 0.1 * sc.milli
    z0 = 50

    w = rf.tlines.find_microstrip_width(er, d, z0)

    ms = rf.Microstrip(er, d, w, verbose=False)

    # Debug
    # print(ms.z0)

    assert round(ms.z0, 0) == 50


if __name__ == "__main__":
    test_teflon_filled_waveguide()
    test_teflon_filled_circular_waveguide()
    test_microstrip_parameters()
    test_microstrip_impedance()
    test_microstrip()
