"""Test the transmission line module."""

import numpy as np
import pytest
import scipy.constants as sc

import rftools as rf 


def test_teflon_filled_waveguide():
    """Example 3.1 on page 112 of Pozar."""

    a = 1.07 * sc.centi
    b = 0.43 * sc.centi
    er = 2.08  # teflon
    wg = rf.RectangularWaveguide(a, b, er=er, verbose=False)

    fc_te10 = wg.cutoff('TE10') / sc.giga
    fc_te20 = wg.cutoff('TE20') / sc.giga
    fc_te01 = wg.cutoff('TE01') / sc.giga
    fc_te11 = wg.cutoff('TE11') / sc.giga
    fc_te21 = wg.cutoff('TE21') / sc.giga

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

    a = 0.5 * sc.centi
    er = 2.08  # teflon
    wg = rf.CircularWaveguide(a, er=er, verbose=False)

    fc_te11 = wg.cutoff('TE11') / sc.giga
    fc_tm01 = wg.cutoff('TM01') / sc.giga

    # print(round(fc_te11, 2))
    # print(round(fc_tm01, 2))

    assert abs(fc_te11 - 12.19) < 0.025
    assert abs(fc_tm01 - 15.92) < 0.025


if __name__ == "__main__":
    test_teflon_filled_waveguide()
    test_teflon_filled_circular_waveguide()
