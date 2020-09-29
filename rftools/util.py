"""Utilities for the RF-tools package."""

import numpy as np 
import scipy.constants as sc 


# Conversions ----------------------------------------------------------------

def dbm_to_w(dbm):
    """Convert dBm to W."""

    return 10 ** (dbm / 10.) * sc.milli


def w_to_dbm(w):
    """Convert W to dBm."""

    return 10 * np.log10(w / sc.milli)
    

# Functions for printing to terminal -----------------------------------------

def pvalf(name, val, units='', comment=''):
    """Print name, value as float, and units to terminal.

    Args:
        name (str): variable name
        val (float): variable value
        units (str): variable units (optional)
        comment (str): comment (optional)

    """

    if units != '':
        units = '\t[' + units + ']'
    if comment != '':
        comment = "  # {0}".format(comment)

    if isinstance(val, complex):
        re = val.real
        im = val.imag
        if val.imag < 0:
            str_tmp = "\t{0:20s}  {1:8.3f} - j{2:7.3f}{3:20s}{4}"
        else:
            str_tmp = "\t{0:20s}  {1:8.3f} + j{2:7.3f}{3:20s}{4}"
        print((str_tmp.format(name, re, abs(im), units, comment)))

    else:
        if val < 1000:
            str_tmp = "\t{0:20s} {1:7.3f}\t{2:20s}{3}"
            print((str_tmp.format(name, val, units, comment)))
        else:
            number = "{:7.3e}".format(val)
            number = number[:5] + ' ' + number[5:].upper()
            str_tmp = "\t{0:20s}   {1}{2:20s}{3}"
            print((str_tmp.format(name, number, units, comment)))

def header(header_string):
    """Print a nice header to the terminal.

    Args:
        header_string (str): Header title to print

    """

    print(("\t" + header_string))
    print(("\t" + "-" * 50))
