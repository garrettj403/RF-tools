RF tools
========

*Tools for designing RF components and networks*

Installation
------------

```bash
# for latest version (from GitHub)
python3 -m pip install git+https://github.com/garrettj403/RF-tools.git

# for lastest release (from PyPI)
python3 -m pip install rftools
```

Examples of the Command Line Tools
----------------------------------

**Calculate the properties of a WR4.3 rectangular waveguide at 230 GHz:**

Input:
```bash
$ waveguide WR4.3 --freq 230
```
Output:
```bash
    Rectangular Waveguide: WR4.3
    --------------------------------------------------

    Dimensions:
    a                      1.092        [mm]
    b                      0.546        [mm]

    Standard frequency range:
    low                  171.553        [GHz]
    mid                  215.471        [GHz]
    high                 259.388        [GHz]

    Cutoff frequencies:
    TE10                 137.242        [GHz]
    TE20                 274.485        [GHz]
    TE01                 274.485        [GHz]
    TE/TM11              306.883        [GHz]

    Properties at 230.0 GHz:
    wavelength             1.624        [mm]
    impedance            469.469        [ohms]
```

**Calculate the attenuation constant of a WR2.8 waveguide at 345 GHz:**

Input
```bash
$ waveguide-att --type WR2.8 --freq 345 --cond 5.85e7
```
Output
```bash
    Rectangular Waveguide: WR2.8
    --------------------------------------------------

    Dimensions:
    a                    711.200        [um]
    b                    355.600        [um]

    Standard frequency range:
    low                  263.457        [GHz]
    mid                  330.901        [GHz]
    high                 398.346        [GHz]

    Properties at 345 GHz:
    wavelength             1.098        [mm]
    impedance            475.852        [ohms]

    Attenuation at 345 GHz:
    conductivity           5.850 E+07   [S/m]
    skin depth           112.030        [nm]
    attenuation            1.976        [Np/m]
                          17.160        [dB/m]
                           0.172        [dB/cm]
``` 
**Calculate the properties of a 0.5 mm radius circular waveguide at 345 GHz:**

Input: 
```bash
$ cwaveguide 0.5 --freq 345
```
Output:
```bash
    Circular Waveguide:
    --------------------------------------------------

    Dimensions:
    radius a               0.500        [mm]

    Cutoff frequencies:
    TE11                 175.681        [GHz]
    TM01                 229.502        [GHz]
    TE21                 291.434        [GHz]
    TE01                 365.676        [GHz]
    TM11                 365.676        [GHz]

    Properties at 345.0 GHz:
    wavelength             1.010        [mm]
    impedance            437.735        [ohms]
```

**Calculate the noise temperature using the Y-factor technique:**

Input:
```bash
$ noisetemp 4.9/2.2 --freq 850 --thot 300 --tcold 20
```
Output:
```bash
    Noise temperature from Y-factor
    --------------------------------------------------

    Physical temperature of black body loads:
    Hot load             300.000        [K]
    Cold load             20.000        [K]

    Equiv. temp. from CW equations (with f=850.0 GHz):
    Hot load             300.462        [K]
    Cold load             26.496        [K]

    Y-factor               2.227

    Noise temperature    196.735        [K]
```

**Calculate the width of a microstrip:**

Input:
```bash
$ 50ohm-line --z0 50 --thickness 15 --er 2.2
```
Output:
```bash
    Microstrip:
    --------------------------------------------------

    Input values:
    desired Z0            50.000        [ohms]
    thickness (t)         15.000        [mil]
    rel. permittivity      2.200

    Output:
    microstrip width       1.174        [mm]
```