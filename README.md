RF tools
========

*Tools to design RF components and networks*

Installation
------------

```bash
# for latest version (from GitHub)
python3 -m pip install git+https://github.com/garrettj403/RF-tools.git

# for lastest release (from PyPI)
python3 -m pip install rftools
```

Examples: Command Line
----------------------

**Analyze a WR4.3 rectangular waveguide:**

Input:
```bash
$ waveguide WR4.3 --freq 230
```
Output:
```bash
    Rectangular Waveguide: WR4.3
    --------------------------------------------------
    a               =   1.092               [mm]
    b               =   0.546               [mm]

    low freq.       = 171.553               [GHz]
    mid freq.       = 215.471               [GHz]
    high freq.      = 259.388               [GHz]

    cutoff TE10     = 137.242               [GHz]
    cutoff TE20     = 274.485               [GHz]
    cutoff TE01     = 274.485               [GHz]
    cutoff TE/TM11  = 306.883               [GHz]

    -> at 230.0 GHz
    wavelength      =   1.624               [mm]
    impedance       = 469.469               [ohms]
```

**Analyze a 0.5 mm radius circular waveguide:**

Input: 
```bash
$ cwaveguide 0.5 --freq 345
```
Output:
```bash
    Circular Waveguide:
    --------------------------------------------------
    a               =   0.500               [mm]

    low freq.       = 219.601               [GHz]
    mid freq.       = 275.819               [GHz]
    high freq.      = 332.037               [GHz]

    cutoff TE11     = 175.681               [GHz]
    cutoff TM01     = 229.502               [GHz]
    cutoff TE21     = 291.434               [GHz]
    cutoff TE01     = 365.676               [GHz]
    cutoff TM11     = 365.676               [GHz]

    -> at 345.0 GHz
    wavelength      =   1.010               [mm]
    impedance       = 437.735               [ohms]
```

**Calculate noise temperature using Y-factor technique:**

Input:
```bash
$ noisetemp 4.9/2.2 --freq 850 --thot 300 --tcold 20
```
Output:
```bash
    Noise temperature from Y-factor
    --------------------------------------------------

    Y-factor        =   2.227

    -> using CW equations with f = 850 GHz
    Hot load        = 300.462       [K]
    Cold load       =  26.496       [K]

    Noise temp.     = 196.735       [K]
```
