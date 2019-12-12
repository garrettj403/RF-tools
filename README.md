RF-tools
========

Tools to help design RF components / networks.

Installation
------------

To install via pip:

```bash
python3 -m pip install git+https://github.com/garrettj403/RF-tools.git
```

If you would like to use the command line utilities, you'll have to clone the repo onto your computer:

```bash
git clone https://github.com/garrettj403/RF-tools.git
ln -s RF-tools/bin/* /usr/local/bin/
```

Example: Command Line
---------------------

**Analyze a rectangular waveguide:**

Input:
```bash
$ waveguide WR4.3 --freq 230
```
Output:
```bash
    Rectangular Waveguide: WR4.3
    --------------------------------------------------
    a               =   1.092       [mm]
    b               =   0.546       [mm]

    cutoff TE10     = 137.242       [GHz]
    cutoff TE20     = 274.485       [GHz]
    cutoff TE01     = 274.485       [GHz]
    cutoff TE/TM11  = 306.883       [GHz]

    -> at 230.0 GHz
    wavelength      =   1.624       [mm]
    impedance       = 469.469       [ohms]
```

**Calculate noise temperature:**

Input:
```bash
$ noisetemp 4.9/2.2 --freq 850 --thot 300 --tcold 20
```
Output:
```bash
    Noise temperature from Y-factor
    --------------------------------------------------

    Y-factor        =   2.227

    -> using CW equations
    Hot load        = 300.462       [K]
    Cold load       =  26.496       [K]

    Noise temp.     = 196.735       [K]
```

