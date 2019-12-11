import numpy as np 


# S-parameters to ... --------------------------------------------------------

def s_to_zparam(sparam, z0):

    # Method 1: 
    # https://en.wikipedia.org/wiki/Impedance_parameters#Relation_to_S-parameters

    # zparam = np.empty_like(sparam)
    # _, _, npts = sparam.shape

    # for idx in range(npts):
        
    #     zsqrt = np.matrix([[np.sqrt(z0[0,idx]), 0],
    #                        [0, np.sqrt(z0[1,idx])]])
    #     s = np.matrix(sparam[:,:,idx])
    #     i = np.matrix([[1., 0.], 
    #                    [0., 1.]])

    #     zparam[:,:,idx] = zsqrt * (i + s) * np.linalg.inv(i - s) * zsqrt

    # Method 2:
    # DOI: 10.1109/22.275248

    z01 = z0[0]
    z02 = z0[1]
    z01c = np.conj(z01)
    z02c = np.conj(z02)
    r01 = z01.real
    r02 = z02.real

    s11 = sparam[0,0]
    s21 = sparam[1,0]
    s12 = sparam[0,1]
    s22 = sparam[1,1]

    denom = (1 - s11) * (1 - s22) - s12 * s21

    z11 = ((z01c + s11 * z01) * (1 - s22) + s12 * s21 * z01) / denom
    z12 = (2 * s12 * (r01 * r02)**0.5) / denom
    z21 = (2 * s21 * (r01 * r02)**0.5) / denom
    z22 = ((z02c + s22 * z02) * (1 - s11) + s12 * s21 * z02) / denom
    
    zparam = np.empty_like(sparam)
    zparam[0,0] = z11
    zparam[1,0] = z21
    zparam[0,1] = z12
    zparam[1,1] = z22

    return zparam


def s_to_tparam(sparam):

    s11 = sparam[0,0]
    s21 = sparam[1,0]
    s12 = sparam[0,1]
    s22 = sparam[1,1]

    t11 = -(s11 * s22 - s12 * s21) / s21 
    t12 =  s11 / s21
    t21 = -s22 / s21
    t22 =  1   / s21 
    
    tparam = np.empty_like(sparam)
    tparam[0,0] = t11
    tparam[1,0] = t21
    tparam[0,1] = t12
    tparam[1,1] = t22

    return tparam


def s_to_zparam(sparam, z0):

    z01 = z0[0]
    z02 = z0[1]
    z01c = np.conj(z01)
    z02c = np.conj(z02)
    r01 = z01.real
    r02 = z02.real

    s11 = sparam[0,0]
    s21 = sparam[1,0]
    s12 = sparam[0,1]
    s22 = sparam[1,1]

    z11 = (z01c + s11 * z01) * (1 - s22) + s12 * s21 * z01
    z12 = 2 * s12 * np.sqrt(r01 * r02)
    z21 = 2 * s21 * np.sqrt(r01 * r02)
    z22 = (z02c + s22 * z02) * (1 - s11) + s12 * s21 * z02

    denom = (1 - s11) * (1 - s22) - s12 * s21

    zparam = np.empty_like(sparam)
    zparam[0, 0] = z11 / denom
    zparam[0, 1] = z12 / denom
    zparam[1, 0] = z21 / denom
    zparam[1, 1] = z22 / denom

    # import matplotlib.pyplot as plt 
    # # plt.plot(denom)
    # plt.plot(np.abs(s12 * s21))
    # plt.plot(np.abs((1 - s11)*(1 - s22)))
    # # plt.plot(rf_cct.f, np.abs(z11))
    # # plt.plot(rf_cct.f, np.abs(z12))
    # # plt.plot(rf_cct.f, np.abs(z21))
    # # plt.plot(rf_cct.f, np.abs(z22))
    # plt.show()

    return zparam 


def s_to_abcd(sparam, z0):

    z01 = z0[0]
    z02 = z0[1]
    z01c = np.conj(z01)
    z02c = np.conj(z02)
    r01 = z01.real
    r02 = z02.real

    s11 = sparam[0,0]
    s21 = sparam[1,0]
    s12 = sparam[0,1]
    s22 = sparam[1,1]

    denom = 2 * s21 * np.sqrt(r01 * r01)
    A = ((z01c + s11 * z01) * (1 - s22) + s12 * s21 * z01) / denom
    B = ((z01c + s11 * z01) * (z02c + s22 * z02)  - s12 * s21 * z01 * z02) / denom
    C = ((1 - s11) * (1 - s22) - s12 * s21) / denom
    D = ((1 - s11) * (z02c +  s22 * z02) + s12 * s21 * z02) / denom

    abcd = np.empty_like(sparam)
    abcd[0, 0] = A 
    abcd[0, 1] = B
    abcd[1, 0] = C
    abcd[1, 1] = D 

    return abcd 


# ... to S-parameters --------------------------------------------------------

def z_to_sparam(zparam, z0):

    z11 = zparam[0, 0]
    z12 = zparam[0, 1]
    z21 = zparam[1, 0]
    z22 = zparam[1, 1]

    z01 = z0[0]
    z02 = z0[1]
    z01c = np.conj(z01)
    z02c = np.conj(z02)
    r01 = z01.real
    r02 = z02.real

    s11 = (z11 - z01c) * (z22 + z02) - z12 * z21 
    s12 = 2 * z12 * np.sqrt(r01 * r02)
    s21 = 2 * z21 * np.sqrt(r01 * r02)
    s22 = (z11 + z01) * (z22 - z02c) - z12 * z21 

    denom = (z11 + z01) * (z22 + z02) - z12 * z21

    sparam = np.empty_like(zparam)
    sparam[0, 0] = s11 / denom
    sparam[0, 1] = s12 / denom
    sparam[1, 0] = s21 / denom
    sparam[1, 1] = s22 / denom

    return sparam


def abcd_to_sparam(abcd, z0):

    z01 = z0[0]
    z02 = z0[1]
    z01c = np.conj(z01)
    z02c = np.conj(z02)
    r01 = z01.real
    r02 = z02.real

    A = abcd[0, 0]
    B = abcd[0, 1]
    C = abcd[1, 0]
    D = abcd[1, 1]

    denom = A * z02 + B + C * z01  * z02 + D * z01
    s11 = (A * z02 + B - C * z01c * z02 - D * z01c) / denom 
    s12 = (2 * (A * D - B * C) * np.sqrt(r01 * r02)) / denom
    s21 = (2 * np.sqrt(r01 * r02)) / denom
    s22 = (-A * z02c + B - C * z01 * z02c + D * z01) / denom

    sparam = np.empty_like(abcd)
    sparam[0,0] = s11
    sparam[0,1] = s12
    sparam[1,0] = s21
    sparam[1,1] = s22

    return sparam


def t_to_sparam(tparam):

    t11 = tparam[0,0]
    t21 = tparam[1,0]
    t12 = tparam[0,1]
    t22 = tparam[1,1]

    s11 =  t12 / t22
    s12 = (t11 * t22 - t12 * t21) / t22
    s21 =  1   / t22
    s22 = -t21 / t22 
    
    sparam = np.empty_like(tparam)
    sparam[0,0] = s11
    sparam[1,0] = s21
    sparam[0,1] = s12
    sparam[1,1] = s22

    return sparam
