"""Class to handle S-parameters."""

import numpy as np 
from numpy.linalg import inv
import matplotlib.pyplot as plt 
from copy import deepcopy as copy
from rftools.parameters import *


class Network(object):
    """Class to handle S-parameters.

    Args:
        filename (str): Touchstone file to load

    Keyword Args:
        comment (str): comment to describe instance

    """

    def __init__(self, filename=None, **kwargs):

        self.comment = kwargs.get('comment', '')

        if filename is not None:
            self.f, self.s, self.z0, self.ports = _read_touchstone(filename)
        else:
            self.f = kwargs.get('f', None) 
            self.s = kwargs.get('s', None) 
            self.z0 = kwargs.get('z0', None) 
            self.ports = kwargs.get('ports', ())

    ### Dunder ###

    def __str__(self):

        msg = '<Network: {} ports, {} points, {}>'
        return msg.format(self.count_ports(), self.count_points(), self.ports)

    def __repr__(self):

        return self.__str__()

    def __mul__(self, network2):
        """Cascade S-parameter matrices."""

        network1 = self.copy()

        assert isinstance(network2, Network), \
            "Only Networks can be multiplied together."
        assert np.all(network1.f == network2.f), \
            "Frequencies must match."
        assert network1._is_2port(), \
            "First network must have 2 ports."
        assert network2._is_1port() or network2._is_2port(), \
            "Second network must have 1 or 2 ports."

        if network2._is_1port():

            _sparam2 = np.ones((2, 2, np.alen(self.f)), dtype=complex) * 1e-15
            _sparam2[0, 0, :] = network2.s[0, 0]

            _z02 = np.ones((2, np.alen(self.f)), dtype=complex) * network2.z0[0]

            abcd1 = s_to_abcd(network1.s, network1.z0)
            abcd2 = s_to_abcd(_sparam2, _z02)

            abcd = np.empty_like(abcd1)
            for i in range(abcd1.shape[2]):
                abcd[:,:,i] = np.dot(abcd1[:,:,i], abcd2[:,:,i])

            z01 = network1.z0[0]
            z02 = _z02[1]
            z0 = np.vstack((z01, z02))
            s = abcd_to_sparam(abcd, z0)

            result = Network()
            result.f = network1.f 
            result.s = np.zeros((1, 1, np.alen(self.f)), dtype=complex)
            result.s[0, 0, :] = s[0, 0]
            result.z0 = np.zeros((1, np.alen(self.f)), dtype=float)
            result.z0[0, :] = network1.z0[0] 
            
            result.ports = (network1.ports[0])

            return result

        if network2._is_2port():

            abcd1 = s_to_abcd(network1.s, network1.z0)
            abcd2 = s_to_abcd(network2.s, network2.z0)

            abcd = np.empty_like(abcd1)
            for i in range(abcd1.shape[2]):
                abcd[:,:,i] = np.dot(abcd1[:,:,i], abcd2[:,:,i])

            z01 = network1.z0[0]
            z02 = network2.z0[1]
            z0 = np.vstack((z01, z02))
            s = abcd_to_sparam(abcd, z0)

            result = Network()
            result.f = network1.f 
            result.s = s 
            result.z0 = z0 
            
            result.ports = (network1.ports[0], network2.ports[1])

            return result

    ### Meta ###

    def copy(self):
        """Do a deep copy."""

        return copy(self)

    def count_ports(self):

        return self.s.shape[0]

    def count_points(self):

        return self.s.shape[2]

    def check_ports(self):

        sn, sm, _ = self.s.shape
        assert sn == sm, \
            "S-parameters not square."
        nports = sn 

        assert nports == self.z0.shape[0]

        assert len(self.ports) == nports

    ### Checks ###

    def _is_1port(self):

        return self.count_ports() == 1

    def _is_2port(self):

        return self.count_ports() == 2

    def _is_normalised(self):

        return np.array_equal(self.z0, np.ones_like(self.z0)*50.)

    def is_reciprocal(self):

        return np.allclose(self.s[0, 1], self.s[1, 0])

    def is_symmetrical(self):
        
        return np.allclose(self.s[0, 0], self.s[1, 1])

    def is_lossless_reciprocal(self):

        test1 = np.allclose(np.abs(self.s[0, 0]), np.abs(self.s[1, 1]))
        test2 = np.allclose(np.abs(self.s[0, 0])**2 + np.abs(self.s[0, 1])**2, np.ones_like(self.f))

        return test1 & test2

    ### Get other parameters ###

    def tparam(self):

        return s_to_tparam(self.s)

    def zparam(self):

        return s_to_zparam(self.s, self.z0)

    def sparam_db(self):

        with np.errstate(divide='ignore'):
            s_db = 20*np.log10(np.abs(self.s))

        return s_db

    ### Port impedance ###

    # def renormalise_ports(self, z0new=None):
    #     # http://qucs.sourceforge.net/tech/node98.html

    #     print "THIS FUNCTION DOESNT WORK"

    #     if z0new is None:
    #         z0new = np.ones_like(self.z0) * 50.

    #     _, _, npts = self.s.shape
    #     nports = self.count_ports()

    #     for i in range(npts):

    #         S = np.matrix(self.s[:,:,i])

    #         Znbefore = self.z0[:,i]
    #         Zn = z0new[:,i]

    #         r = (Zn - Znbefore) / (Zn + Znbefore)
    #         R = np.diag(r)
    #         print R

    #         a = np.sqrt(Zn / Znbefore) / (Zn + Znbefore)
    #         A = np.diag(a)

    #         self.s[:,:,i] = inv(A) * (S - R) * inv(np.identity(nports) - R * S) * A
    #         self.z0[:,i] = Zn

    ### Manipulate ports ###

    def port_number(self, port_name):

        for i in range(len(self.ports)):
            if self.ports[i] == port_name:
                return i 
        raise ValueError

    def delete_port(self, port):

        port_num = self.port_number(port)

        mask = np.arange(len(self.ports)) != port_num
        s = self.s[mask]
        self.s = s[:,mask]
        self.z0 = self.z0[mask]

        ports = list(self.ports)
        ports.remove(port)
        self.ports = tuple(ports)

    def rename_port(self, old_name, new_name):

        ports = list(self.ports)

        new_ports = []
        for port in ports:
            if port == old_name:
                new_ports.append(new_name)
            else:
                new_ports.append(port)

        self.ports = tuple(new_ports)

    def flip_ports(self):

        npts = self.count_points()
        assert self._is_2port()

        # Flip s-parameters
        for i in range(npts):
            mat = np.array([[0., 1.], [1., 0.]])
            self.s[:,:,i] = np.dot(np.dot(mat, self.s[:,:,i]), mat)

        # Flip Z0
        z01 = self.z0[0] 
        z02 = self.z0[1]
        self.z0 = np.vstack((z02, z01))

        # Flip port names
        self.ports = (self.ports[1], self.ports[0])

    ### Truncate data ###

    def freq_range(self, fmin, fmax):

        mask = (fmin <= self.f) & (self.f <= fmax)

        self.f = self.f[mask]
        self.s = self.s[:, :, mask]
        self.z0 = self.z0[:, mask]

    ### Get property from given port name ### 

    def get_s(self, porta, portb):

        a = self.port_number(porta)
        b = self.port_number(portb)

        return self.s[a, b]

    def get_s_db(self, porta, portb, sigma=None):

        a = self.port_number(porta)
        b = self.port_number(portb)

        s_db = self.sparam_db()[a, b]

        if sigma is not None:
            sigma = sigma / (self.f[1] - self.f[0])
            s_db = _gauss_conv(s_db, sigma)

        return s_db 

    def get_s_mag(self, porta, portb, sigma=None):

        a = self.port_number(porta)
        b = self.port_number(portb)

        s_mag = np.abs(self.s[a, b])

        if sigma is not None:
            sigma = sigma / (self.f[1] - self.f[0])
            s_mag = _gauss_conv(s_mag, sigma)

        return s_mag

    def get_p_mag(self, porta, portb, sigma=None):

        a = self.port_number(porta)
        b = self.port_number(portb)

        s_mag = np.abs(self.s[a, b])

        if sigma is not None:
            sigma = sigma / (self.f[1] - self.f[0])
            s_mag = _gauss_conv(s_mag, sigma)

        return s_mag**2

    def get_gain(self, porta, portb, sigma=None):

        a = self.port_number(porta)
        b = self.port_number(portb)

        s_mag = np.abs(self.s[a, b])

        if sigma is not None:
            sigma = sigma / (self.f[1] - self.f[0])
            s_mag = _gauss_conv(s_mag, sigma)

        return s_mag**2

    def get_z0(self, port):

        a = self.port_number(port)

        return self.z0[a]

    def get_zin(self, port, sigma=None):
        
        a = self.port_number(port)
        zin = (1 + self.s[a, a]) / (1 - self.s[a, a]) * self.z0[a]

        if sigma is not None:
            sigma = sigma / (self.f[1] - self.f[0])
            zreal = _gauss_conv(zin.real, sigma)
            zimag = _gauss_conv(zin.imag, sigma)
            zin = zreal + 1j * zimag

        return zin

    ### Get value at a certain frequency ###

    def get_zin_at_f(self, f, port, sigma=None):
        
        a = self.port_number(port)
        zin = (1 + self.s[a, a]) / (1 - self.s[a, a]) * self.z0[a]

        if sigma is not None:
            sigma = sigma / (self.f[1] - self.f[0])
            zreal = _gauss_conv(zin.real, sigma)
            zimag = _gauss_conv(zin.imag, sigma)
            zin = zreal + 1j * zimag

        return np.interp(f, self.f, zin)

    def get_z0_at_f(self, f, port):
        
        a = self.port_number(port)
        z0 = self.z0[a]

        return np.interp(f, self.f, z0)

    def get_s_mag_at_f(self, f, porta, portb, sigma=None):

        s_mag = self.get_s_mag(porta, portb, sigma=None)

        return np.interp(f, self.f, s_mag)

    def get_gain_at_f(self, f, porta=None, portb=None, sigma=None):

        if porta is None or portb is None:
            porta = self.ports[1]
            portb = self.ports[0]

        il = self.get_s_mag(porta, portb, sigma=None)

        return np.interp(f, self.f, il)**2

    ### Plot network properties ###

    def plot_sparam(self, filename=None, ax=None, **params):

        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.get_figure()

        for i, key1 in enumerate(self.ports):
            for j, key2 in enumerate(self.ports):
                if i == 1:
                    ax.plot(self.f, self.get_s_db(key1, key2), label="S({},{})".format(key1, key2), ls='--')
                else:
                    ax.plot(self.f, self.get_s_db(key1, key2), label="S({},{})".format(key1, key2))
        ax.legend()
        ax.set(**params)
        ax.set_xlabel('Frequency (GHz)')
        ax.set_ylabel('Magnitude (dB)')

        if filename is not None:
            fig.savefig(fig_name, bbox_inches='tight')
            plt.close(fig)
            return
        else:
            return ax

    def plot_zin(self, port, filename=None, ax=None, **params):
    
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.get_figure()
            
        zin = self.get_zin(port)
        ax.plot(self.f, zin.real, 'b', label='Real')
        ax.plot(self.f, zin.imag, 'r', label='Imaginary') 
        ax.set_xlabel(r'Frequency (GHz)')
        ax.set_ylabel(r'Impedance ($\Omega$)')
        ax.legend(title=r'$Z_\mathrm{{in}}$ at {}'.format(port))
        ax.set(**params)

        if filename is not None:
            fig.savefig(fig_name, bbox_inches='tight')
            plt.close(fig)
            return
        else:
            return ax

    def plot(self, filename=None, ax=None, **params):

        return self.plot_sparam(filename, ax=None, **params)


# Read touchstone ------------------------------------------------------------

def _read_touchstone(filename, max_size=10000, z0=None):
    """Read touchstone from HFSS."""

    with open(filename, 'r') as fin:
        data = fin.readlines()

    # Get header
    for i, line in enumerate(data):
        if line[0] != '!' and line[0] != '#':
            break
    header = data[:i]
    data = data[i:]

    # # Get data type
    # assert data[0] == '# GHZ S MA\n', "Only GHz, S, MA supported"
    # data = data[1:]

    # # Divide between header and data
    # for i, line in enumerate(data[:10]):
    #     if line[:1] != '!' and line[:1] != '#':
    #         break

    # header = data[:i]
    
    # # # Check if normalised
    # # print header[0].split()[-2:]

    # data = data[i+1:]

    for i in range(len(header)):
        header[i] = header[i].rstrip()
    for i in range(len(data)):
        data[i] = data[i].rstrip()

    # READ HEADER ------------------------------------------------------------

    # Get ports
    ports = []
    for line in header:
        if 'Port[' in line:
            port_name = line.split()[-1]
            port_name = port_name.split(':')[0]
            ports.append(port_name)
    ports = tuple(ports)
    nports = len(ports) 
    # print(nports)

    # Initialize S-parameters
    f = np.zeros(max_size, dtype=float)
    s = np.zeros((nports, nports, max_size), dtype=complex)

    # READ DATA --------------------------------------------------------------

    # Remove carriage returns in data
    new_data = []
    for line in data:
        if line == '' or line[0] != ' ':
            new_data.append(line)
        elif line[0] == ' ':
            new_data[-1] += line 
    data = new_data

    # Read data lines
    idx = 0
    for line in data:
        if len(line) > 0 and line[0] != '!':
            dat = [float(entry) for entry in line.split()]
            f[idx] = dat[0]
            dat = dat[1:]
            mag = dat[::2]
            arg = dat[1::2]
            for i in range(nports):
                for j in range(nports):
                    s[i,j,idx] = mag.pop(0) * np.exp(1j * _deg_to_rad(arg.pop(0)))
            idx += 1

    # Read other data (z0)
    if z0 is None:
        idx = 0
        z0 = np.zeros((nports, max_size))
        for line in data:
            if len(line) > 0 and line[0:10] == '! Port Imp':
                # TODO: read imaginary component as well!
                z0[:, idx] = [float(entry) for entry in line[16:].split()[::2] if _is_float(entry)]
                idx += 1
    else:
        z0 = z0 * np.ones((nports, max_size))

    # Format output data
    mask = f != 0
    f = f[mask]
    s = s[:,:,mask]
    z0 = z0[:,mask]

    # Check matrix dimensions
    nports1, nports2, nfreq = s.shape
    assert nports1 == nports2, "S-matrix must be square."
    assert nfreq == np.alen(f)
    assert nfreq == np.alen(z0[0])
    assert nports1 == np.alen(z0[:,0])
    
    # If 2-port, exchange s12 and s21 (historical reasons...)
    if nports1 == 2:
        s[0, 1], s[1, 0] = s[1, 0], s[0, 1]

    return f, s, z0, ports


# Helper functions -----------------------------------------------------------

def _is_float(entry):
    try:
        float(entry)
        return True
    except ValueError:
        return False 


def _deg_to_rad(deg):

    return deg * np.pi / 180.   


# Filters -------------------------------------------------------------------

def _gauss_conv(x, sigma=10, ext_x=3):
    """Gaussian convolution filter

    Args:
        x (ndarray): noisy data
        sigma (float): standard deviation of Gaussian curve
        ext_x (float): Gaussian curve will extend from ext_x * sigma in each
                       direction

    Returns:
        ndarray: filtered data

    """

    wind = _gauss(sigma, ext_x)
    wlen = np.alen(wind)

    assert wlen <= np.alen(x), "Window size must be smaller than data size"
    assert sigma * ext_x >= 1, \
        "Window size must be larger than 1. Increase ext_x."

    s = np.r_[x[wlen - 1:0:-1], x, x[-2:-wlen - 1:-1]]
    y_out = np.convolve(wind / wind.sum(), s, mode='valid')
    y_out = y_out[wlen // 2:-wlen // 2 + 1]

    return y_out


def _gauss(sigma, n_sigma=3):
    """Discrete, normalized Gaussian centered on zero. Used for filtering data.

    Args:
        sigma (float): standard deviation
        n_sigma (float): extend x in each direction by ext_x * sigma

    Returns:
        ndarray: discrete Gaussian curve

    """

    x_range = n_sigma * sigma
    x = np.arange(-x_range, x_range + 1e-5, 1, dtype=float)

    y = 1 / (sigma * np.sqrt(2 * np.pi)) * np.exp(-0.5 * (x / sigma)**2)

    return y


# Main -----------------------------------------------------------------------

if __name__ == "__main__":
    
    network = Network('../workbooks/data//simple-waveguide.s2p')

    # data = Network('../tests/data/hfss-sparam.s3p')
    # print(data) 
    # print('No. ports: ', data.count_ports())
    # print('No. points: ', data.count_points())
    # data.check_ports()
    # print("")

    # data.delete_port('IF')
    # print(data)
    # print('2 port: ', data._is_2port())
    # print('Normalised: ', data._is_normalised())
    # print('J port number: ', data.port_number('J'))
    # # data.plot_sparam(xlim=[180, 260], ylim=[-30, 0])
    # # data.plot_zin('J', ylim=[-20, 30], xlim=[200, 260])

    # data1 = data.copy()
    # data1.rename_port('J', 'J1')
    # data1.rename_port('WG', 'WG1')

    # data2 = data.copy()
    # data2.flip_ports()
    # data2.rename_port('J', 'J2')
    # data2.rename_port('WG', 'WG2')

    # new = data1 * data2 
    # # new.plot_sparam(xlim=[200, 260], ylim=[-15, 0])

    # # data.renormalise_ports()
    # # data.plot_sparam()
    # # print data.z0
