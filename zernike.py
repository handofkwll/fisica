# Zernike routines written by Tim van Werkhoven (werkhoven@strw.leidenuniv.nl)
# Subset of his module libtim-py, code reformatted.

import numpy as np
from scipy.misc import factorial as fac

def zernike_rad(m, n, rho):
    """
    Make radial Zernike polynomial on coordinate grid **rho**.
    @param [in] m Radial Zernike index
    @param [in] n Azimuthal Zernike index
    @param [in] rho Radial coordinate grid
    @return Radial polynomial with identical shape as **rho**
    """
    if (np.mod(n-m, 2) == 1):
        return rho*0.0

    wf = rho*0.0
    for k in range((n-m)/2+1):
        wf += rho**(n-2.0*k) * (-1.0)**k * fac(n-k) / \
          ( fac(k) * fac( (n+m)/2.0 - k ) * fac( (n-m)/2.0 - k ) )

    return wf

def zernike(m, n, rho, phi, norm=True):
    """
    Calculate Zernike mode (m,n) on grid **rho** and **phi**.
    **rho** and **phi** should be radial and azimuthal coordinate grids of identical shape, respectively.
    @param [in] m Radial Zernike index
    @param [in] n Azimuthal Zernike index
    @param [in] rho Radial coordinate grid
    @param [in] phi Azimuthal coordinate grid
    @param [in] norm Normalize modes to unit variance
    @return Zernike mode (m,n) with identical shape as rho, phi
    @see <http://research.opt.indiana.edu/Library/VSIA/VSIA-2000_taskforce/TOPS4_2.html> and <http://research.opt.indiana.edu/Library/HVO/Handbook.html>.
    """
    nc = 1.0
    if (norm):
        nc = (2*(n+1)/(1+(m==0)))**0.5
    if (m > 0):
        return nc*zernike_rad(m, n, rho) * np.cos(m * phi)
    if (m < 0):
        return nc*zernike_rad(-m, n, rho) * np.sin(-m * phi)

    return nc*zernike_rad(0, n, rho)
