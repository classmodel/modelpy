import numpy as np
from numba import jit

@jit(nopython=True)
def psim(zeta):
    if(zeta <= 0):
        x     = (1. - 16. * zeta)**(0.25)
        psim  = 3.14159265 / 2. - 2. * np.arctan(x) + np.log((1. + x)**2. * (1. + x**2.) / 8.)
        #x     = (1. + 3.6 * abs(zeta) ** (2./3.)) ** (-0.5)
        #psim = 3. * np.log( (1. + 1. / x) / 2.)
    else:
        psim  = -2./3. * (zeta - 5./0.35) * np.exp(-0.35 * zeta) - zeta - (10./3.) / 0.35
    return psim

@jit(nopython=True)
def psih(zeta):
    if(zeta <= 0):
        x     = (1. - 16. * zeta)**(0.25)
        psih  = 2. * np.log( (1. + x*x) / 2.)
        #x     = (1. + 7.9 * abs(zeta) ** (2./3.)) ** (-0.5)
        #psih  = 3. * np.log( (1. + 1. / x) / 2.)
    else:
        psih  = -2./3. * (zeta - 5./0.35) * np.exp(-0.35 * zeta) - (1. + (2./3.) * zeta) ** (1.5) - (10./3.) / 0.35 + 1.
    return psih

@jit(nopython=True)
def ribtol(Rib, zsl, z0m, z0h):
    if(Rib > 0.):
        L    = 1.
        L0   = 2.
    else:
        L  = -1.
        L0 = -2.

    while (abs(L - L0) > 0.001):
        L0      = L
        fx      = Rib - zsl / L * (np.log(zsl / z0h) - psih(zsl / L) + psih(z0h / L)) / (np.log(zsl / z0m) - psim(zsl / L) + psim(z0m / L))**2.
        Lstart  = L - 0.001*L
        Lend    = L + 0.001*L
        fxdif   = ( (- zsl / Lstart * (np.log(zsl / z0h) - psih(zsl / Lstart) + psih(z0h / Lstart)) / \
                                      (np.log(zsl / z0m) - psim(zsl / Lstart) + psim(z0m / Lstart))**2.) \
                  - (-zsl /  Lend   * (np.log(zsl / z0h) - psih(zsl / Lend  ) + psih(z0h / Lend  )) / \
                                      (np.log(zsl / z0m) - psim(zsl / Lend  ) + psim(z0m / Lend  ))**2.) ) / (Lstart - Lend)
        L       = L - fx / fxdif

        if(abs(L) > 1e15):
            break

    return L




