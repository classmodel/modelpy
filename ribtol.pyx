#cython: boundscheck=False
#cython: wraparound=False

from libc.math cimport atan, log, exp, fabs

cdef double psim(double zeta):
    cdef double x, psim

    if(zeta <= 0):
        x     = (1. - 16. * zeta)**(0.25)
        psim  = 3.14159265 / 2. - 2. * atan(x) + log((1. + x)**2. * (1. + x**2.) / 8.)
    else:
        psim  = -2./3. * (zeta - 5./0.35) * exp(-0.35 * zeta) - zeta - (10./3.) / 0.35
    return psim
      
cdef double psih(double zeta):
    if(zeta <= 0):
        x     = (1. - 16. * zeta)**(0.25)
        psih  = 2. * log( (1. + x*x) / 2.)
    else:
        psih  = -2./3. * (zeta - 5./0.35) * exp(-0.35 * zeta) - (1. + (2./3.) * zeta) ** (1.5) - (10./3.) / 0.35 + 1.
    return psih

def ribtol(double Rib, double zsl, double z0m, double z0h): 
    cdef double L, L0, fx, Lstart, Lend, fxdif

    if(Rib > 0.):
        L    = 1.
        L0   = 2.
    else:
        L  = -1.
        L0 = -2.
    
    while (fabs(L - L0) > 0.001):
        L0      = L
        fx      = Rib - zsl / L * (log(zsl / z0h) - psih(zsl / L) + psih(z0h / L)) / (log(zsl / z0m) - psim(zsl / L) + psim(z0m / L))**2.
        Lstart  = L - 0.001*L
        Lend    = L + 0.001*L
        fxdif   = ( (- zsl / Lstart * (log(zsl / z0h) - psih(zsl / Lstart) + psih(z0h / Lstart)) / \
                                      (log(zsl / z0m) - psim(zsl / Lstart) + psim(z0m / Lstart))**2.) \
                  - (-zsl /  Lend   * (log(zsl / z0h) - psih(zsl / Lend  ) + psih(z0h / Lend  )) / \
                                      (log(zsl / z0m) - psim(zsl / Lend  ) + psim(z0m / Lend  ))**2.) ) / (Lstart - Lend)
        L       = L - fx / fxdif

        if(fabs(L) > 1e15):
            break

    return L
