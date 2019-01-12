#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 12 10:46:20 2018

@author: vsc42247
"""



# purpose of calc_cm_ch: calculate momentum and thermal turbulent diffusion coefficients of the surface layer with a non-iterative procedure (Wouters et al., 2012)

# input:

# zrib = bulk Richardson number = (g/T)* DT * z/(ua^2)
#   with:
#     g = 9.81 m/s2 the gravitational acceleration
#     z = height (in meters) of the surface layer under consideration 
#     T = (reference) temperature (in Kelvin) at height z 
#     DT = (T - T_s) = temperature (in Kelvin) gradient between the surface and height z 
#     u_a^2 = u^2 +  v^2 is the squared horizontal absolute wind speed 
# zzz0m = ratio z/z0 between the height z and the momentum roughness length z0m
# zkbm = ln(z0m/z0h), with z0m, z0h the momentum and thermal roughness length, respectively.

# output: diffusion coefficients (CM and CH) which cna be used to determine surface-layer turbulent transport
# u'w' = - CM ua^2.
# w'T' = - CH ua DT 


# Reference:
# Wouters, H., De Ridder, K., and Lipzig, N. P. M.: Comprehensive
# Parametrization of Surface-Layer Transfer Coefficients for Use
# in Atmospheric Numerical Models, Bound.-Lay. Meteorol., 145,
# 539â550, doi:10.1007/s10546-012-9744-3, 2012.

import numpy as np

def calc_cm_ch (zeta,zzz0m,zkbm):
    krm = 0.4

    #ZETA = zeta_hs2(zrib,zzz0m,zkbm)
    FUNM,FUNH = funcsche(ZETA,zzz0m,zkbm)
    CM = krm**2.0/FUNM/FUNM
    CH = krm**2.0/FUNM/FUNH

    # FUNMn,FUNHn = funcsche(0.,zzz0m,zkbm)
    # CMn = krm**2.0/FUNMn/FUNMn
    # CHn = krm**2.0/FUNMn/FUNHn

    # print ZETA,FUNM,FUNH
    # print 'CMCMN',CM/CMn
    # print 'CHCHN',CH/CHn

    return CM,CH


def zeta_hs2(RiB,zzz0m,kBmin1):
    #print(RiB,zzz0m,kBmin1)
    mum=2.59
    muh=0.95
    nu=0.5
    lam=1.5

    betah = 5.0

    zzz0h = zzz0m*np.exp(kBmin1)
    zzzs = zzz0m*0.06 # to be changed!! r. 101 nog bekijken!!

    L0M = np.log(zzz0m)
    L0H = np.log(zzz0h)
    facM = np.log(1.+lam/mum/zzzs)*np.exp(-mum*zzzs)/lam
    facH = np.log(1.+lam/muh/zzzs)*np.exp(-muh*zzzs)/lam
    L0Ms = L0M + facM 
    L0Hs = L0H + facH

    if RiB < 0.:
        p = np.log(1.-RiB)
        Q = -0.486 +0.219*p - 0.0331*p**2-4.93*np.exp(-L0H) - 3.65/L0H +\
            0.38*p/L0H+ 14.8/L0H/L0H-0.946*p/L0H/L0H-10.0/L0H**3+ \
            0.392*L0M/L0H-0.084*p*L0M/L0H+0.368*L0M/L0H/L0H
        # print 'p: ',p
        # print 'Q: ',Q
        zeta = (1. + p*Q)* L0Ms**2/L0Hs * RiB
    else:
        betam = 4.76+7.03/zzz0m +0.24*zzz0m/zzz0h # to be changed
        # betam = 5.0 + 1.59*10.**(-5.)*(np.exp(13.0-L0M)-1.0) \
        #         +0.24*(np.exp(-kBmin1)-1.0) # to be changed!!
        # print('betam',betam)
        lL0M = np.log(L0M)
        S0Ms = 1.-1./zzz0m + (1.+nu/mum/zzzs)*facM
        S0Hs = 1.-1./zzz0h + (1.+nu/muh/zzzs)*facH
        zetat = -0.316-0.515*np.exp(-L0H) + 25.8 *np.exp(-2.*L0H) + 4.36/L0H \
                -6.39/L0H/L0H+0.834*lL0M - 0.0267*lL0M**2
        # print('zetat',zetat)
        RiBt = zetat *(L0Hs+ S0Hs*betah*zetat)/(L0Ms+S0Ms*betam*zetat)**2 
        # print('RiBt',RiBt)

        if (RiB > RiBt):
            D = (L0Ms+S0Ms*betam*zetat)**3/\
                (L0Ms*L0Hs+zetat*(2.*S0Hs * betah * L0Ms - S0Ms*betam*L0Hs))
            zeta = zetat + D*(RiB-RiBt)
        else:
            r = RiB - S0Hs*betah/(S0Ms*betam)**2
            B = S0Ms*betam*L0Hs- 2.*S0Hs*betah*L0Ms
            C = 4.*(S0Ms*betam)**2 * L0Ms *(S0Hs*betah*L0Ms-S0Ms*betam*L0Hs)
            zeta = - L0Ms / S0Ms/betam - B*C/(4.*(S0Ms*betam)**3 *(B**2+abs(C*r)))
            if r != 0:
                zeta = zeta + (B-np.sqrt(B**2+C*r) + B*C*r/(2.*(B**2+abs(C*r))))/(2.*(S0Ms*betam)**3*r)
    # print('zeta',zeta)
    return zeta

def funcsche(zeta,zzz0,kBmin1):


    mum=2.5
    muh=0.9
    nu=0.5
    lam=1.5
    
    p2=3.141592/2.
    
    lnzzz0=np.log(zzz0)
    zzzs=zzz0*0.06
    zetamcorr=(1.+nu/(mum*zzzs))*zeta
    zetam0=zeta/zzz0
    zetahcorr=(1.+nu/(muh*zzzs))*zeta
    zetah0=zeta/(zzz0*np.exp(kBmin1))
    
    if (zeta <= 0.):
    
        gamma=15.2
        alfam=0.25
        xx=(1.-gamma*zeta)**alfam
        psim=2.*np.log((1.+xx)/2.)+np.log((1.+xx**2.)/2.)-2.*np.arctan(xx)+p2
        xx0=(1.-gamma*zetam0)**alfam
        psim0=2.*np.log((1.+xx0)/2.)+np.log((1.+xx0**2.)/2.)-2.*np.arctan(xx0)+p2
        phimcorr=(1.-gamma*zetamcorr)**(-alfam)
        
        alfah=0.5
        yy=(1.-gamma*zeta)**alfah
        psih=2.*np.log((1.+yy)/2.)
        yy0=(1.-gamma*zetah0)**alfah
        psih0=2.*np.log((1.+yy0)/2.)
        phihcorr=(1.-gamma*zetahcorr)**(-alfah)
    else: 
    
        aa=6.1
        bb=2.5
        psim=-aa*np.log(zeta+(1.+zeta**bb)**(1./bb))
        psim0=-aa*np.log(zetam0+(1.+zetam0**bb)**(1./bb))
        phimcorr=1.+aa*(zetamcorr+zetamcorr**bb*(1.+zetamcorr**bb)**((1.-bb)/bb))/(zetamcorr+(1.+zetamcorr**bb)**(1./bb))
        
        cc=5.3
        dd=1.1
        psih=-cc*np.log(zeta+(1.+zeta**dd)**(1./dd))
        psih0=-cc*np.log(zetah0+(1.+zetah0**dd)**(1./dd))
        phihcorr=1.+cc*(zetahcorr+zetahcorr**dd*(1.+zetahcorr**dd)**((1.-dd)/dd))/(zetahcorr+(1.+zetahcorr**dd)**(1./dd))
    
    psistrm=phimcorr*(1./lam)*np.log(1.+lam/(mum*zzzs))*np.exp(-mum*zzzs)
    psistrh=phihcorr*(1./lam)*np.log(1.+lam/(muh*zzzs))*np.exp(-muh*zzzs)
    
    funm=lnzzz0-psim+psim0 +psistrm
    funh=lnzzz0+kBmin1-psih+psih0 +psistrh
    return funm,funh

