#
# Example of how to run the Python code, and access the output
# This case is identical to the default setup of CLASS (the version with interface) 
#

import matplotlib.pyplot as pl
import numpy as np
import copy

from model import *

""" 
Create empty model_input and set up ARM-SGP case
"""
arm = model_input()

arm.dt         = 60.       # time step [s]
arm.runtime    = 12*3600    # total run time [s]

# mixed-layer input
arm.sw_ml      = True      # mixed-layer model switch
arm.sw_shearwe = False     # shear growth mixed-layer switch
arm.sw_fixft   = False     # Fix the free-troposphere switch
arm.h          = 140.      # initial ABL height [m]
arm.Ps         = 97000     # surface pressure [Pa]
arm.divU       = 0.        # horizontal large-scale divergence of wind [s-1]
arm.fc         = 1.e-4     # Coriolis parameter [m s-1]

arm.theta      = 301.4     # initial mixed-layer potential temperature [K]
arm.dtheta     = 0.4       # initial temperature jump at h [K]
arm.gammatheta = None      # free atmosphere potential temperature lapse rate [K m-1]
arm.advtheta   = 0.        # advection of heat [K s-1]
arm.beta       = 0.15      # entrainment ratio for virtual heat [-]
arm.wtheta     = None      # surface kinematic heat flux [K m s-1]

arm.q          = 15.3e-3   # initial mixed-layer specific humidity [kg kg-1]
arm.dq         = -0.2e-3   # initial specific humidity jump at h [kg kg-1]
arm.gammaq     = None      # free atmosphere specific humidity lapse rate [kg kg-1 m-1]
arm.advq       = 0.        # advection of moisture [kg kg-1 s-1]
arm.wq         = None      # surface kinematic moisture flux [kg kg-1 m s-1]

arm.CO2        = 422.      # initial mixed-layer CO2 [ppm]
arm.dCO2       = -44.      # initial CO2 jump at h [ppm]
arm.gammaCO2   = 0.        # free atmosphere CO2 lapse rate [ppm m-1]
arm.advCO2     = 0.        # advection of CO2 [ppm s-1]
arm.wCO2       = 0.        # surface kinematic CO2 flux [ppm m s-1]

arm.sw_wind    = False     # prognostic wind switch
arm.u          = 6.        # initial mixed-layer u-wind speed [m s-1]
arm.du         = 4.        # initial u-wind jump at h [m s-1]
arm.gammau     = 0.        # free atmosphere u-wind speed lapse rate [s-1]
arm.advu       = 0.        # advection of u-wind [m s-2]

arm.v          = -4.0      # initial mixed-layer u-wind speed [m s-1]
arm.dv         = 4.0       # initial u-wind jump at h [m s-1]
arm.gammav     = 0.        # free atmosphere v-wind speed lapse rate [s-1]
arm.advv       = 0.        # advection of v-wind [m s-2]

arm.sw_sl      = False     # surface layer switch
arm.ustar      = 0.3       # surface friction velocity [m s-1]
arm.z0m        = 0.02      # roughness length for momentum [m]
arm.z0h        = 0.002     # roughness length for scalars [m]

arm.sw_rad     = False     # radiation switch
arm.lat        = 51.97     # latitude [deg]
arm.lon        = -4.93     # longitude [deg]
arm.doy        = 268.      # day of the year [-]
arm.tstart     = 12.5      # time of the day [h UTC]
arm.cc         = 0.0       # cloud cover fraction [-]
arm.Q          = 400.      # net radiation [W m-2] 
arm.dFz        = 0.        # cloud top radiative divergence [W m-2] 

arm.sw_ls      = False     # land surface switch
arm.ls_type    = 'js'      # land-surface parameterization ('js' for Jarvis-Stewart or 'ags' for A-Gs)
arm.wg         = 0.21      # volumetric water content top soil layer [m3 m-3]
arm.w2         = 0.21      # volumetric water content deeper soil layer [m3 m-3]
arm.cveg       = 0.85      # vegetation fraction [-]
arm.Tsoil      = 285.      # temperature top soil layer [K]
arm.T2         = 286.      # temperature deeper soil layer [K]
arm.a          = 0.219     # Clapp and Hornberger retention curve parameter a
arm.b          = 4.90      # Clapp and Hornberger retention curve parameter b
arm.p          = 4.        # Clapp and Hornberger retention curve parameter c
arm.CGsat      = 3.56e-6   # saturated soil conductivity for heat

arm.wsat       = 0.472     # saturated volumetric water content ECMWF config [-]
arm.wfc        = 0.323     # volumetric water content field capacity [-]
arm.wwilt      = 0.171     # volumetric water content wilting point [-]

arm.C1sat      = 0.132     
arm.C2ref      = 1.8

arm.LAI        = 2.        # leaf area index [-]
arm.gD         = 0.0       # correction factor transpiration for VPD [-]
arm.rsmin      = 110.      # minimum resistance transpiration [s m-1]
arm.rssoilmin  = 50.       # minimun resistance soil evaporation [s m-1]
arm.alpha      = 0.25      # surface albedo [-]

arm.Ts         = 290.      # initial surface temperature [K]

arm.Wmax       = 0.0002    # thickness of water layer on wet vegetation [m]
arm.Wl         = 0.0000    # equivalent water layer depth for wet vegetation [m]

arm.Lambda     = 5.9       # thermal diffusivity skin layer [-]

arm.c3c4       = 'c3'      # Plant type ('c3' or 'c4')

arm.sw_cu      = True      # Cumulus parameterization switch
arm.dz_h       = 150.      # Transition layer thickness [m]

# Time dependent surface variables; linearly interpolated by the model
# Note the time offset, as the mixed-layer model starts one hour later than LES!
time   = np.array([0., 4, 6.5,  7.5,  10, 12.5, 14.5])-1
H      = np.array([-30.,  90., 140., 140., 100., -10.,  -10])
LE     = np.array([  5., 250., 450., 500., 420., 180.,    0])
rho    = arm.Ps / (287. * arm.theta * (1. + 0.61 * arm.q))
wtheta = H  / (rho*1005.)
wq     = LE / (rho*2.5e6)
time   = time*3600.

arm.timedep    = {'wtheta': (time, wtheta),
                  'wq':     (time, wq)}

# Binned height dependent lapse rates
z1         = np.array([0, 700, 5000])
gammatheta = np.array([3.4e-3, 5.7e-3])

z2         = np.array([0, 650, 1300, 5000])
gammaq     = np.array([-0.6e-6, -2e-6, -8.75e-6])

arm.heightdep  = {'gammatheta': (z1, gammatheta),
                  'gammaq':     (z2, gammaq)}


"""
Init and run the model
"""
# With cloud parameterisation
r1 = model(arm)
r1.run()

# Without cloud parameterisation
arm.sw_cu = False
r2 = model(arm)
r2.run()

"""
Plot output
"""
pl.close('all')

pl.figure()
pl.subplot(331)
pl.plot(r1.out.t, r1.out.h, label='ARM')
pl.plot(r2.out.t, r2.out.h, label='ARM no-clouds')
pl.xlabel('time [h]')
pl.ylabel('h [m]')
pl.legend()

pl.subplot(332)
pl.plot(r1.out.t, r1.out.theta)
pl.plot(r2.out.t, r2.out.theta)
pl.xlabel('time [h]')
pl.ylabel('theta [K]')

pl.subplot(333)
pl.plot(r1.out.t, r1.out.q*1000.)
pl.plot(r2.out.t, r2.out.q*1000.)
pl.xlabel('time [h]')
pl.ylabel('q [g kg-1]')

pl.subplot(334)
pl.plot(r1.out.t, r1.out.dtheta)
pl.plot(r2.out.t, r2.out.dtheta)
pl.xlabel('time [h]')
pl.ylabel('dtheta [K]')

pl.subplot(335)
pl.plot(r1.out.t, r1.out.dq*1000.)
pl.plot(r2.out.t, r2.out.dq*1000.)
pl.xlabel('time [h]')
pl.ylabel('dq [g kg-1]')

pl.subplot(336)
pl.plot(r1.out.t, r1.out.ac)
pl.plot(r2.out.t, r2.out.ac)
pl.xlabel('time [h]')
pl.ylabel('cloud core frac [-]')

pl.subplot(337)
pl.plot(r1.out.t, r1.out.wtheta, label='surface')
pl.plot(r1.out.t, r1.out.wthetae, label='entrainment')
pl.xlabel('time [h]')
pl.ylabel('wtheta [K m s-1]')
pl.legend()

pl.subplot(338)
pl.plot(r1.out.t, r1.out.wq*1000, label='surface')
pl.plot(r1.out.t, r1.out.wqe*1000, label='entrainment')
pl.xlabel('time [h]')
pl.ylabel('wq [g kg-1 m s-1]')
