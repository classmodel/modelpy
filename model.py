# 
# CLASS
# Copyright (c) 2010-2015 Meteorology and Air Quality section, Wageningen University and Research centre
# Copyright (c) 2011-2015 Jordi Vila-Guerau de Arellano
# Copyright (c) 2011-2015 Chiel van Heerwaarden
# Copyright (c) 2011-2015 Bart van Stratum
# Copyright (c) 2011-2015 Kees van den Dries
# 
# This file is part of CLASS
# 
# CLASS is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# CLASS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with CLASS.  If not, see <http://www.gnu.org/licenses/>.
# 
# Open issues:
# sensible heat flux always peaks too early in the morning, maybe other Jarvis functions?
#

import copy as cp
import numpy as np
import sys
#import ribtol

def esat(T):
    return 0.611e3 * np.exp(17.2694 * (T - 273.16) / (T - 35.86))

def qsat(T,p):
    return 0.622 * esat(T) / p

class model:
    def __init__(self, model_input):
        # initialize the different components of the model
        self.input = cp.deepcopy(model_input)
  
    def run(self):
        # initialize model variables
        self.init()
  
        # time integrate model 
        for self.t in range(self.tsteps):
          
            # time integrate components
            self.timestep()
  
        # delete unnecessary variables from memory
        self.exitmodel()
    
    def init(self):
        # assign variables from input data
        # initialize constants
        self.Lv         = 2.45e6                # heat of vaporization [J kg-1]
        self.cp         = 1005.                 # specific heat of dry air [J kg-1 K-1]
        self.rho        = 1.2                   # density of air [kg m-3]
        self.k          = 0.4                   # Von Karman constant [-]
        self.g          = 9.81                  # gravity acceleration [m s-2]
        self.Rd         = 287.                  # gas constant for dry air [J kg-1 K-1]
        self.Rv         = 461.5                 # gas constant for moist air [J kg-1 K-1]
        self.bolz       = 5.67e-8               # Bolzman constant [-]
        self.rhow       = 1000.                 # density of water [kg m-3]
        self.S0         = 1368.                 # solar constant [W m-2]

        # A-Gs constants and settings
        # Plant type:       -C3-     -C4-
        self.CO2comp298 =  [68.5,    4.3    ]   # CO2 compensation concentration [mg m-3]
        self.Q10CO2     =  [1.5,     1.5    ]   # function parameter to calculate CO2 compensation concentration [-]
        self.gm298      =  [7.0,     17.5   ]   # mesophyill conductance at 298 K [mm s-1]
        self.Ammax298   =  [2.2,     1.7    ]   # CO2 maximal primary productivity [mg m-2 s-1]
        self.Q10gm      =  [2.0,     2.0    ]   # function parameter to calculate mesophyll conductance [-]
        self.T1gm       =  [278.,    286.   ]   # reference temperature to calculate mesophyll conductance gm [K]
        self.T2gm       =  [301.,    309.   ]   # reference temperature to calculate mesophyll conductance gm [K]
        self.Q10Am      =  [2.0,     2.0    ]   # function parameter to calculate maximal primary profuctivity Ammax
        self.T1Am       =  [281.,    286.   ]   # reference temperature to calculate maximal primary profuctivity Ammax [K]
        self.T2Am       =  [311.,    311.   ]   # reference temperature to calculate maximal primary profuctivity Ammax [K]
        self.f0         =  [0.89,    0.85   ]   # maximum value Cfrac [-]
        self.ad         =  [0.07,    0.15   ]   # regression coefficient to calculate Cfrac [kPa-1]
        self.alpha0     =  [0.017,   0.014  ]   # initial low light conditions [mg J-1]
        self.Kx         =  [0.7,     0.7    ]   # extinction coefficient PAR [-]
        self.gmin       =  [0.25e-3, 0.25e-3]   # cuticular (minimum) conductance [mm s-1]

        self.mco2       =  44.;                 # molecular weight CO2 [g mol -1]
        self.mair       =  28.9;                # molecular weight air [g mol -1]
        self.nuco2q     =  1.6;                 # ratio molecular viscosity water to carbon dioxide

        self.Cw         =  0.0016;              # constant water stress correction (eq. 13 Jacobs et al. 2007) [-]
        self.wmax       =  0.55;                # upper reference value soil water [-]
        self.wmin       =  0.005;               # lower reference value soil water [-]
        self.R10        =  0.23;                # respiration at 10 C [mg CO2 m-2 s-1]
        self.E0         =  53.3e3;              # activation energy [53.3 kJ kmol-1]

        # Read switches
        self.sw_ml      = self.input.sw_ml      # mixed-layer model switch
        self.sw_shearwe = self.input.sw_shearwe # shear growth ABL switch
        self.sw_wind    = self.input.sw_wind    # prognostic wind switch
        self.sw_sl      = self.input.sw_sl      # surface layer switch
        self.sw_rad     = self.input.sw_rad     # radiation switch
        self.sw_ls      = self.input.sw_ls      # land surface switch
        self.ls_type    = self.input.ls_type    # land surface paramaterization (js or ags)
        self.sw_cu      = self.input.sw_cu      # cumulus parameterization switch
  
        # initialize mixed-layer
        self.h          = self.input.h          # initial ABL height [m]
        self.Ps         = self.input.Ps         # surface pressure [Pa]
        self.ws         = self.input.ws         # large scale vertical velocity [m s-1]
        self.fc         = self.input.fc         # coriolis parameter [s-1]
        self.we         = -1.                   # entrainment velocity [m s-1]
       
         # Temperature 
        self.theta      = self.input.theta      # initial mixed-layer potential temperature [K]
        self.dtheta     = self.input.dtheta     # initial temperature jump at h [K]
        self.gammatheta = self.input.gammatheta # free atmosphere potential temperature lapse rate [K m-1]
        self.advtheta   = self.input.advtheta   # advection of heat [K s-1]
        self.beta       = self.input.beta       # entrainment ratio for virtual heat [-]
        self.wtheta     = self.input.wtheta     # surface kinematic heat flux [K m s-1]
        self.wthetae    = None                  # entrainment kinematic heat flux [K m s-1]
 
        self.wstar      = 0.                    # convective velocity scale [m s-1]
 
        # 2m diagnostic variables 
        self.T2m        = None                  # 2m temperature [K]
        self.q2m        = None                  # 2m specific humidity [kg kg-1]
        self.e2m        = None                  # 2m vapor pressure [Pa]
        self.esat2m     = None                  # 2m saturated vapor pressure [Pa]
        self.u2m        = None                  # 2m u-wind [m s-1]
        self.v2m        = None                  # 2m v-wind [m s-1]
 
        # Surface variables 
        self.thetasurf  = self.input.theta      # surface potential temperature [K]
        self.thetavsurf = None                  # surface virtual potential temperature [K]
        self.qsurf      = None                  # surface specific humidity [g kg-1]

        # Mixed-layer top variables
        self.P_h        = None                  # Mixed-layer top pressure [pa]
        self.T_h        = None                  # Mixed-layer top absolute temperature [K]
        self.q2_h       = None                  # Mixed-layer top specific humidity variance [kg2 kg-2]
        self.RH_h       = None                  # Mixed-layer top relavtive humidity [-]
        self.dz_h       = None                  # Transition layer thickness [-]

        # Virtual temperatures and fluxes
        self.thetav     = None                  # initial mixed-layer potential temperature [K]
        self.dthetav    = None                  # initial virtual temperature jump at h [K]
        self.wthetav    = None                  # surface kinematic virtual heat flux [K m s-1]
        self.wthetave   = None                  # entrainment kinematic virtual heat flux [K m s-1]
       
        # Moisture 
        self.q          = self.input.q          # initial mixed-layer specific humidity [kg kg-1]
        self.dq         = self.input.dq         # initial specific humidity jump at h [kg kg-1]
        self.gammaq     = self.input.gammaq     # free atmosphere specific humidity lapse rate [kg kg-1 m-1]
        self.advq       = self.input.advq       # advection of moisture [kg kg-1 s-1]
        self.wq         = self.input.wq         # surface kinematic moisture flux [kg kg-1 m s-1]
        self.wqe        = None                  # entrainment moisture flux [kg kg-1 m s-1]
        self.wqM        = None                  # moisture cumulus mass flux [kg kg-1 m s-1]
  
        self.qsat       = None                  # mixed-layer saturated specific humidity [kg kg-1]
        self.esat       = None                  # mixed-layer saturated vapor pressure [Pa]
        self.e          = None                  # mixed-layer vapor pressure [Pa]
        self.qsatsurf   = None                  # surface saturated specific humidity [g kg-1]
        self.dqsatdT    = None                  # slope saturated specific humidity curve [g kg-1 K-1]
      
        # CO2
        fac = self.mair / (self.rho*self.mco2)  # Conversion factor mgC m-2 s-1 to ppm m s-1
        self.CO2        = self.input.CO2        # initial mixed-layer CO2 [ppm]
        self.dCO2       = self.input.dCO2       # initial CO2 jump at h [ppm]
        self.gammaCO2   = self.input.gammaCO2   # free atmosphere CO2 lapse rate [ppm m-1]
        self.advCO2     = self.input.advCO2     # advection of CO2 [ppm s-1]
        self.wCO2       = self.input.wCO2 * fac # surface kinematic CO2 flux [ppm m s-1]
        self.wCO2A      = 0                     # surface assimulation CO2 flux [ppm m s-1]
        self.wCO2R      = 0                     # surface respiration CO2 flux [ppm m s-1]
        self.wCO2e      = None                  # entrainment CO2 flux [ppm m s-1]
        self.wCO2M      = 0                     # CO2 mass flux [ppm m s-1]
       
        # Wind 
        self.u          = self.input.u          # initial mixed-layer u-wind speed [m s-1]
        self.du         = self.input.du         # initial u-wind jump at h [m s-1]
        self.gammau     = self.input.gammau     # free atmosphere u-wind speed lapse rate [s-1]
        self.advu       = self.input.advu       # advection of u-wind [m s-2]
        
        self.v          = self.input.v          # initial mixed-layer u-wind speed [m s-1]
        self.dv         = self.input.dv         # initial u-wind jump at h [m s-1]
        self.gammav     = self.input.gammav     # free atmosphere v-wind speed lapse rate [s-1]
        self.advv       = self.input.advv       # advection of v-wind [m s-2]
 
        # Tendencies 
        self.htend      = None                  # tendency of CBL [m s-1]
        self.thetatend  = None                  # tendency of mixed-layer potential temperature [K s-1]
        self.dthetatend = None                  # tendency of potential temperature jump at h [K s-1]
        self.qtend      = None                  # tendency of mixed-layer specific humidity [kg kg-1 s-1]
        self.dqtend     = None                  # tendency of specific humidity jump at h [kg kg-1 s-1]
        self.CO2tend    = None                  # tendency of CO2 humidity [ppm]
        self.dCO2tend   = None                  # tendency of CO2 jump at h [ppm s-1]
        self.utend      = None                  # tendency of u-wind [m s-1 s-1]
        self.dutend     = None                  # tendency of u-wind jump at h [m s-1 s-1]
        self.vtend      = None                  # tendency of v-wind [m s-1 s-1]
        self.dvtend     = None                  # tendency of v-wind jump at h [m s-1 s-1]
  
        # initialize surface layer
        self.ustar      = self.input.ustar      # surface friction velocity [m s-1]
        self.uw         = None                  # surface momentum flux in u-direction [m2 s-2]
        self.vw         = None                  # surface momentum flux in v-direction [m2 s-2]
        self.z0m        = self.input.z0m        # roughness length for momentum [m]
        self.z0h        = self.input.z0h        # roughness length for scalars [m]
        self.Cm         = 1e12                  # drag coefficient for momentum [-]
        self.Cs         = 1e12                  # drag coefficient for scalars [-]
        self.L          = None                  # Obukhov length [m]
        self.Rib        = None                  # bulk Richardson number [-]
        self.ra         = None                  # aerodynamic resistance [s m-1]
  
        # initialize radiation
        self.lat        = self.input.lat        # latitude [deg]
        self.lon        = self.input.lon        # longitude [deg]
        self.doy        = self.input.doy        # day of the year [-]
        self.tstart     = self.input.tstart     # time of the day [-]
        self.cc         = self.input.cc         # cloud cover fraction [-]
        self.Swin       = None                  # incoming short wave radiation [W m-2]
        self.Swout      = None                  # outgoing short wave radiation [W m-2]
        self.Lwin       = None                  # incoming long wave radiation [W m-2]
        self.Lwout      = None                  # outgoing long wave radiation [W m-2]
        self.Q          = self.input.Q          # net radiation [W m-2]
  
        # initialize land surface
        self.wg         = self.input.wg         # volumetric water content top soil layer [m3 m-3]
        self.w2         = self.input.w2         # volumetric water content deeper soil layer [m3 m-3]
        self.Tsoil      = self.input.Tsoil      # temperature top soil layer [K]
        self.T2         = self.input.T2         # temperature deeper soil layer [K]
                           
        self.a          = self.input.a          # Clapp and Hornberger retention curve parameter a [-]
        self.b          = self.input.b          # Clapp and Hornberger retention curve parameter b [-]
        self.p          = self.input.p          # Clapp and Hornberger retention curve parameter p [-]
        self.CGsat      = self.input.CGsat      # saturated soil conductivity for heat
                           
        self.wsat       = self.input.wsat       # saturated volumetric water content ECMWF config [-]
        self.wfc        = self.input.wfc        # volumetric water content field capacity [-]
        self.wwilt      = self.input.wwilt      # volumetric water content wilting point [-]
                           
        self.C1sat      = self.input.C1sat      
        self.C2ref      = self.input.C2ref      
        
        self.LAI        = self.input.LAI        # leaf area index [-]
        self.gD         = self.input.gD         # correction factor transpiration for VPD [-]
        self.rsmin      = self.input.rsmin      # minimum resistance transpiration [s m-1]
        self.rssoilmin  = self.input.rssoilmin  # minimum resistance soil evaporation [s m-1]
        self.alpha      = self.input.alpha      # surface albedo [-]
  
        self.rs         = 1.e6                  # resistance transpiration [s m-1]
        self.rssoil     = 1.e6                  # resistance soil [s m-1]
                           
        self.Ts         = self.input.Ts         # surface temperature [K]
                           
        self.cveg       = self.input.cveg       # vegetation fraction [-]
        self.Wmax       = self.input.Wmax       # thickness of water layer on wet vegetation [m]
        self.Wl         = self.input.Wl         # equivalent water layer depth for wet vegetation [m]
        self.cliq       = None                  # wet fraction [-]
                          
        self.Lambda     = self.input.Lambda     # thermal diffusivity skin layer [-]
  
        self.Tsoiltend  = None                  # soil temperature tendency [K s-1]
        self.wgtend     = None                  # soil moisture tendency [m3 m-3 s-1]
        self.Wltend     = None                  # equivalent liquid water tendency [m s-1]
  
        self.H          = None                  # sensible heat flux [W m-2]
        self.LE         = None                  # evapotranspiration [W m-2]
        self.LEliq      = None                  # open water evaporation [W m-2]
        self.LEveg      = None                  # transpiration [W m-2]
        self.LEsoil     = None                  # soil evaporation [W m-2]
        self.LEpot      = None                  # potential evaporation [W m-2]
        self.LEref      = None                  # reference evaporation using rs = rsmin / LAI [W m-2]
        self.G          = None                  # ground heat flux [W m-2]

        # initialize A-Gs surface scheme
        self.c3c4       = self.input.c3c4       # plant type ('c3' or 'c4')
 
        # initialize cumulus parameterization
        self.sw_cu      = self.input.sw_cu      # Cumulus parameterization switch
        self.dz_h       = self.input.dz_h       # Transition layer thickness [m]
        self.ac         = 0.                    # Cloud core fraction [-]
        self.M          = 0.                    # Cloud core mass flux [m s-1] 
        self.wqM        = 0.                    # Cloud core moisture flux [kg kg-1 m s-1] 
  
        # initialize time variables
        self.tsteps = int(np.floor(self.input.runtime / self.input.dt))
        self.dt     = self.input.dt
        self.t      = 0
  
        # initialize output
        self.out = model_output(self.tsteps)
 
        self.statistics()
  
        # calculate initial diagnostic variables
        if(self.sw_rad):
            self.run_radiation()
 
        if(self.sw_sl):
            for i in range(10): 
                self.run_surface_layer()
  
        if(self.sw_ls):
            self.run_land_surface()
       
        if(self.sw_cu):
            self.run_mixed_layer()
            self.run_cumulus()
        
        if(self.sw_ml):
            self.run_mixed_layer()

    def timestep(self):
        self.statistics()

        # run radiation model
        if(self.sw_rad):
            self.run_radiation()
  
        # run surface layer model
        if(self.sw_sl):
            self.run_surface_layer()
        
        # run land surface model
        if(self.sw_ls):
            self.run_land_surface()
 
        # run cumulus parameterization
        if(self.sw_cu):
            self.run_cumulus()
   
        # run mixed-layer model
        if(self.sw_ml):
            self.run_mixed_layer()
 
        # store output before time integration
        self.store()
  
        # time integrate land surface model
        if(self.sw_ls):
            self.integrate_land_surface()
  
        # time integrate mixed-layer model
        if(self.sw_ml):
            self.integrate_mixed_layer()
  
    def statistics(self):
        # Calculate virtual temperatures 
        self.thetav   = self.theta  + 0.61 * self.theta * self.q
        self.wthetav  = self.wtheta + 0.61 * self.theta * self.wq
        self.dthetav  = (self.theta + self.dtheta) * (1. + 0.61 * (self.q + self.dq)) - self.theta * (1. + 0.61 * self.q)

        # Mixed-layer top properties
        self.P_h    = self.Ps - self.rho * self.g * self.h
        self.T_h    = self.theta - self.g/self.cp * self.h

        self.P_h    = self.Ps / np.exp((self.g * self.h)/(self.Rd * self.theta))
        self.T_h    = self.theta / (self.Ps / self.P_h)**(self.Rd/self.cp)

        self.RH_h   = self.q / qsat(self.T_h, self.P_h)

    def run_cumulus(self):
        # Calculate mixed-layer top relative humidity variance (Neggers et. al 2006/7)
        wqe         = -self.we * self.dq
        self.q2_h   = -(wqe+self.wqM) * self.dq * self.h / (self.dz_h * self.wstar)

        # calculate cloud core fraction (ac), mass flux (M) and moisture flux (wqM)
        self.ac     = max(0., 0.5 + (0.36 * np.arctan(1.55 * ((self.q - qsat(self.T_h, self.P_h)) / self.q2_h**0.5))))
        self.M      = self.ac * self.wstar
        self.wqM    = self.M * self.q2_h**0.5

    def run_mixed_layer(self):
        if(not self.sw_sl):
            # decompose ustar along the wind components
            self.uw = - np.sign(self.u) * (self.ustar ** 4. / (self.v ** 2. / self.u ** 2. + 1.)) ** (0.5)
            self.vw = - np.sign(self.v) * (self.ustar ** 4. / (self.u ** 2. / self.v ** 2. + 1.)) ** (0.5)
       
        # calculate convective velocity scale w* 
        if(self.wthetav > 0.):
          self.wstar = ((self.g * self.h * self.wthetav) / self.thetav)**(1./3.)
        else:
          self.wstar  = 1e-6;
      
        # Virtual heat entrainment flux 
        self.wthetave    = -self.beta * self.wthetav 
        
        # compute mixed-layer tendencies
        if(self.sw_shearwe):
            self.we    = (-self.wthetave + 5. * self.ustar ** 3. * self.thetav / (self.g * self.h)) / self.dthetav
        else:
            self.we    = -self.wthetave / self.dthetav

        # Don't allow boundary layer shrinking if wtheta < 0 
        if(self.we < 0):
            self.we = 0.
 
        # Calculate entrainment fluxes
        self.wthetae     = -self.we * self.dtheta
        self.wqe         = -self.we * self.dq
        self.wCO2e       = -self.we * self.dCO2
  
        self.htend       = self.we + self.ws - self.M
       
        self.thetatend   = (self.wtheta - self.wthetae           ) / self.h + self.advtheta 
        self.qtend       = (self.wq     - self.wqe     - self.wqM) / self.h + self.advq
        self.CO2tend     = (self.wCO2   - self.wCO2e             ) / self.h + self.advCO2
        
        self.dthetatend  = self.gammatheta * (self.we - self.M) - self.thetatend
        self.dqtend      = self.gammaq     * (self.we - self.M) - self.qtend
        self.dCO2tend    = self.gammaCO2   * (self.we - self.M) - self.CO2tend
     
        # assume u + du = ug, so ug - u = du
        if(self.sw_wind):
            self.utend       = -self.fc * self.dv + (self.uw + self.we * self.du)  / self.h + self.advu
            self.vtend       =  self.fc * self.du + (self.vw + self.we * self.dv)  / self.h + self.advv
  
            self.dutend      = self.gammau * (self.we - self.M) - self.utend
            self.dvtend      = self.gammav * (self.we - self.M) - self.vtend
   
    def integrate_mixed_layer(self):
        # set values previous time step
        h0      = self.h
        
        theta0  = self.theta
        dtheta0 = self.dtheta
        q0      = self.q
        dq0     = self.dq
        CO20    = self.CO2
        dCO20   = self.dCO2
        
        u0      = self.u
        du0     = self.du
        v0      = self.v
        dv0     = self.dv
  
        # integrate mixed-layer equations
        self.h        = h0      + self.dt * self.htend
  
        self.theta    = theta0  + self.dt * self.thetatend
        self.dtheta   = dtheta0 + self.dt * self.dthetatend
        self.q        = q0      + self.dt * self.qtend
        self.dq       = dq0     + self.dt * self.dqtend
        self.CO2      = CO20    + self.dt * self.CO2tend
        self.dCO2     = dCO20   + self.dt * self.dCO2tend
  
        if(self.sw_wind):
            self.u        = u0      + self.dt * self.utend
            self.du       = du0     + self.dt * self.dutend
            self.v        = v0      + self.dt * self.vtend
            self.dv       = dv0     + self.dt * self.dvtend
 
    def run_radiation(self):
        sda    = 0.409 * np.cos(2. * np.pi * (self.doy - 173.) / 365.)
        sinlea = np.sin(2. * np.pi * self.lat / 360.) * np.sin(sda) - np.cos(2. * np.pi * self.lat / 360.) * np.cos(sda) * np.cos(2. * np.pi * (self.t * self.dt + self.tstart * 3600.) / 86400. - 2. * np.pi * self.lon / 360.)
        sinlea = max(sinlea, 0.0001)
        
        Ta  = self.theta * ((self.Ps - 0.1 * self.h * self.rho * self.g) / self.Ps ) ** (self.Rd / self.cp)
  
        Tr  = (0.6 + 0.2 * sinlea) * (1. - 0.4 * self.cc)
  
        self.Swin  = self.S0 * Tr * sinlea
        self.Swout = self.alpha * self.S0 * Tr * sinlea
        self.Lwin  = 0.8 * self.bolz * Ta ** 4.
        self.Lwout = self.bolz * self.Ts ** 4.
          
        self.Q     = self.Swin - self.Swout + self.Lwin - self.Lwout
  
    def run_surface_layer(self):
        ueff           = max(0.01, np.sqrt(self.u**2. + self.v**2. + self.wstar**2.))
        self.thetasurf = self.theta + self.wtheta / (self.Cs * ueff)
        qsatsurf       = qsat(self.thetasurf, self.Ps)
        cq             = (1. + self.Cs * ueff * self.rs) ** -1.
        self.qsurf     = (1. - cq) * self.q + cq * qsatsurf

        self.thetavsurf = self.thetasurf * (1. + 0.61 * self.qsurf)
  
        zsl       = 0.1 * self.h
        self.Rib  = self.g / self.thetav * zsl * (self.thetav - self.thetavsurf) / ueff**2.
        self.Rib  = min(self.Rib, 0.2)

        self.L     = self.ribtol(self.Rib, zsl, self.z0m, self.z0h)  # Slow python iteration
        #self.L    = ribtol.ribtol(self.Rib, zsl, self.z0m, self.z0h) # Fast C++ iteration
 
        self.Cm   = self.k**2. / (np.log(zsl / self.z0m) - self.psim(zsl / self.L) + self.psim(self.z0m / self.L)) ** 2.
        self.Cs   = self.k**2. / (np.log(zsl / self.z0m) - self.psim(zsl / self.L) + self.psim(self.z0m / self.L)) / (np.log(zsl / self.z0h) - self.psih(zsl / self.L) + self.psih(self.z0h / self.L))
  
        self.ustar = np.sqrt(self.Cm) * ueff
        self.uw    = - self.Cm * ueff * self.u
        self.vw    = - self.Cm * ueff * self.v
 
        # diagnostic meteorological variables
        self.T2m    = self.thetasurf - self.wtheta / self.ustar / self.k * (np.log(2. / self.z0h) - self.psih(2. / self.L) + self.psih(self.z0h / self.L))
        self.q2m    = self.qsurf     - self.wq     / self.ustar / self.k * (np.log(2. / self.z0h) - self.psih(2. / self.L) + self.psih(self.z0h / self.L))
        self.u2m    =                - self.uw     / self.ustar / self.k * (np.log(2. / self.z0m) - self.psim(2. / self.L) + self.psim(self.z0m / self.L))
        self.v2m    =                - self.vw     / self.ustar / self.k * (np.log(2. / self.z0m) - self.psim(2. / self.L) + self.psim(self.z0m / self.L))
        self.esat2m = 0.611e3 * np.exp(17.2694 * (self.T2m - 273.16) / (self.T2m - 35.86))
        self.e2m    = self.q2m * self.Ps / 0.622
     
    def ribtol(self, Rib, zsl, z0m, z0h): 
        if(Rib > 0.):
            L    = 1.
            L0   = 2.
        else:
            L  = -1.
            L0 = -2.
        
        while (abs(L - L0) > 0.001):
            L0      = L
            fx      = Rib - zsl / L * (np.log(zsl / z0h) - self.psih(zsl / L) + self.psih(z0h / L)) / (np.log(zsl / z0m) - self.psim(zsl / L) + self.psim(z0m / L))**2.
            Lstart  = L - 0.001*L
            Lend    = L + 0.001*L
            fxdif   = ( (- zsl / Lstart * (np.log(zsl / z0h) - self.psih(zsl / Lstart) + self.psih(z0h / Lstart)) / \
                                          (np.log(zsl / z0m) - self.psim(zsl / Lstart) + self.psim(z0m / Lstart))**2.) \
                      - (-zsl /  Lend   * (np.log(zsl / z0h) - self.psih(zsl / Lend  ) + self.psih(z0h / Lend  )) / \
                                          (np.log(zsl / z0m) - self.psim(zsl / Lend  ) + self.psim(z0m / Lend  ))**2.) ) / (Lstart - Lend)
            L       = L - fx / fxdif

            if(abs(L) > 1e15):
                break

        return L
      
    def psim(self, zeta):
        if(zeta <= 0):
            x     = (1. - 16. * zeta)**(0.25)
            psim  = 3.14159265 / 2. - 2. * np.arctan(x) + np.log((1. + x)**2. * (1. + x**2.) / 8.)
            #x     = (1. + 3.6 * abs(zeta) ** (2./3.)) ** (-0.5)
            #psim = 3. * np.log( (1. + 1. / x) / 2.)
        else:
            psim  = -2./3. * (zeta - 5./0.35) * np.exp(-0.35 * zeta) - zeta - (10./3.) / 0.35
        return psim
      
    def psih(self, zeta):
        if(zeta <= 0):
            x     = (1. - 16. * zeta)**(0.25)
            psih  = 2. * np.log( (1. + x*x) / 2.)
            #x     = (1. + 7.9 * abs(zeta) ** (2./3.)) ** (-0.5)
            #psih  = 3. * np.log( (1. + 1. / x) / 2.)
        else:
            psih  = -2./3. * (zeta - 5./0.35) * np.exp(-0.35 * zeta) - (1. + (2./3.) * zeta) ** (1.5) - (10./3.) / 0.35 + 1.
        return psih
 
    def jarvis_stewart(self):
        # calculate surface resistances using Jarvis-Stewart model
        if(self.sw_rad):
            f1 = 1. / min(1.,((0.004 * self.Swin + 0.05) / (0.81 * (0.004 * self.Swin + 1.))))
        else:
            f1 = 1.
  
        if(self.w2 > self.wwilt):# and self.w2 <= self.wfc):
            f2 = (self.wfc - self.wwilt) / (self.w2 - self.wwilt)
        else:
            f2 = 1.e8
 
        # Limit f2 in case w2 > wfc, where f2 < 1
        f2 = max(f2, 1.);
 
        f3 = 1. / np.exp(- self.gD * (self.esat - self.e) / 100.)
        f4 = 1./ (1. - 0.0016 * (298.0-self.theta)**2.)
  
        self.rs = self.rsmin / self.LAI * f1 * f2 * f3 * f4

    def factorial(self,k):
        factorial = 1
        for n in range(2,k+1):
            factorial = factorial * float(n)
        return factorial;

    def E1(self,x):
        E1sum = 0
        for k in range(1,100):
            E1sum += pow((-1.),(k + 0.0)) * pow(x,(k + 0.0)) / ((k + 0.0) * self.factorial(k))
        return -0.57721566490153286060 - np.log(x) - E1sum
 
    def ags(self):
        # Select index for plant type
        if(self.c3c4 == 'c3'):
            c = 0
        elif(self.c3c4 == 'c4'):
            c = 1
        else:
            sys.exit('option \"%s\" for \"c3c4\" invalid'%self.c3c4)

        # calculate CO2 compensation concentration
        CO2comp       = self.CO2comp298[c] * self.rho * pow(self.Q10CO2[c],(0.1 * (self.thetasurf - 298.)))  

        # calculate mesophyll conductance
        gm            = self.gm298[c] *  pow(self.Q10gm[c],(0.1 * (self.thetasurf-298.))) \
                          / ( (1. + np.exp(0.3 * (self.T1gm[c] - self.thetasurf))) * (1. + np.exp(0.3 * (self.thetasurf - self.T2gm[c]))))
        gm            = gm / 1000. # conversion from mm s-1 to m s-1
  
        # calculate CO2 concentration inside the leaf (ci)
        fmin0         = self.gmin[c] / self.nuco2q - 1. / 9. * gm
        fmin          = -fmin0 + pow((pow(fmin0,2.) + 4 * self.gmin[c]/self.nuco2q * gm),0.5) / (2. * gm)
  
        Ds            = (esat(self.Ts) - self.e) / 1000. # kPa
        D0            = (self.f0[c] - fmin) / self.ad[c]
  
        cfrac         = self.f0[c] * (1. - (Ds / D0)) + fmin * (Ds / D0)
        co2abs        = self.CO2 * (self.mco2 / self.mair) * self.rho # conversion mumol mol-1 (ppm) to mgCO2 m3
        ci            = cfrac * (co2abs - CO2comp) + CO2comp
  
        # calculate maximal gross primary production in high light conditions (Ag)
        Ammax         = self.Ammax298[c] *  pow(self.Q10Am[c],(0.1 * (self.thetasurf - 298.))) / ( (1. + np.exp(0.3 * (self.T1Am[c] - self.thetasurf))) * (1. + np.exp(0.3 * (self.thetasurf - self.T2Am[c]))))
  
        # calculate effect of soil moisture stress on gross assimilation rate
        betaw         = max(1e-3, min(1.,(self.wg - self.wwilt)/(self.wfc - self.wwilt)))
  
        # calculate stress function
        fstr          = betaw;
  
        # calculate gross assimilation rate (Am)
        Am           = Ammax * (1. - np.exp(-(gm * (ci - CO2comp) / Ammax)))
        Rdark        = (1. / 9.) * Am
        PAR          = 0.5 * max(1e-1,self.Swin * self.cveg)
  
        # calculate  light use efficiency
        alphac       = self.alpha0[c] * (co2abs - CO2comp) / (co2abs + 2. * CO2comp)
  
        # calculate gross primary productivity
        Ag           = (Am + Rdark) * (1 - np.exp(alphac * PAR / (Am + Rdark)))
  
        # 1.- calculate upscaling from leaf to canopy: net flow CO2 into the plant (An)
        y            =  alphac * self.Kx[c] * PAR / (Am + Rdark)
        An           = (Am + Rdark) * (1. - 1. / (self.Kx[c] * self.LAI) * (self.E1(y * np.exp(-self.Kx[c] * self.LAI)) - self.E1(y)))
  
        # 2.- calculate upscaling from leaf to canopy: CO2 conductance at canopy level
        a1           = 1. / (1. - self.f0[c])
        Dstar        = D0 / (a1 * (self.f0[c] - fmin))
  
        gcco2        = self.LAI * (self.gmin[c] / self.nuco2q + a1 * fstr * An / ((co2abs - CO2comp) * (1. + Ds / Dstar)))
  
        # calculate surface resistance for moisture and carbon dioxide
        self.rs      = 1. / (1.6 * gcco2)
        rsCO2        = 1. / gcco2
  
        # calculate net flux of CO2 into the plant (An)
        An           = -(co2abs - ci) / (self.ra + rsCO2)
  
        # CO2 soil surface flux
        fw           = self.Cw * self.wmax / (self.wg + self.wmin)
        Resp         = self.R10 * (1. - fw) * np.exp(self.E0 / (283.15 * 8.314) * (1. - 283.15 / (self.Tsoil)))
  
        # CO2 flux
        self.wCO2A   = An   * (self.mair / (self.rho * self.mco2))
        self.wCO2R   = Resp * (self.mair / (self.rho * self.mco2))
        self.wCO2    = self.wCO2A + self.wCO2R
 
    def run_land_surface(self):
        # compute ra
        ueff = np.sqrt(self.u ** 2. + self.v ** 2. + self.wstar**2.)

        if(self.sw_sl):
          self.ra = (self.Cs * ueff)**-1.
        else:
          self.ra = ueff / max(1.e-3, self.ustar)**2.

        # first calculate essential thermodynamic variables
        self.esat    = esat(self.theta)
        self.qsat    = qsat(self.theta, self.Ps)
        desatdT      = self.esat * (17.2694 / (self.theta - 35.86) - 17.2694 * (self.theta - 273.16) / (self.theta - 35.86)**2.)
        self.dqsatdT = 0.622 * desatdT / self.Ps
        self.e       = self.q * self.Ps / 0.622

        if(self.ls_type == 'js'): 
            self.jarvis_stewart() 
        elif(self.ls_type == 'ags'):
            self.ags()
        else:
            sys.exit('option \"%s\" for \"ls_type\" invalid'%self.ls_type)

        # recompute f2 using wg instead of w2
        if(self.wg > self.wwilt):# and self.w2 <= self.wfc):
          f2          = (self.wfc - self.wwilt) / (self.wg - self.wwilt)
        else:
          f2        = 1.e8
        self.rssoil = self.rssoilmin * f2 
 
        Wlmx = self.LAI * self.Wmax
        self.cliq = min(1., self.Wl / Wlmx) 
     
        # calculate skin temperature implictly
        self.Ts   = (self.Q  + self.rho * self.cp / self.ra * self.theta \
            + self.cveg * (1. - self.cliq) * self.rho * self.Lv / (self.ra + self.rs    ) * (self.dqsatdT * self.theta - self.qsat + self.q) \
            + (1. - self.cveg)             * self.rho * self.Lv / (self.ra + self.rssoil) * (self.dqsatdT * self.theta - self.qsat + self.q) \
            + self.cveg * self.cliq        * self.rho * self.Lv /  self.ra                * (self.dqsatdT * self.theta - self.qsat + self.q) + self.Lambda * self.Tsoil) \
            / (self.rho * self.cp / self.ra + self.cveg * (1. - self.cliq) * self.rho * self.Lv / (self.ra + self.rs) * self.dqsatdT \
            + (1. - self.cveg) * self.rho * self.Lv / (self.ra + self.rssoil) * self.dqsatdT + self.cveg * self.cliq * self.rho * self.Lv / self.ra * self.dqsatdT + self.Lambda)

        esatsurf      = esat(self.Ts)
        self.qsatsurf = qsat(self.Ts, self.Ps)

        self.LEveg  = (1. - self.cliq) * self.cveg * self.rho * self.Lv / (self.ra + self.rs) * (self.dqsatdT * (self.Ts - self.theta) + self.qsat - self.q)
        self.LEliq  = self.cliq * self.cveg * self.rho * self.Lv / self.ra * (self.dqsatdT * (self.Ts - self.theta) + self.qsat - self.q)
        self.LEsoil = (1. - self.cveg) * self.rho * self.Lv / (self.ra + self.rssoil) * (self.dqsatdT * (self.Ts - self.theta) + self.qsat - self.q)
  
        self.Wltend      = - self.LEliq / (self.rhow * self.Lv)
  
        self.LE     = self.LEsoil + self.LEveg + self.LEliq
        self.H      = self.rho * self.cp / self.ra * (self.Ts - self.theta)
        self.G      = self.Lambda * (self.Ts - self.Tsoil)
        self.LEpot  = (self.dqsatdT * (self.Q - self.G) + self.rho * self.cp / self.ra * (self.qsat - self.q)) / (self.dqsatdT + self.cp / self.Lv)
        self.LEref  = (self.dqsatdT * (self.Q - self.G) + self.rho * self.cp / self.ra * (self.qsat - self.q)) / (self.dqsatdT + self.cp / self.Lv * (1. + self.rsmin / self.LAI / self.ra))
        
        CG          = self.CGsat * (self.wsat / self.w2)**(self.b / (2. * np.log(10.)))
  
        self.Tsoiltend   = CG * self.G - 2. * np.pi / 86400. * (self.Tsoil - self.T2)
   
        d1          = 0.1
        C1          = self.C1sat * (self.wsat / self.wg) ** (self.b / 2. + 1.)
        C2          = self.C2ref * (self.w2 / (self.wsat - self.w2) )
        wgeq        = self.w2 - self.wsat * self.a * ( (self.w2 / self.wsat) ** self.p * (1. - (self.w2 / self.wsat) ** (8. * self.p)) )
        self.wgtend = - C1 / (self.rhow * d1) * self.LEsoil / self.Lv - C2 / 86400. * (self.wg - wgeq)
  
        # calculate kinematic heat fluxes
        self.wtheta   = self.H  / (self.rho * self.cp)
        self.wq       = self.LE / (self.rho * self.Lv)
  
    def integrate_land_surface(self):
        # integrate soil equations
        Tsoil0        = self.Tsoil
        wg0           = self.wg
        Wl0           = self.Wl
  
        self.Tsoil    = Tsoil0  + self.dt * self.Tsoiltend
        self.wg       = wg0     + self.dt * self.wgtend
        self.Wl       = Wl0     + self.dt * self.Wltend
  
    # store model output
    def store(self):
        t                      = self.t
        self.out.t[t]          = t * self.dt / 3600. + self.tstart
        self.out.h[t]          = self.h
        
        self.out.theta[t]      = self.theta
        self.out.thetav[t]     = self.thetav
        self.out.dtheta[t]     = self.dtheta
        self.out.dthetav[t]    = self.dthetav
        self.out.wtheta[t]     = self.wtheta
        self.out.wthetav[t]    = self.wthetav
        self.out.wthetae[t]    = self.wthetae
        self.out.wthetave[t]   = self.wthetave
        
        self.out.q[t]          = self.q
        self.out.dq[t]         = self.dq
        self.out.wq[t]         = self.wq
        self.out.wqe[t]        = self.wqe
        self.out.wqM[t]        = self.wqM
      
        self.out.qsat[t]       = self.qsat
        self.out.e[t]          = self.e
        self.out.esat[t]       = self.esat
      
        fac = (self.rho*self.mco2)/self.mair
        self.out.CO2[t]        = self.CO2
        self.out.dCO2[t]       = self.dCO2
        self.out.wCO2[t]       = self.wCO2  * fac
        self.out.wCO2e[t]      = self.wCO2e * fac
        self.out.wCO2R[t]      = self.wCO2R * fac
        self.out.wCO2A[t]      = self.wCO2A * fac

        self.out.u[t]          = self.u
        self.out.du[t]         = self.du
        self.out.uw[t]         = self.uw
        
        self.out.v[t]          = self.v
        self.out.dv[t]         = self.dv
        self.out.vw[t]         = self.vw
        
        self.out.T2m[t]        = self.T2m
        self.out.q2m[t]        = self.q2m
        self.out.u2m[t]        = self.u2m
        self.out.v2m[t]        = self.v2m
        self.out.e2m[t]        = self.e2m
        self.out.esat2m[t]     = self.esat2m
        
        self.out.thetasurf[t]  = self.thetasurf
        self.out.thetavsurf[t] = self.thetavsurf
        self.out.qsurf[t]      = self.qsurf
        self.out.ustar[t]      = self.ustar
        self.out.Cm[t]         = self.Cm
        self.out.Cs[t]         = self.Cs
        self.out.L[t]          = self.L
        self.out.Rib[t]        = self.Rib
  
        self.out.Swin[t]       = self.Swin
        self.out.Swout[t]      = self.Swout
        self.out.Lwin[t]       = self.Lwin
        self.out.Lwout[t]      = self.Lwout
        self.out.Q[t]          = self.Q
  
        self.out.ra[t]         = self.ra
        self.out.rs[t]         = self.rs
        self.out.H[t]          = self.H
        self.out.LE[t]         = self.LE
        self.out.LEliq[t]      = self.LEliq
        self.out.LEveg[t]      = self.LEveg
        self.out.LEsoil[t]     = self.LEsoil
        self.out.LEpot[t]      = self.LEpot
        self.out.LEref[t]      = self.LEref
        self.out.G[t]          = self.G

        self.out.ac[t]         = self.ac
        self.out.M[t]          = self.M
  
    # delete class variables to facilitate analysis in ipython
    def exitmodel(self):
        del(self.Lv)
        del(self.cp)
        del(self.rho)
        del(self.k)
        del(self.g)
        del(self.Rd)
        del(self.Rv)
        del(self.bolz)
        del(self.S0)
        del(self.rhow)
  
        del(self.t)
        del(self.dt)
        del(self.tsteps)
         
        del(self.h)          
        del(self.Ps)        
        del(self.fc)        
        del(self.ws)
        del(self.we)
        
        del(self.theta)
        del(self.dtheta)
        del(self.gammatheta)
        del(self.advtheta)
        del(self.beta)
        del(self.wtheta)
    
        del(self.T2m)
        del(self.q2m)
        del(self.e2m)
        del(self.esat2m)
        del(self.u2m)
        del(self.v2m)
        
        del(self.thetasurf)
        del(self.qsatsurf)
        del(self.thetav)
        del(self.dthetav)
        del(self.thetavsurf)
        del(self.qsurf)
        del(self.wthetav)
        
        del(self.q)
        del(self.qsat)
        del(self.dqsatdT)
        del(self.e)
        del(self.esat)
        del(self.dq)
        del(self.gammaq)
        del(self.advq)
        del(self.wq)
        
        del(self.u)
        del(self.du)
        del(self.gammau)
        del(self.advu)
        
        del(self.v)
        del(self.dv)
        del(self.gammav)
        del(self.advv)
  
        del(self.htend)
        del(self.thetatend)
        del(self.dthetatend)
        del(self.qtend)
        del(self.dqtend)
        del(self.utend)
        del(self.dutend)
        del(self.vtend)
        del(self.dvtend)
     
        del(self.Tsoiltend) 
        del(self.wgtend)  
        del(self.Wltend) 
  
        del(self.ustar)
        del(self.uw)
        del(self.vw)
        del(self.z0m)
        del(self.z0h)        
        del(self.Cm)         
        del(self.Cs)
        del(self.L)
        del(self.Rib)
        del(self.ra)
  
        del(self.lat)
        del(self.lon)
        del(self.doy)
        del(self.tstart)
   
        del(self.Swin)
        del(self.Swout)
        del(self.Lwin)
        del(self.Lwout)
        del(self.cc)
  
        del(self.wg)
        del(self.w2)
        del(self.cveg)
        del(self.cliq)
        del(self.Tsoil)
        del(self.T2)
        del(self.a)
        del(self.b)
        del(self.p)
        del(self.CGsat)
  
        del(self.wsat)
        del(self.wfc)
        del(self.wwilt)
  
        del(self.C1sat)
        del(self.C2ref)
  
        del(self.LAI)
        del(self.rs)
        del(self.rssoil)
        del(self.rsmin)
        del(self.rssoilmin)
        del(self.alpha)
        del(self.gD)
  
        del(self.Ts)
  
        del(self.Wmax)
        del(self.Wl)
  
        del(self.Lambda)
        
        del(self.Q)
        del(self.H)
        del(self.LE)
        del(self.LEliq)
        del(self.LEveg)
        del(self.LEsoil)
        del(self.LEpot)
        del(self.LEref)
        del(self.G)
  
        del(self.sw_ls)
        del(self.sw_rad)
        del(self.sw_sl)
        del(self.sw_wind)
        del(self.sw_shearwe)

# class for storing mixed-layer model output data
class model_output:
    def __init__(self, tsteps):
        self.t          = np.zeros(tsteps)    # time [s]

        # mixed-layer variables
        self.h          = np.zeros(tsteps)    # ABL height [m]
        
        self.theta      = np.zeros(tsteps)    # initial mixed-layer potential temperature [K]
        self.thetav     = np.zeros(tsteps)    # initial mixed-layer virtual potential temperature [K]
        self.dtheta     = np.zeros(tsteps)    # initial potential temperature jump at h [K]
        self.dthetav    = np.zeros(tsteps)    # initial virtual potential temperature jump at h [K]
        self.wtheta     = np.zeros(tsteps)    # surface kinematic heat flux [K m s-1]
        self.wthetav    = np.zeros(tsteps)    # surface kinematic virtual heat flux [K m s-1]
        self.wthetae    = np.zeros(tsteps)    # entrainment kinematic heat flux [K m s-1]
        self.wthetave   = np.zeros(tsteps)    # entrainment kinematic virtual heat flux [K m s-1]
        
        self.q          = np.zeros(tsteps)    # mixed-layer specific humidity [kg kg-1]
        self.dq         = np.zeros(tsteps)    # initial specific humidity jump at h [kg kg-1]
        self.wq         = np.zeros(tsteps)    # surface kinematic moisture flux [kg kg-1 m s-1]
        self.wqe        = np.zeros(tsteps)    # entrainment kinematic moisture flux [kg kg-1 m s-1]
        self.wqM        = np.zeros(tsteps)    # cumulus mass-flux kinematic moisture flux [kg kg-1 m s-1]

        self.qsat       = np.zeros(tsteps)    # mixed-layer saturated specific humidity [kg kg-1]
        self.e          = np.zeros(tsteps)    # mixed-layer vapor pressure [Pa]
        self.esat       = np.zeros(tsteps)    # mixed-layer saturated vapor pressure [Pa]

        self.CO2        = np.zeros(tsteps)    # mixed-layer CO2 [ppm]
        self.dCO2       = np.zeros(tsteps)    # initial CO2 jump at h [ppm]
        self.wCO2       = np.zeros(tsteps)    # surface total CO2 flux [mgC m-2 s-1]
        self.wCO2A      = np.zeros(tsteps)    # surface assimilation CO2 flux [mgC m-2 s-1]
        self.wCO2R      = np.zeros(tsteps)    # surface respiration CO2 flux [mgC m-2 s-1]
        self.wCO2e      = np.zeros(tsteps)    # entrainment CO2 flux [mgC m-2 s-1]
        self.wCO2M      = np.zeros(tsteps)    # CO2 mass flux [mgC m-2 s-1]
        
        self.u          = np.zeros(tsteps)    # initial mixed-layer u-wind speed [m s-1]
        self.du         = np.zeros(tsteps)    # initial u-wind jump at h [m s-1]
        self.uw         = np.zeros(tsteps)    # surface momentum flux u [m2 s-2]
        
        self.v          = np.zeros(tsteps)    # initial mixed-layer u-wind speed [m s-1]
        self.dv         = np.zeros(tsteps)    # initial u-wind jump at h [m s-1]
        self.vw         = np.zeros(tsteps)    # surface momentum flux v [m2 s-2]

        # diagnostic meteorological variables
        self.T2m        = np.zeros(tsteps)    # 2m temperature [K]   
        self.q2m        = np.zeros(tsteps)    # 2m specific humidity [kg kg-1]
        self.u2m        = np.zeros(tsteps)    # 2m u-wind [m s-1]    
        self.v2m        = np.zeros(tsteps)    # 2m v-wind [m s-1]    
        self.e2m        = np.zeros(tsteps)    # 2m vapor pressure [Pa]
        self.esat2m     = np.zeros(tsteps)    # 2m saturated vapor pressure [Pa]

        # surface-layer variables
        self.thetasurf  = np.zeros(tsteps)    # surface potential temperature [K]
        self.thetavsurf = np.zeros(tsteps)    # surface virtual potential temperature [K]
        self.qsurf      = np.zeros(tsteps)    # surface specific humidity [kg kg-1]
        self.ustar      = np.zeros(tsteps)    # surface friction velocity [m s-1]
        self.z0m        = np.zeros(tsteps)    # roughness length for momentum [m]
        self.z0h        = np.zeros(tsteps)    # roughness length for scalars [m]
        self.Cm         = np.zeros(tsteps)    # drag coefficient for momentum []
        self.Cs         = np.zeros(tsteps)    # drag coefficient for scalars []
        self.L          = np.zeros(tsteps)    # Obukhov length [m]
        self.Rib        = np.zeros(tsteps)    # bulk Richardson number [-]

        # radiation variables
        self.Swin       = np.zeros(tsteps)    # incoming short wave radiation [W m-2]
        self.Swout      = np.zeros(tsteps)    # outgoing short wave radiation [W m-2]
        self.Lwin       = np.zeros(tsteps)    # incoming long wave radiation [W m-2]
        self.Lwout      = np.zeros(tsteps)    # outgoing long wave radiation [W m-2]
        self.Q          = np.zeros(tsteps)    # net radiation [W m-2]

        # land surface variables
        self.ra         = np.zeros(tsteps)    # aerodynamic resistance [s m-1]
        self.rs         = np.zeros(tsteps)    # surface resistance [s m-1]
        self.H          = np.zeros(tsteps)    # sensible heat flux [W m-2]
        self.LE         = np.zeros(tsteps)    # evapotranspiration [W m-2]
        self.LEliq      = np.zeros(tsteps)    # open water evaporation [W m-2]
        self.LEveg      = np.zeros(tsteps)    # transpiration [W m-2]
        self.LEsoil     = np.zeros(tsteps)    # soil evaporation [W m-2]
        self.LEpot      = np.zeros(tsteps)    # potential evaporation [W m-2]
        self.LEref      = np.zeros(tsteps)    # reference evaporation at rs = rsmin / LAI [W m-2]
        self.G          = np.zeros(tsteps)    # ground heat flux [W m-2]

        # cumulus variables
        self.ac         = np.zeros(tsteps)    # cloud core fraction [-]
        self.M          = np.zeros(tsteps)    # cloud core mass flux [m s-1]

# class for storing mixed-layer model input data
class model_input:
    def __init__(self):
        # general model variables
        self.runtime    = None  # duration of model run [s]
        self.dt         = None  # time step [s]

        # mixed-layer variables
        self.sw_ml      = None  # mixed-layer model switch
        self.sw_shearwe = None  # Shear growth ABL switch
        self.h          = None  # initial ABL height [m]
        self.Ps         = None  # surface pressure [Pa]
        self.ws         = None  # large scale vertical velocity [m s-1]
        self.fc         = None  # Coriolis parameter [s-1]
        
        self.theta      = None  # initial mixed-layer potential temperature [K]
        self.dtheta     = None  # initial temperature jump at h [K]
        self.gammatheta = None  # free atmosphere potential temperature lapse rate [K m-1]
        self.advtheta   = None  # advection of heat [K s-1]
        self.beta       = None  # entrainment ratio for virtual heat [-]
        self.wtheta     = None  # surface kinematic heat flux [K m s-1]
        
        self.q          = None  # initial mixed-layer specific humidity [kg kg-1]
        self.dq         = None  # initial specific humidity jump at h [kg kg-1]
        self.gammaq     = None  # free atmosphere specific humidity lapse rate [kg kg-1 m-1]
        self.advq       = None  # advection of moisture [kg kg-1 s-1]
        self.wq         = None  # surface kinematic moisture flux [kg kg-1 m s-1]

        self.CO2        = None  # initial mixed-layer potential temperature [K]
        self.dCO2       = None  # initial temperature jump at h [K]
        self.gammaCO2   = None  # free atmosphere potential temperature lapse rate [K m-1]
        self.advCO2     = None  # advection of heat [K s-1]
        self.wCO2       = None  # surface kinematic heat flux [K m s-1]
        
        self.sw_wind    = None  # prognostic wind switch
        self.u          = None  # initial mixed-layer u-wind speed [m s-1]
        self.du         = None  # initial u-wind jump at h [m s-1]
        self.gammau     = None  # free atmosphere u-wind speed lapse rate [s-1]
        self.advu       = None  # advection of u-wind [m s-2]

        self.v          = None  # initial mixed-layer u-wind speed [m s-1]
        self.dv         = None  # initial u-wind jump at h [m s-1]
        self.gammav     = None  # free atmosphere v-wind speed lapse rate [s-1]
        self.advv       = None  # advection of v-wind [m s-2]

        # surface layer variables
        self.sw_sl      = None  # surface layer switch
        self.ustar      = None  # surface friction velocity [m s-1]
        self.z0m        = None  # roughness length for momentum [m]
        self.z0h        = None  # roughness length for scalars [m]
        self.Cm         = None  # drag coefficient for momentum [-]
        self.Cs         = None  # drag coefficient for scalars [-]
        self.L          = None  # Obukhov length [-]
        self.Rib        = None  # bulk Richardson number [-]

        # radiation parameters
        self.sw_rad     = None  # radiation switch
        self.lat        = None  # latitude [deg]
        self.lon        = None  # longitude [deg]
        self.doy        = None  # day of the year [-]
        self.tstart     = None  # time of the day [h UTC]
        self.cc         = None  # cloud cover fraction [-]

        # land surface parameters
        self.sw_ls      = None  # land surface switch
        self.ls_type    = None  # land-surface parameterization ('js' for Jarvis-Stewart or 'ags' for A-Gs)
        self.wg         = None  # volumetric water content top soil layer [m3 m-3]
        self.w2         = None  # volumetric water content deeper soil layer [m3 m-3]
        self.Tsoil      = None  # temperature top soil layer [K]
        self.T2         = None  # temperature deeper soil layer [K]
        
        self.a          = None  # Clapp and Hornberger retention curve parameter a
        self.b          = None  # Clapp and Hornberger retention curve parameter b
        self.p          = None  # Clapp and Hornberger retention curve parameter p 
        self.CGsat      = None  # saturated soil conductivity for heat
        
        self.wsat       = None  # saturated volumetric water content ECMWF config [-]
        self.wfc        = None  # volumetric water content field capacity [-]
        self.wwilt      = None  # volumetric water content wilting point [-]
        
        self.C1sat      = None 
        self.C2ref      = None
        
        self.LAI        = None  # leaf area index [-]
        self.gD         = None  # correction factor transpiration for VPD [-]
        self.rsmin      = None  # minimum resistance transpiration [s m-1]
        self.rssoilmin  = None  # minimum resistance soil evaporation [s m-1]
        self.alpha      = None  # surface albedo [-]
        
        self.Ts         = None  # initial surface temperature [K]
        
        self.cveg       = None  # vegetation fraction [-]
        self.Wmax       = None  # thickness of water layer on wet vegetation [m]
        self.Wl         = None  # equivalent water layer depth for wet vegetation [m]
        
        self.Lambda     = None  # thermal diffusivity skin layer [-]

        # A-Gs parameters
        self.c3c4       = None  # Plant type ('c3' or 'c4')

        # Cumulus parameters
        self.sw_cu      = None  # Cumulus parameterization switch
        self.dz_h       = None  # Transition layer thickness [m]

class empty:
    def __init__(self):
        pass

if(__name__ == "__main__"):
    
    r1in = model_input()
    
    r1in.dt         = 60.       # time step [s]
    r1in.runtime    = 12*3600    # total run time [s]
    
    # mixed-layer input
    r1in.sw_ml      = True      # mixed-layer model switch
    r1in.sw_shearwe = False     # shear growth mixed-layer switch
    r1in.h          = 200.      # initial ABL height [m]
    r1in.Ps         = 101300.   # surface pressure [Pa]
    r1in.ws         = 0.        # large scale vertical velocity [m s-1]
    r1in.fc         = 1.e-4     # Coriolis parameter [m s-1]
    
    r1in.theta      = 288.      # initial mixed-layer potential temperature [K]
    r1in.dtheta     = 1.        # initial temperature jump at h [K]
    r1in.gammatheta = 0.006     # free atmosphere potential temperature lapse rate [K m-1]
    r1in.advtheta   = 0.        # advection of heat [K s-1]
    r1in.beta       = 0.2       # entrainment ratio for virtual heat [-]
    r1in.wtheta     = 0.1       # surface kinematic heat flux [K m s-1]
    
    r1in.q          = 0.008     # initial mixed-layer specific humidity [kg kg-1]
    r1in.dq         = -0.001    # initial specific humidity jump at h [kg kg-1]
    r1in.gammaq     = 0.        # free atmosphere specific humidity lapse rate [kg kg-1 m-1]
    r1in.advq       = 0.        # advection of moisture [kg kg-1 s-1]
    r1in.wq         = 0.1e-3    # surface kinematic moisture flux [kg kg-1 m s-1]
   
    r1in.CO2        = 422.      # initial mixed-layer CO2 [ppm]
    r1in.dCO2       = -44.      # initial CO2 jump at h [ppm]
    r1in.gammaCO2   = 0.        # free atmosphere CO2 lapse rate [ppm m-1]
    r1in.advCO2     = 0.        # advection of CO2 [ppm s-1]
    r1in.wCO2       = 10.       # surface kinematic CO2 flux [ppm m s-1]
    
    r1in.sw_wind    = False      # prognostic wind switch
    r1in.u          = 6.        # initial mixed-layer u-wind speed [m s-1]
    r1in.du         = 4.        # initial u-wind jump at h [m s-1]
    r1in.gammau     = 0.        # free atmosphere u-wind speed lapse rate [s-1]
    r1in.advu       = 0.        # advection of u-wind [m s-2]
    
    r1in.v          = -4.0      # initial mixed-layer u-wind speed [m s-1]
    r1in.dv         = 4.0       # initial u-wind jump at h [m s-1]
    r1in.gammav     = 0.        # free atmosphere v-wind speed lapse rate [s-1]
    r1in.advv       = 0.        # advection of v-wind [m s-2]
    
    # surface layer input
    r1in.sw_sl      = False      # surface layer switch
    r1in.ustar      = 0.3       # surface friction velocity [m s-1]
    r1in.z0m        = 0.02      # roughness length for momentum [m]
    r1in.z0h        = 0.002     # roughness length for scalars [m]
    
    # radiation parameters
    r1in.sw_rad     = True      # radiation switch
    r1in.lat        = 51.97     # latitude [deg]
    r1in.lon        = -4.93     # longitude [deg]
    r1in.doy        = 268.      # day of the year [-]
    r1in.tstart     = 6.8       # time of the day [h UTC]
    r1in.cc         = 0.0       # cloud cover fraction [-]
    r1in.Q          = 400.      # net radiation [W m-2] 
    
    # land surface parameters
    r1in.sw_ls      = True      # land surface switch
    r1in.ls_type    = 'ags'     # land-surface parameterization ('js' for Jarvis-Stewart or 'ags' for A-Gs)
    r1in.wg         = 0.21      # volumetric water content top soil layer [m3 m-3]
    r1in.w2         = 0.21      # volumetric water content deeper soil layer [m3 m-3]
    r1in.cveg       = 0.85      # vegetation fraction [-]
    r1in.Tsoil      = 285.      # temperature top soil layer [K]
    r1in.T2         = 286.      # temperature deeper soil layer [K]
    r1in.a          = 0.219     # Clapp and Hornberger retention curve parameter a
    r1in.b          = 4.90      # Clapp and Hornberger retention curve parameter b
    r1in.p          = 4.        # Clapp and Hornberger retention curve parameter c
    r1in.CGsat      = 3.56e-6   # saturated soil conductivity for heat
    
    r1in.wsat       = 0.472     # saturated volumetric water content ECMWF config [-]
    r1in.wfc        = 0.323     # volumetric water content field capacity [-]
    r1in.wwilt      = 0.171     # volumetric water content wilting point [-]
    
    r1in.C1sat      = 0.132     
    r1in.C2ref      = 1.8
    
    r1in.LAI        = 2.        # leaf area index [-]
    r1in.gD         = 0.0       # correction factor transpiration for VPD [-]
    r1in.rsmin      = 110.      # minimum resistance transpiration [s m-1]
    r1in.rssoilmin  = 50.       # minimun resistance soil evaporation [s m-1]
    r1in.alpha      = 0.25      # surface albedo [-]
    
    r1in.Ts         = 290.      # initial surface temperature [K]
    
    r1in.Wmax       = 0.0002    # thickness of water layer on wet vegetation [m]
    r1in.Wl         = 0.0000    # equivalent water layer depth for wet vegetation [m]
    
    r1in.Lambda     = 5.9       # thermal diffusivity skin layer [-]

    # A-Gs parameters
    r1in.c3c4       = 'c3'      # Plant type ('c3' or 'c4')

    # Cloud parameters
    r1in.sw_cu      = False      # Cumulus parameterization switch
    r1in.dz_h       = 150.      # Transition layer thickness [m]
    
    r1 = model(r1in)
    r1.run()

    from read_class import *
    import sys

    rr = ReadCLASS('/Users/m300241/Desktop/run1.csv') 

    from pylab import *
    close('all')

    if(False):
        figure()
        subplot(331)
        plot(r1.out.t, r1.out.h, 'b-', label='h python')
        plot(rr.timeUTC, rr.h, 'g-', label='CLASS')
        legend(frameon=False)

        subplot(332)
        plot(r1.out.t, r1.out.theta, 'b-', label='theta python')
        plot(rr.timeUTC, rr.th, 'g-', label='CLASS')
        legend(frameon=False)

        subplot(333)
        plot(r1.out.t, r1.out.q*1000, 'b-', label='q python')
        plot(rr.timeUTC, rr.q*1000, 'g-', label='CLASS')
        legend(frameon=False)

        subplot(334)
        plot(r1.out.t, r1.out.Swin, 'b-', label='swin python')
        plot(rr.timeUTC, rr.Swin, 'g-', label='CLASS')
        plot(r1.out.t, r1.out.Swout, 'b--', label='swout python')
        plot(rr.timeUTC, rr.Swout, 'g--', label='CLASS')
        legend(frameon=False)

        #subplot(335)
        #plot(r1.out.t, r1.out.CO2, 'b-', label='CO2 python')
        #plot(rr.timeUTC, rr.CO2-1., 'g-', label='CLASS')
        #legend(frameon=False)

        #subplot(336)
        #plot(r1.out.t, r1.out.wCO2, 'b-', label='wCO2 python')
        #plot(rr.timeUTC, rr.wCO2s, 'g-', label='CLASS')
        #legend(frameon=False)

        #subplot(337)
        #plot(r1.out.t, r1.out.wCO2A, 'b-', label='Assym python')
        #plot(rr.timeUTC, rr.wCO2A, 'g-', label='CLASS')
        #legend(frameon=False)

        #subplot(338)
        #plot(r1.out.t, r1.out.wCO2R, 'b-', label='Reso python')
        #plot(rr.timeUTC, rr.wCO2R, 'g-', label='CLASS')
        #legend(frameon=False)

        #subplot(339)
        #plot(r1.out.t, r1.out.wCO2e, 'b-', label='wCO2e python')
        #plot(rr.timeUTC, rr.wCO2e, 'g-', label='CLASS')
        #legend(frameon=False)

    # A-Gs comparison
    if(True):
        figure()
        subplot(321)
        plot(r1.out.t,   r1.out.CO2, 'g-', label='python')
        plot(rr.timeUTC, rr.CO2, 'k-', label='CLASS')
        ylabel('CO2')
        legend(frameon=False)

        subplot(322)
        plot(r1.out.t,   r1.out.dCO2, 'g-', label='python')
        plot(rr.timeUTC, rr.dCO2, 'k-', label='CLASS')
        ylabel('dCO2')

        subplot(323)
        plot(r1.out.t,   r1.out.wCO2, 'g-', label='python')
        plot(rr.timeUTC, rr.wCO2s, 'k-', label='CLASS')
        ylabel('wCO2s')

        subplot(324)
        plot(r1.out.t,   r1.out.wCO2e, 'g-', label='python')
        plot(rr.timeUTC, rr.wCO2e, 'k-', label='CLASS')
        ylabel('wCO2e')

        subplot(325)
        plot(r1.out.t,   r1.out.wCO2A, 'g-', label='assym python')
        plot(rr.timeUTC, rr.wCO2A, 'k-', label='CLASS')
        plot(r1.out.t,   r1.out.wCO2R, 'g--', label='assym python')
        plot(rr.timeUTC, rr.wCO2R, 'k--', label='CLASS')
        ylabel('wCO2A/R')