#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  7 10:51:03 2017

@author: vsc42247

Purpose: Set surface conditions for the CLASS boundary-layer model
"""


import netCDF4 as nc4
import numpy as np
import datetime as dt
#you can install with
import pynacolada as pcd
import pandas as pd

def get_class4gl_ground(class_settings,**kwargs):   
    
    key = "IGBPDIS"
    if ((kwargs == {}) or ((key in kwargs.keys()) and (kwargs[key]))):
    
        
        input_fn = "/user/data/gent/gvo000/gvo00090/EXT/data/IGBP-DIS/wsat.nc"
        print('reading soil water saturation from '+input_fn)

        input_nc = nc4.Dataset(input_fn,'r')
        ilat = np.where(input_nc.variables['lat'][:] >= class_settings.lat)[0][-1]
        ilon = np.where(input_nc.variables['lon'][:] >= class_settings.lon)[0][0]
        class_settings.__dict__['wsat'] = input_nc.variables['wsat'][ilon,ilat]
        input_nc.close()

        input_fn = "/user/data/gent/gvo000/gvo00090/EXT/data/IGBP-DIS/wfc.nc"
        print('reading soil water field capacity from '+input_fn)
    
        input_nc = nc4.Dataset(input_fn,'r')
        ilat = np.where(input_nc.variables['lat'][:] >= class_settings.lat)[0][-1]
        ilon = np.where(input_nc.variables['lon'][:] >= class_settings.lon)[0][0]
        class_settings.__dict__['wfc'] = input_nc.variables['wfc'][ilon,ilat]
        input_nc.close()
        
        input_fn = "/user/data/gent/gvo000/gvo00090/EXT/data/IGBP-DIS/wwp.nc"
        print('reading soil wilting point from '+input_fn)
        input_nc = nc4.Dataset(input_fn,'r')
        ilat = np.where(input_nc.variables['lat'][:] >= class_settings.lat)[0][-1]
        ilon = np.where(input_nc.variables['lon'][:] >= class_settings.lon)[0][0]
        class_settings.__dict__['wwilt'] = input_nc.variables['wwp'][ilon,ilat]
        input_nc.close()
        
    key = "GLEAM"
    if ((kwargs == {}) or ((key in kwargs.keys()) and (kwargs[key]))):
        
        #INPUT_gleam = gleam() 
        #INPUT_gleam.path = "/kyukon/data/gent/gvo000/gvo00090/GLEAM/data/v3.1a/" 
        
        gleam_path = "/user/data/gent/gvo000/gvo00090/GLEAM/data/v3.1a/"
        print('reading soil-water content for "+str(class_settings,datetime.year)+" from '+gleam_path)
        
        gleam_files = {}
        
        gleam_vars = ['SMroot','SMsurf']
        
        for VAR in gleam_vars:
            gleam_files[VAR] = nc4.Dataset(gleam_path+'/'+str(class_settings.datetime.year)+'/'+VAR+'_'+str(class_settings.datetime.year)+'_GLEAM_v3.1a.nc','r')
        

        year = class_settings.datetime.year
        day = class_settings.datetime.day
        hour = class_settings.datetime.hour
  
        ilat = np.where(gleam_files['SMsurf'].variables['lat'][:] >= class_settings.lat)[0][-1]
        ilon = np.where(gleam_files['SMsurf'].variables['lon'][:] >= class_settings.lon)[0][0]
        
        VAR = 'SMsurf'; class_settings.wg = gleam_files[VAR].variables[VAR][day-1,ilon,ilat]
        VAR = 'SMroot'; class_settings.w2 = gleam_files[VAR].variables[VAR][day-1,ilon,ilat]
        
        for VAR in gleam_vars:
            gleam_files[VAR].close()
    
    key = "MOD44B"
    if ((kwargs == {}) or ((key in kwargs.keys()) and (kwargs[key]))):
    
        input_fn = "/user/data/gent/gvo000/gvo00090/EXT/data/MOD44B/fv.nc"
        print('initializing vegetation fraction from '+input_fn)
        var = 'cveg'
        input_nc = nc4.Dataset(input_fn,'r')
        ilat = np.where(input_nc.variables['lat'][:] >= class_settings.lat)[0][-1]
        ilon = np.where(input_nc.variables['lon'][:] >= class_settings.lon)[0][0]
        class_settings.__dict__[var] = input_nc.variables['fv'][ilon,ilat]
        input_nc.close()
        
    key = "DSMW"
    if ((kwargs == {}) or ((key in kwargs.keys()) and (kwargs[key]))):
         # Procedure of the thermal properties:
         # 1. determine soil texture from DSMW
         # 2. soil type with look-up table (according to DWD/EXTPAR)
         # 3. Thermal properties used in the force-restore method (Clapp and Hornberger, 1987) 
         #    with parameter look-up table from Noilhan and Planton (1989). 
         #    Note: The look-up table is inspired on DWD/COSMO
                 
       
        #preparing for soil thermal properties
        input_fn = "/user/data/gent/gvo000/gvo00090/EXT/data/DSMW/FAO_DSMW_DP.nc"
        
        print("deriving soil thermal properties for the force-restore methodes from the soil texture file "+ input_fn)
        
        input_nc = nc4.Dataset(input_fn,'r')
        ilat = np.where(input_nc.variables['lat'][:] >= class_settings.lat)[0][-1]
        ilon = np.where(input_nc.variables['lon'][:] >= class_settings.lon)[0][0]
        DSMW = input_nc.variables['DSMW'][ilat,ilon]
        
        
        #EXTPAR: zfine   = soil_texslo(soil_unit)%tex_fine
        SP = {}; SPKEYS = ['tex_coarse', 'tex_medium', 'tex_fine', 'code']
        for SPKEY in SPKEYS: 
            SP[SPKEY] = np.array(input_nc.variables[SPKEY][DSMW])
        input_nc.close()
        
        SP['texture'] = (0.5*SP['tex_medium']+1.0*SP['tex_coarse']) /(SP['tex_coarse']+SP['tex_medium']+SP['tex_fine'])
        
        if pd.isnull(SP['texture']):
            print('Warning, texture is invalid> Setting to Ocean')
            SP['itex'] = 9
        
        else:
            SP['itex'] = int(SP['texture']*100)
        
        #adopted from mo_agg_soil.f90 (EXTPAR3.0)
        SP['isoil'] = np.zeros_like(SP['itex'],dtype=np.int)
        LOOKUP = [
                  [0 ,7],# fine textured, clay (soil type 7)
                  [20,6],# medium to fine textured, loamy clay (soil type 6)
                  [40,5],# medium textured, loam (soil type 5)
                  [60,4],# coarse to medium textured, sandy loam (soil type 4)
                  [80,3],# coarse textured, sand (soil type 3)
                ]
        for iitex,iisoil in LOOKUP: 
            SP['isoil'][SP['itex'] >= iitex ] = iisoil 
        
        #adopted from mo_agg_soil.f90 (EXTPAR3.0)
        LOOKUP = [
                  [9001, 1 ], # ice, glacier (soil type 1) 
                  [9002, 2 ], # rock, lithosols (soil type 2)
                  [9003, 3 ], # salt, set soiltype to sand (soil type 3)
                  [9004, 8 ], # histosol, e.g. peat (soil type 8)
                  [9,    9 ], # undefined (ocean)
                  [9005, 3 ], # shifting sands or dunes, set soiltype to sand (soil type 3)
                  [9000, 9 ], # undefined (inland lake)
                  [9009, 5 ], #  default_soiltype ! undefined (nodata), set soiltype to loam (soil type )
                  [9012, 5 ], #  default_soiltype undefined (dominant part undefined), set soiltype to loam (soil type 5)
                ]
        # EXTPAR: soil_code = soil_texslo(soil_unit)%dsmw_code # the legend has some special cases for the "soil_code"
        for icode,iisoil in LOOKUP: 
            SP['isoil'][SP['code'] == icode] = iisoil 
        
        #adopted from data_soil.f90 (COSMO5.0)
        SP_LOOKUP = { 
          # soil type:         ice        rock       sand        sandy      loam         clay        clay        peat        sea        sea  
          # (by index)                                           loam                    loam                                water      ice
          'cporv'  : [ np.nan, 1.E-10   , 1.E-10   , 0.364     , 0.445     , 0.455     , 0.475     , 0.507     , 0.863     , 1.E-10   , 1.E-10   ],
          'cfcap'  : [ np.nan, 1.E-10   , 1.E-10   , 0.196     , 0.260     , 0.340     , 0.370     , 0.463     , 0.763     , 1.E-10   , 1.E-10   ],
          'cpwp'   : [ np.nan, 0.0      , 0.0      , 0.042     , 0.100     , 0.110     , 0.185     , 0.257     , 0.265     , 0.0      ,  0.0     ],
          'cadp'   : [ np.nan, 0.0      , 0.0      , 0.012     , 0.030     , 0.035     , 0.060     , 0.065     , 0.098     , 0.0      ,  0.0     ],
          'crhoc'  : [ np.nan, 1.92E6   , 2.10E6   , 1.28E6    , 1.35E6    , 1.42E6    , 1.50E6    , 1.63E6    , 0.58E6    , 4.18E6   , 1.92E6   ],
          'cik2'   : [ np.nan, 0.0      , 0.0      , 0.0035    , 0.0023    , 0.0010    , 0.0006    , 0.0001    , 0.0002    , 0.0      ,  0.0     ],
          'ckw0'   : [ np.nan, 0.0      , 0.0      , 479.E-7   , 943.E-8   , 531.E-8   , 764.E-9   , 17.E-9    , 58.E-9    , 0.0      ,  0.0     ],
          'ckw1'   : [ np.nan, 0.0      , 0.0      , -19.27    , -20.86    , -19.66    , -18.52    , -16.32    , -16.48    , 0.0      ,  0.0     ],
          'cdw0'   : [ np.nan, 0.0      , 0.0      , 184.E-7   , 346.E-8   , 357.E-8   , 118.E-8   , 442.E-9   , 106.E-9   , 0.0      ,  0.0     ],
          'cdw1'   : [ np.nan, 0.0      , 0.0      , -8.45     , -9.47     , -7.44     , -7.76     , -6.74     , -5.97     , 0.0      ,  0.0     ],
          'crock'  : [ np.nan, 0.0      , 0.0      , 1.0       , 1.0       , 1.0       , 1.0       , 1.0       , 1.0       , 0.0      ,  0.0     ],
          'cala0'  : [ np.nan, 2.26     , 2.41     , 0.30      , 0.28      , 0.25      , 0.21      , 0.18      , 0.06      , 1.0      ,  2.26    ],
          'cala1'  : [ np.nan, 2.26     , 2.41     , 2.40      , 2.40      , 1.58      , 1.55      , 1.50      , 0.50      , 1.0      ,  2.26    ],
          'csalb'  : [ np.nan, 0.70     , 0.30     , 0.30      , 0.25      , 0.25      , 0.25      , 0.25      , 0.20      , 0.07     ,  0.70    ],
          'csalbw' : [ np.nan, 0.00     , 0.00     , 0.44      , 0.27      , 0.24      , 0.23      , 0.22      , 0.10      , 0.00     ,  0.00    ],
          'ck0di'  : [ np.nan, 1.E-4    , 1.E-4    , 2.E-4     , 2.E-5     , 6.E-6     , 2.E-6     , 1.E-6     , 1.5E-6    , 0.00     ,  0.00    ],
          'cbedi'  : [ np.nan, 1.00     , 1.00     , 3.5       , 4.8       , 6.1       , 8.6       , 10.0      , 9.0       , 0.00     ,  0.00    ],
          'csandf' : [ np.nan, 0.0      , 0.0      , 90.       , 65.       , 40.       , 35.       , 15.       , 90.       , 0.00     ,  0.00    ],
          'cclayf' : [ np.nan, 0.0      , 0.0      , 5.0       , 10.       , 20.       , 35.       , 70.       , 5.0       , 0.00     ,  0.00    ],
          #supplement Noihhan andf Planton 1989 soil texture parameters for the force-restore method.
          'b'      : [ np.nan, np.nan   , np.nan   , 4.05      , 4.90      , 5.39      , 8.52      , 11.40     , np.nan    , np.nan   ,  np.nan  ],
          #error in table 2 of NP89: values need to be multiplied by e-6
          'CGsat'  : [ np.nan, np.nan   , np.nan   , 3.222e-6     , 3.560e-6     , 4.111e-6     , 3.995e-6     , 3.600e-6     , np.nan    , np.nan   ,  np.nan  ],
          'p'  :     [ np.nan, np.nan   , np.nan   , 4.        , 4.        , 6.        , 10.       , 12.       , np.nan    , np.nan   ,  np.nan  ],
          'a'  :     [ np.nan, np.nan   , np.nan   , 0.387     , 0.219     , 0.148     , 0.084     , 0.083     , np.nan    , np.nan   ,  np.nan  ],
          'C1sat'  : [ np.nan, np.nan   , np.nan   , 0.082     , 0.132     , 0.191     , 0.227     , 0.342     , np.nan    , np.nan   ,  np.nan  ],
          'C2ref'  : [ np.nan, np.nan   , np.nan   , 3.9       , 1.8       , 0.8       , 0.6       , 0.3       , np.nan    , np.nan   ,  np.nan  ],
        }
        
        for SPKEY in SP_LOOKUP.keys(): 
            SP[SPKEY] = np.zeros_like(SP['isoil'],dtype=np.float)
        
        for i in range(11):
            SELECT = (SP['isoil'] == i)
            for SPKEY in SP_LOOKUP.keys(): 
                SP[SPKEY][SELECT] = SP_LOOKUP[SPKEY][i]
        
        for SPKEY in list(SP_LOOKUP.keys())[-6:]: 
            var = SPKEY
            class_settings.__dict__[var] = np.float(SP[SPKEY])
            
        # only print the last parameter value in the plot
        
        #inputs.append(cp.deepcopy(class_settings))
        #var = 'cala'
        #class_settings.__dict__[var] = np.float(SP['cala0'])
        #valnew = class_settings.__dict__[var]
        #labels.append(var+': '+format(valold,"0.2g")+'->'+format(valnew,"0.2g"))
        
        #inputs.append(cp.deepcopy(class_settings))
        #var = 'crhoc'
        #class_settings.__dict__[var] = np.float(SP['crhoc'])
        #valnew = class_settings.__dict__[var]
        #labels.append(var+': '+format(valold,"0.2g")+'->'+format(valnew,"0.2g"))
        
    key = "CERES"
    if ((kwargs == {}) or ((key in kwargs.keys()) and (kwargs[key]))):    
        
        CERES_start_date = dt.datetime(2000,3,1)
        DT_CERES_START = (CERES_start_date + dt.timedelta(days=(int((class_settings.datetime - CERES_start_date ).days/61) * 61)))
        DT_CERES_END   = DT_CERES_START +dt.timedelta(days=60)
        
        input_fn = "/user/data/gent/gvo000/gvo00090/EXT/data/CERES/CERES_SYN1deg-1H_Terra-Aqua-MODIS_Ed4A_Subset_"+DT_CERES_START.strftime("%Y%m%d")+"-"+DT_CERES_END.strftime("%Y%m%d")+".nc"
        print("Reading afternoon cloud cover for "+str(class_settings.datetime)+" from "+input_fn)
            
        var = 'cc'
        
        input_nc = nc4.Dataset(input_fn,'r')
        
        idatetime = np.where(np.array(pcd.ncgetdatetime(input_nc))  >= class_settings.datetime)[0][0]
        
        ilat = np.where(input_nc.variables['lat'][:] >= class_settings.lat)[0][-1]
        ilon = np.where(input_nc.variables['lon'][:] >= class_settings.lon)[0][0]
        print(class_settings.lat,class_settings.lon)
        
        class_settings.__dict__[var] = np.nanmean(input_nc.variables['cldarea_total_1h'][idatetime:(idatetime+class_settings.runtime),ilat,ilon])/100.
   
        input_nc.close()
    
    key = "GIMMS"
    if ((kwargs == {}) or ((key in kwargs.keys()) and (kwargs[key]))):    
       
    
        input_fn = "/user/data/gent/gvo000/gvo00090/EXT/data/GIMMS/v2/LAI/gimms-3g.v2.lai.1981-2015_monmean.nc"
        print("Reading Leag Area Index from "+input_fn)
        var = 'LAI'
        
        #plt.plot
        
        input_nc = nc4.Dataset(input_fn,'r')
        
        #idatetime = np.where(np.array(pcd.ncgetdatetime(input_nc))  >= class_settings.datetime)[0][0]
        idatetime = np.where(np.array(pcd.ncgetdatetime(input_nc))  >= class_settings.datetime)[0][0]
        
        ilatitude = np.where(input_nc.variables['lat'][:] >= class_settings.lat)[0][-1]
        ilongitude = np.where(input_nc.variables['lon'][:] >= class_settings.lon)[0][0]
        
        # divide by cveg, since it only reflects the LAI for the vegetation fraction and not for the entire (satellite) grid cell
        
        print('Warning! Dividing by cveg, which is: '+str(class_settings.cveg))
        tarray = np.array(input_nc.variables['LAI'][:,ilatitude,ilongitude])/class_settings.cveg
        
        if np.isnan(tarray[idatetime]):
            print("interpolating GIMMS cveg nan value")
            
            mask = np.isnan(tarray)
            if np.where(mask)[0].shape[0] < 0.25*mask.shape[0]:
                tarray[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), tarray[~mask])
            else:
                print("Warning. Could not interpolate GIMMS cveg nan value")
                
        class_settings.__dict__[var] = tarray[idatetime]
        
        input_nc.close()
 
    key = "IGBPDIS_ALPHA"
    if ((kwargs == {}) or ((key in kwargs.keys()) and (kwargs[key]))):          
       
        var = 'alpha'
        
        input_fn = "/user/data/gent/gvo000/gvo00090/EXT/data/IGBP-DIS/FRACTIONS_GLEAMv31a.nc"
        print("Reading albedo from "+input_fn)
    
        input_nc = nc4.Dataset(input_fn,'r')
        ilat = np.where(input_nc.variables['lat'][:] >= class_settings.lat)[0][-1]
        ilon = np.where(input_nc.variables['lon'][:] >= class_settings.lon)[0][0]
        
        
        landfr = {}
        for ltype in ['W','B','H','TC']:   
            landfr[ltype] = input_nc.variables['f'+ltype][0,ilon,ilat]
        
        aweights = {'W':0.075,'TC':0.15,'H':0.22,'B':0.30}
        
        alpha=0.
        for ltype in landfr.keys():
            alpha += landfr[ltype]*aweights[ltype]
        
        
        class_settings.__dict__[var] = alpha
        input_nc.close()        
        
        
    key = "ERAINT_ST"
    if ((kwargs == {}) or ((key in kwargs.keys()) and (kwargs[key]))):          
       
        input_fn = '/user/data/gent/gvo000/gvo00090/EXT/data/ERA-INTERIM/by_var_nc/stl1_3hourly/stl1_'+str(class_settings.datetime.year)+"_3hourly.nc"
        print("Reading soil temperature from "+input_fn)
        
        var = 'Tsoil'
        input_nc = nc4.Dataset(input_fn,'r')
        
        idatetime = np.where(np.array(pcd.ncgetdatetime(input_nc))  >= class_settings.datetime)[0][0]
        
        ilatitude = np.where(input_nc.variables['latitude'][:] >= class_settings.lat)[0][-1]
        ilongitude = np.where(input_nc.variables['longitude'][:] >= class_settings.lon)[0][0]
        
        
        class_settings.__dict__[var] = input_nc.variables['stl1'][idatetime,ilatitude,ilongitude]
        
        input_fn = '/user/data/gent/gvo000/gvo00090/EXT/data/ERA-INTERIM/by_var_nc/stl2_3hourly/stl2_'+str(class_settings.datetime.year)+"_3hourly.nc"
        var = 'T2'
        
        input_nc = nc4.Dataset(input_fn,'r')
        
        idatetime = np.where(np.array(pcd.ncgetdatetime(input_nc))  >= class_settings.datetime)[0][0]
        
        ilatitude = np.where(input_nc.variables['latitude'][:] >= class_settings.lat)[0][-1]
        ilongitude = np.where(input_nc.variables['longitude'][:] >= class_settings.lon)[0][0]
        
        
        class_settings.__dict__[var] = input_nc.variables['stl2'][idatetime,ilatitude,ilongitude]
        
        
        input_nc.close()
        
        
    
    #inputs.append(cp.deepcopy(class_settings))
    #var = 'T2'
    #valold = class_settings.__dict__[var]
    #
    #class_settings.__dict__[var] = 305.
    #class_settings.__dict__['Tsoil'] = 302.
    #valnew = class_settings.__dict__[var]
    #labels.append(var+': '+format(valold,"0.2g")+'->'+format(valnew,"0.2g"))
    
    
    
    #inputs.append(cp.deepcopy(class_settings))
    #
    #var = 'Lambda'
    #valold = class_settings.__dict__[var]
    
    ## I presume that the skin layer conductivity scales with both LAI and vegetation fraction, which seems ~ valid according to table 10.6 in CLASS-book. 
    ## I need to ask Chiel.
    ## I extrapolate from Lambda value of grass with Lambda = 5.9 W m-2 K-1, LAI = 2 and cveg = 0.85
    #
    #valnew = 5.9 / 2. / 0.85 * class_settings.__dict__['LAI'] * class_settings.__dict__['cveg'] 
    #class_settings.__dict__[var] = valnew
    #labels.append(var+': '+format(valold,"0.2g")+'->'+format(valnew,"0.2g"))
    
    
    
    key = "GLAS"
    if ((kwargs == {}) or ((key in kwargs.keys()) and (kwargs[key]))):          
       
        input_fn = "/user/data/gent/gvo000/gvo00090/EXT/data/GLAS/global_canopy_height_0.25.nc"
        print("Reading canopy height for determining roughness length from "+input_fn)
        var = 'z0m'
    
        
        #plt.plot
        
        input_nc = nc4.Dataset(input_fn,'r')
        
        ilat = np.where(input_nc.variables['lat'][:] >= class_settings.lat)[0][0]
        ilon = np.where(input_nc.variables['lon'][:] >= class_settings.lon)[0][0]
        
        testval = np.float64(input_nc.variables['Band1'][ilat,ilon])/10.
        
        lowerlimit = 0.01
        if testval < lowerlimit:
            print('forest canopy height very very small. We take a value of '+str(lowerlimit))
            class_settings.__dict__[var] = lowerlimit
        else:
            class_settings.__dict__[var] = testval
        
        class_settings.__dict__['z0h'] =  class_settings.__dict__['z0m']/10.
        
        
        input_nc.close()
        
