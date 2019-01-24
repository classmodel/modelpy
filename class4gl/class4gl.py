# -*- coding: utf-8 -*-

"""

Created on Mon Jan 29 12:33:51 2018

Module file for class4gl, which  extents the class-model to be able to take
global air profiles as input. It exists of:

CLASSES:
    - an input object, namely class4gl_input. It includes:
        - a function to read Wyoming sounding data from a yyoming stream object
        - a function to read global data from a globaldata library object 
    - the model object: class4gl
    - ....    

DEPENDENCIES:
    - xarray
    - numpy
    - data_global
    - Pysolar
    - yaml

@author: Hendrik Wouters

"""



""" Setup of envirnoment """

# Standard modules of the stand class-boundary-layer model
from model import model
from model import model_output as class4gl_output
from model import model_input
from model import qsat
#from data_soundings import wyoming 

import importlib
spam_loader = importlib.find_loader('Pysolar')
found = spam_loader is not None
if found:
    import Pysolar
    import Pysolar.util as Pysolarutil
    GetSunriseSunset = Pysolarutil.GetSunriseSunset
    GetAzimuth = Pysolarutil.solar.GetAzimuth
    GetAltitude = Pysolarutil.solar.GetAltitude
else:
    import pysolar as Pysolar
    Pysolarutil = Pysolar.util
    GetSunriseSunset = Pysolarutil.get_sunrise_sunset
    GetAzimuth = Pysolar.solar.get_azimuth
    GetAltitude = Pysolar.solar.get_altitude
import yaml
import logging
import warnings
import pytz

#formatter = logging.Formatter()
logging.basicConfig(format='%(asctime)s - \
                               %(name)s - \
                               %(levelname)s - \
                               %(message)s')


# Generic Python Packages
import numpy as np
import datetime as dt
import pandas as pd
import xarray as xr
import io
#from skewt.thermodynamics import TempK,DewPoint,MixR2VaporPress,GammaW,degCtoK, Rs_da, Cp_da,VaporPressure,MixRatio
from data_global import data_global
grav = 9.81

# this is just a generic input object
class generic_input(object):
    def __init__(self):
        self.init = True


# all units from all variables in CLASS(4GL) should be defined here!
units = {
         'h':'m',
         'theta':'K', 
         'q':'kg/kg',
         'cc': '-',
         'cveg': '-',
         'wg': 'm3 m-3',
         'w2': 'm3 m-3',
         #'wg': 'kg/kg',
         'Tsoil': 'K',
         'T2': 'K',
         'z0m': 'm',
         'alpha': '-',
         'LAI': '-',
         'dhdt':'m/h',
         'dthetadt':'K/h',
         'dqdt':'kg/kg/h',
         'BR': '-',
         'EF': '-',
         'advt_x': 'K/s',
         'advt_y': 'K/s',
}

class class4gl_input(object):
    """
    this is the class4gl_input. It extends the model_input, which is now
    assigned to self.pars. It now also includes initial profiles as pandas
    Dataframes:
        self.air_balloon: raw profile input for profile of u,v,theta,q (not used)
        self.air_ap : the same as self.air_balloonm, but for which a mixed
                      layer is fitted. Thi profile is used as input.
        self.air_ac : atmospheric circulation profiles for advection and
                      subsidence

    # FYI this was the way it was defined in an early version:
    #    class4gl_input = type('class4gl_input', (model_input,gl_input,gl_dia), dict(c='c'))
    """

    def __init__(self,set_pars_defaults=True,debug_level=logging.WARNING):

        """ set up logger (see: https://docs.python.org/2/howto/logging.html)
        """

        self.logger = logging.getLogger('class4gl_input')
        if debug_level is not None:
            self.logger.setLevel(debug_level)

            # # create logger
            # #self.logger = logging.getLogger('class4gl_input')
            # #self.logger.setLevel(debug_level)

            # # create console handler and set level to debug
            # ch = logging.StreamHandler()
            # ch.setLevel(debug_level)

            # # create formatter
            # formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
            # # add formatter to ch
            # ch.setFormatter(formatter)
     
            # # add ch to logger
            # self.logger.addHandler(ch)
            # # print("TESTTESTSETSTETETS")
            # # self.logger.warning("testsetsetsttets")
            # #stop

            # # """ end set up logger """



        # these are the standard model input single-value parameters for class
        self.pars = model_input()

        # diagnostic parameters of the initial profile
        self.diag = dict()

        # In this variable, we keep track of the different parameters from where it originates from. 
        self.sources = {}

        if set_pars_defaults:
            self.set_pars_defaults()

    def set_pars_defaults(self):

        """ 
        Create empty model_input and set up case
        """
        defaults = dict( 
        dt         = 60.    , # time step [s] 
        runtime    = 6*3600 ,  # total run time [s]
        
        # mixed-layer input
        sw_ml      = True   ,  # mixed-layer model switch
        sw_shearwe = False  ,  # shear growth mixed-layer switch
        sw_fixft   = False  ,  # Fix the free-troposphere switch
        h          = 200.   ,  # initial ABL height [m]
        Ps         = 101300.,  # surface pressure [Pa]
        divU       = 0.     ,  # horizontal large-scale divergence of wind [s-1]
        #fc         = 1.e-4  ,  # Coriolis parameter [m s-1]
        
        theta      = 288.   ,  # initial mixed-layer potential temperature [K]
        dtheta     = 1.     ,  # initial temperature jump at h [K]
        gammatheta = 0.006  ,  # free atmosphere potential temperature lapse rate [K m-1]
        advtheta   = 0.     ,  # advection of heat [K s-1]
        beta       = 0.2    ,  # entrainment ratio for virtual heat [-]
        wtheta     = 0.1    ,  # surface kinematic heat flux [K m s-1]
        
        q          = 0.008  ,  # initial mixed-layer specific humidity [kg kg-1]
        dq         = -0.001 ,  # initial specific humidity jump at h [kg kg-1]
        gammaq     = 0.     ,  # free atmosphere specific humidity lapse rate [kg kg-1 m-1]
        advq       = 0.     ,  # advection of moisture [kg kg-1 s-1]
        wq         = 0.1e-3 ,  # surface kinematic moisture flux [kg kg-1 m s-1]
        
        CO2        = 422.   ,  # initial mixed-layer CO2 [ppm]
        dCO2       = -44.   ,  # initial CO2 jump at h [ppm]
        gammaCO2   = 0.     ,  # free atmosphere CO2 lapse rate [ppm m-1]
        advCO2     = 0.     ,  # advection of CO2 [ppm s-1]
        wCO2       = 0.     ,  # surface kinematic CO2 flux [ppm m s-1]
        sw_wind    = True  ,  # prognostic wind switch
        u          = 0.     ,  # initial mixed-layer u-wind speed [m s-1]
        du         = 0.     ,  # initial u-wind jump at h [m s-1]
        gammau     = 0.     ,  # free atmosphere u-wind speed lapse rate [s-1]
        advu       = 0.     ,  # advection of u-wind [m s-2]
        v          = 0.0    , # initial mixed-layer u-wind speed [m s-1]
        dv         = 0.0    ,  # initial u-wind jump at h [m s-1]
        gammav     = 0.     ,  # free atmosphere v-wind speed lapse rate [s-1]
        advv       = 0.     ,  # advection of v-wind [m s-2]
        sw_sl      = True   , # surface layer switch
        ustar      = 0.3    ,  # surface friction velocity [m s-1]
        z0m        = 0.02   ,  # roughness length for momentum [m]
        z0h        = 0.02* 0.1 ,  # roughness length for scalars [m]
        sw_rad     = True   , # radiation switch
        lat        = 51.97  ,  # latitude [deg]
        lon        = -4.93  ,  # longitude [deg]
        doy        = 268.   ,  # day of the year [-]
        tstart     = 6.8    ,  # time of the day [h UTC]
        cc         = 0.0    ,  # cloud cover fraction [-]
        Q          = 400.   ,  # net radiation [W m-2] 
        dFz        = 0.     ,  # cloud top radiative divergence [W m-2] 
        ls_type    = 'js'   ,  # land-surface parameterization ('js' for Jarvis-Stewart or 'ags' for A-Gs)
        wg         = 0.21   ,  # volumetric water content top soil layer [m3 m-3]
        w2         = 0.21   ,  # volumetric water content deeper soil layer [m3 m-3]
        cveg       = 0.85   ,  # vegetation fraction [-]
        Tsoil      = 295.   ,  # temperature top soil layer [K]
        Ts         = 295.   ,    # initial surface temperature [K]
        T2         = 296.   ,  # temperature deeper soil layer [K]
        a          = 0.219  ,  # Clapp and Hornberger retention curve parameter a
        b          = 4.90   ,  # Clapp and Hornberger retention curve parameter b
        p          = 4.     ,  # Clapp and Hornberger retention curve parameter c
        CGsat      = 3.56e-6,  # saturated soil conductivity for heat
        wsat       = 0.472  ,  # saturated volumetric water content ECMWF config [-]
        wfc        = 0.323  ,  # volumetric water content field capacity [-]
        wwilt      = 0.171  ,  # volumetric water content wilting point [-]
        C1sat      = 0.132  ,  
        C2ref      = 1.8    ,
        LAI        = 2.     ,  # leaf area index [-]
        gD         = 0.0    ,  # correction factor transpiration for VPD [-]
        rsmin      = 110.   ,  # minimum resistance transpiration [s m-1]
        rssoilmin  = 50.    ,  # minimun resistance soil evaporation [s m-1]
        alpha      = 0.25   ,  # surface albedo [-]
        Wmax       = 0.0012 ,  # thickness of water layer on wet vegetation [m]
        Wl         = 0.0000 ,  # equivalent water layer depth for wet vegetation [m]
        Lambda     = 5.9    ,  # thermal diffusivity skin layer [-]
        c3c4       = 'c3'   ,  # Plant type ('c3' or 'c4')
        sw_cu      = False  ,  # Cumulus parameterization switch
        dz_h       = 150.   ,  # Transition layer thickness [m]
        cala       = None   ,  # soil heat conductivity [W/(K*m)]
        crhoc      = None   ,  # soil heat capacity  [J/K*m**3]
        sw_ls      = True   ,
        sw_ap      = True  ,   # switch that tells to initialize with fitted Air Profiles (eg., from balloon soundings) as input
        sw_ac      = None  ,   # switch that tells to use large-scale gridded Air Circulation (advection and subsindence) fields  as input from eg., ERA-INTERIM
        sw_lit     = False,
        )
        pars = model_input()
        for key in defaults:
            pars.__dict__[key] = defaults[key]
        
        self.update(source='defaults',pars=pars)
        
    def clear(self):
        """ this procudure clears the class4gl_input """

        for key in list(self.__dict__.keys()):
            del(self.__dict__[key])
        self.__init__()

    def dump(self,file):
        """ this procedure dumps the class4gl_input object into a yaml file
            
            Input: 
                - self.__dict__ (internal): the dictionary from which we read 
            Output:
                - file: All the parameters in self.__init__() are written to
                the yaml file, including pars, air_ap, sources etc.
        """
        file.write('---\n')
        index = file.tell()
        file.write('# CLASS4GL input; format version: 0.1\n')

        # write out the position of the current record
        yaml.dump({'index':index}, file, default_flow_style=False)

        # we do not include the none values
        for key,data in self.__dict__.items():
            #if ((type(data) == model_input) or (type(class4gl_input):
            if key == 'pars':

                pars = {'pars' : self.__dict__['pars'].__dict__}
                parsout = {}
                for key in pars.keys():
                    if pars[key] is not None:
                        parsout[key] = pars[key]

                yaml.dump(parsout, file, default_flow_style=False)
            elif type(data) == dict:
                if key == 'sources':
                    # in case of sources, we want to have a
                    # condensed list format as well, so we leave out
                    # 'default_flow_style=False'
                    yaml.dump({key : data}, file)
                else: 
                    yaml.dump({key : data}, file,
                              default_flow_style=False)
            elif type(data) == pd.DataFrame:
                # in case of dataframes (for profiles), we want to have a
                # condensed list format as well, so we leave out
                # 'default_flow_style=False'
                yaml.dump({key: data.to_dict(orient='list')},file)

                # # these are trials to get it into a more human-readable
                # fixed-width format, but it is too complex
                #stream = yaml.dump({key : False},width=100, default_flow_style=False)
                #file.write(stream)
                
                # workaround. I don't know how to put a table in a readable format by using yaml. So I do it manually here
                #file.write(key+': !!str |\n')
                #file.write(str(data)+'\n')
       
    def load_yaml_dict(self,yaml_dict,reset=True):
        """ this procedure loads class4gl_input data from a dictionary obtained from yaml
            
            Input: 
                - yaml_dict: the dictionary from which we read 
                - reset: reset data before reading        
            Output:
                - All the parameters in self, eg., (pars, air_ap, sources etc.,).
        """
        
        if reset:
            for key in list(self.__dict__.keys()):
                del(self.__dict__[key])
            self.__init__()

        for key,data in yaml_dict.items():
            if key == 'pars':
                self.__dict__[key] = model_input()
                self.__dict__[key].__dict__ = data
            elif key in ['air_ap','air_balloon','air_ac','air_ach']:
                self.__dict__[key] = pd.DataFrame(data)
            elif key == 'sources':
                self.__dict__[key] = data
            elif key == 'diag':
                self.__dict__[key] = data
            else: 
                warnings.warn("Key '"+key+"' may not be implemented.")
                self.__dict__[key] = data

    def update(self,source,**kwargs):
        """ this procedure is to make updates of input parameters and tracking
        of their source more convenient. It implements the assignment of
        parameter source/sensitivity experiment IDs ('eg.,
        'defaults', 'sounding balloon', any satellite information, climate
        models, sensitivity tests etc.). These are all stored in a convenient
        way with as class4gl_input.sources.  This way, the user can always consult with
        from where parameters data originates from.  
        
        Input:
            - source:    name of the underlying dataset
            - **kwargs: a dictionary of data input, for which the key values
            refer to the class4gl data type ('pars', 'air_ap', 'air_balloon', etc.) and
            the values is a again a dictionary/dataframe of datakeys/columns
            ('wg','PRES','datetime', ...) and datavalues (either single values,
            profiles ...), eg., 

                pars = {'wg': 0.007  , 'w2', 0.005}
                pars = {pd.Dataframe('PRES': [1005.,9523,...]  , 'THTA': [295.,
                                     300.,...]}
            
        Output:
            - self.__dict__[datatype] : object to which the parameters are
                                        assigned. They can be consulted with
                                        self.pars, self.profiles, etc.
                                        
            - self.sources[source] : It supplements the overview overview of
                                     data sources can be consulted with
                                     self.sources. The structure is as follows:
                                     as:
                self.sources = { 
                'wyoming': ['pars:datetime','air_balloon:PRES','air_ap:QABS', ...],
                'GLEAM' :  ['pars:wg','pars:w2', ...],
                 ...
                }
        
        """

        #print(source,kwargs)

        for key,data in kwargs.items():
            #print('update',key,data)

            #print(key)
            # if the key is not in class4gl_input object, then just add it. In
            # that case, the update procedures below will just overwrite it 
            if key not in self.__dict__:
                self.__dict__[key] = data


            

            #... we do an additional check to see whether there is a type
            # match. I not then raise a key error
            if (type(data) != type(self.__dict__[key]) \
                # we allow dict input for model_input pars
                and not ((key == 'pars') and (type(data) == dict) and \
                (type(self.__dict__[key]) == model_input))):

                raise TypeError('input key '+key+' is not of the same type as the one in the class4gl_object')


            # This variable keeps track of the added data that is supplemented
            # by the current source. We add this to class4gl_input.sources
            datakeys = []

            #... and we update the class4gl_input data, and this depends on the
            # data type

            if type(self.__dict__[key]) == pd.DataFrame:
                # If the data type is a dataframe, then we update the columns
                for column in list(data.columns):
                    #print(column)
                    self.__dict__[key][column] = data[column]
                    datakeys.append(column)
                    

            elif type(self.__dict__[key]) == model_input:
                # if the data type is a model_input (pars), then we update its internal
                # dictionary of parameters
                if type(data) == model_input:
                    self.__dict__[key].__dict__ = {**self.__dict__[key].__dict__, \
                                                   **data.__dict__}
                    datakeys = list(data.__dict__.keys())
                elif type(data) == dict:
                    datakeys = list(data.keys())
                    datavalues = list(data.values())
                    for idatavalue,datavalue in enumerate(datavalues):

                        # convert numpy to native python value types, so that
                        # we get clean output in the yaml file
                        if type(datavalue).__module__ == 'numpy':
                            datavalues[idatavalue] = datavalue.item()

                    self.__dict__[key].__dict__ = \
                            {**self.__dict__[key].__dict__, \
                            **dict(zip(datakeys, datavalues))}
                else:
                    raise TypeError('input key '+key+' is not of the same type\
                                    as the one in the class4gl_object')



            elif type(self.__dict__[key]) == dict:
                # if the data type is a dictionary, we update the
                # dictionary 
               # print('before update', self.__dict__[key] , data)
                self.__dict__[key] = {self.__dict__[key] , data}
               # print('after update',self.__dict__[key] )
                datakeys = list(data.keys())


            # if source entry is not existing yet, we add it
            if source not in self.sources.keys():
                self.sources[source] = []


            # self.logger.debug('updating section "'+\
            #                  key+' ('+' '.join(datakeys)+')'\
            #                  '" from source \
            #                  "'+source+'"')

            # Update the source dictionary: add the provided data keys to the
            # specified source list
            for datakey in datakeys:
                # At first, remove the occurences of the keys in the other
                # source lists
                for sourcekey,sourcelist in self.sources.items():
                    if key+':'+datakey in sourcelist:
                        self.sources[sourcekey].remove(key+':'+datakey)
                # Afterwards, add it to the current source list
                self.sources[source].append(key+':'+datakey)


        # # in case the datatype is a class4gl_input_pars, we update its keys
        # # according to **kwargs dictionary
        # if type(self.__dict__[datatype]) == class4gl_input_pars:
        #     # add the data parameters to the datatype object dictionary of the
        #     # datatype
        #     self.__dict__[datatype].__dict__ = {**self.__dict__[datatype].__dict__ ,
        #                                        **kwargs}
        # # in case, the datatype reflects a dataframe, we update the columns according
        # # to the *args list
        # elif type(self.__dict__[datatype]) == pd.DataFrame:
        #     for dataframe in args:
        #         for column in list(dataframe.columns):
        #             self.__dict__[datatype][column] = dataframe[column]
        

    def get_profile(self,IOBJ, *args, **argv):
        # if type(IOBJ) == wyoming:
        self.get_profile_wyoming(IOBJ,*args,**argv)
        # else:
        #     raise TypeError('Type '+str(type(IOBJ))+' is not supported')
        
    def get_profile_wyoming(self,wy_strm,air_ap_mode = 'b'):
        """ 
            Purpose: 
                This procedure assigns wyoming air profiles and parameters to the class4gl_input object.

            Input:
                1. wy_strm   = wyoming html (beautifulsoup) stream object. The
                function will take the profile at the stream's current
                position. 
                2. air_ap_mode: which air profile do we take? 
                    - b : best
                    - l : according to lower limit for the mixed-layer height
                            estimate
                    - u : according to upper limit for the mixed-layer height
                            estimate


            Output:
                1. all single-value parameters are stored in the
                   class4gl_input.pars object
                2. the souding profiles are stored in the in the
                   class4gl_input.air_balloon dataframe
                3. modified sounding profiles for which the mixed layer height
                   is fitted
                4. ...

        """


        # Raise an error in case the input stream is not the correct object
        # if type(wy_strm) is not wyoming:
        #    raise TypeError('Not a wyoming type input stream')

        # Let's tell the class_input object that it is a Wyoming fit type
        self.air_ap_type = 'wyoming'
        # ... and which mode of fitting we apply
        self.air_ap_mode = air_ap_mode

        """ Temporary variables used for output """
        # single value parameters derived from the sounding profile
        dpars = dict()
        # profile values
        air_balloon = pd.DataFrame()
        # fitted profile values
        air_ap = pd.DataFrame()
        
        string = wy_strm.current.find_next('pre').text
        string = string.split('\n')[:-1]
        string =  '\n'.join(string)
        
        columns = [ 'PRES', 'HGHT', 'TEMP', 'DWPT', 'RELH', 'MIXR', 'DRCT','SKNT' , 'THTA','THTE', 'THTV']             
        air_balloon = pd.read_fwf(io.StringIO(str(string)),widths=[7]*11,names=columns,skiprows=5,dtype=np.float,skipfooter=0)#.iloc[5:-1]
        #ONE_COLUMN = pd.read_table(io.StringIO(str(string)),sep=r"\s*",skiprows=[0,1,3,4])
        
        #string =  soup.pre.next_sibling.next_sibling
        
        string = wy_strm.current.find_next('pre').find_next('pre').text
        
        # this crazy long line just loads the sounding parameter table into parameters object (using amongst others the pandas internal engine to detect the right value types (int, float, np.Datetime64 etc.)).
        dpars = {**dpars,
                **pd.read_fwf(io.StringIO(str(string)),widths=[43,1,20],names=['descr','dummy','value']).iloc[1:-1].drop("dummy",1).set_index("descr").T.convert_objects(convert_numeric=True).iloc[0].to_dict()
               }
        
        # we get weird output when it's a numpy Timestamp, so we convert it to
        # pd.datetime type

        dpars['datetime'] = pytz.utc.localize(dt.datetime.strptime(dpars['Observation time'], "%y%m%d/%H%M"))
        dpars['STNID'] = dpars['Station number']

        # altitude above ground level
        air_balloon['z'] = air_balloon.HGHT -dpars['Station elevation']
        # absolute humidity in g/kg
        air_balloon['q']= (air_balloon.MIXR/1000.) \
                              / \
                             (air_balloon.MIXR/1000.+1.)
        # convert wind speed from knots to m/s
        air_balloon['WSPD'] = 0.51444 * air_balloon.SKNT
        angle_x = (90.-air_balloon.DRCT)/180.*np.pi # assuming that wind in direction of the south is 0 degrees.
        
        air_balloon['u'] = air_balloon.WSPD * np.sin(angle_x)
        air_balloon['v'] = air_balloon.WSPD * np.cos(angle_x)

        

        cp         = 1005.                 # specific heat of dry air [J kg-1 K-1]
        Rd         = 287.                  # gas constant for dry air [J kg-1 K-1]
        Rv         = 461.5                 # gas constant for moist air [J kg-1 K-1]

        air_balloon['R'] = (Rd*(1.-air_balloon.q) + Rv*air_balloon.q)
        air_balloon['p'] = air_balloon.PRES*100.


        # Therefore, determine the sounding that are valid for 'any' column 
        is_valid = ~np.isnan(air_balloon).any(axis=1) & (air_balloon.z >= 0)
        #is_valid = (air_balloon.z >= 0)
        # # this is an alternative pipe/numpy method
        # (~np.isnan(air_balloon).any(axis=1) & (air_balloon.z >= 0)).pipe(np.where)[0]
        valid_indices = air_balloon.index[is_valid].values
        #print(valid_indices)

        dpars['Ps'] = air_balloon.p.iloc[valid_indices[0]]

        air_balloon['t'] = air_balloon['TEMP']+273.15
        air_balloon['theta'] = (air_balloon.t) * \
                   (dpars['Ps']/(air_balloon.PRES*100.))**(air_balloon['R']/cp)
        air_balloon['thetav']   = air_balloon['theta']*(1. + 0.61 * air_balloon['q'])

        if len(valid_indices) > 0:
            #calculated mixed-layer height considering the critical Richardson number of the virtual temperature profile
            dpars['h'],dpars['h_u'],dpars['h_l'] = blh(air_balloon.z,air_balloon.thetav,air_balloon.WSPD)
            
            dpars['h_b'] = np.max((dpars['h'],10.))
            dpars['h_u'] = np.max((dpars['h_u'],10.)) #upper limit of mixed layer height
            dpars['h_l'] = np.max((dpars['h_l'],10.)) #low limit of mixed layer height
            dpars['h_e'] = np.abs( dpars['h_u'] - dpars['h_l']) # error of mixed-layer height
            
            # the final mixed-layer height that will be used by class. We round it
            # to 1 decimal so that we get a clean yaml output format
            dpars['h'] = np.round(dpars['h_'+air_ap_mode],1)
        else:
            dpars['h_u'] =np.nan
            dpars['h_l'] =np.nan
            dpars['h_e'] =np.nan
            dpars['h'] =np.nan


        if np.isnan(dpars['h']):
            dpars['Ps'] = np.nan

        if ~np.isnan(dpars['h']):
            # determine mixed-layer properties (moisture, potential temperature...) from profile
            
            # ... and those of the mixed layer
            is_valid_below_h = is_valid & (air_balloon.z < dpars['h'])
            valid_indices_below_h =  air_balloon.index[is_valid_below_h].values
            if len(valid_indices) > 1:
                if len(valid_indices_below_h) >= 3.:
                    ml_mean = air_balloon[is_valid_below_h].mean()
                else:
                    ml_mean = air_balloon.iloc[valid_indices[0]:valid_indices[1]].mean()
            elif len(valid_indices) == 1:
                ml_mean = (air_balloon.iloc[0:1]).mean()
            else:
                temp =  pd.DataFrame(air_balloon)
                temp.iloc[0] = np.nan
                ml_mean = temp
                       
            dpars['theta']= ml_mean.theta
            dpars['q']    = ml_mean.q
            dpars['u']    = ml_mean.u 
            dpars['v']    = ml_mean.v 
        else:
            dpars['theta'] = np.nan
            dpars['q'] = np.nan
            dpars['u'] = np.nan
            dpars['v'] = np.nan
            
        # First 3 data points of the mixed-layer fit. We create a empty head
        # first
        air_ap_head = air_balloon[0:0] #pd.DataFrame(columns = air_balloon.columns)
        # All other  data points above the mixed-layer fit
        air_ap_tail = air_balloon[air_balloon.z > dpars['h']]
        
        #calculate mixed-layer jump ( this should be larger than 0.1)
        
        air_ap_head['z'] = pd.Series(np.array([2.,dpars['h'],dpars['h']]))
        air_ap_head['HGHT'] = air_ap_head['z'] \
                                + \
                                np.round(dpars[ 'Station elevation'],1)
        
        # make a row object for defining the jump
        jump = air_ap_head.iloc[0] * np.nan
            
        if air_ap_tail.shape[0] > 1:

            # we originally used THTA, but that has another definition than the
            # variable theta that we need which should be the temperature that
            # one would have if brought to surface (NOT reference) pressure.
            for column in ['theta','q','u','v']:
               
               # initialize the profile head with the mixed-layer values
               air_ap_head[column] = ml_mean[column]
               # calculate jump values at mixed-layer height, which will be
               # added to the third datapoint of the profile head
               jump[column] = (air_ap_tail[column].iloc[1]\
                               -\
                               air_ap_tail[column].iloc[0])\
                              /\
                              (air_ap_tail.z.iloc[1]\
                               - air_ap_tail.z.iloc[0])\
                              *\
                              (dpars['h']- air_ap_tail.z.iloc[0])\
                              +\
                              air_ap_tail[column].iloc[0]\
                              -\
                              ml_mean[column] 
               if column == 'theta':
                  # for potential temperature, we need to set a lower limit to
                  # avoid the model to crash
                  jump.theta = np.max((0.1,jump.theta))
        
               air_ap_head[column][2] += jump[column]
        
        air_ap_head.WSPD = np.sqrt(air_ap_head.u**2 +air_ap_head.v**2)



        # make theta increase strong enough to avoid numerical
        # instability
        air_ap_tail_orig = pd.DataFrame(air_ap_tail)
        air_ap_tail = pd.DataFrame()
        #air_ap_tail = air_ap_tail.append(air_ap_tail_orig.iloc[0],ignore_index=True)
        #air_ap_tail = air_ap_tail.append(air_ap_tail_orig.iloc[0],ignore_index=True)
        theta_low = dpars['theta']
        z_low =     dpars['h']
        ibottom = 0
        for itop in range(0,len(air_ap_tail_orig)):
            theta_mean = air_ap_tail_orig.theta.iloc[ibottom:(itop+1)].mean()
            z_mean =     air_ap_tail_orig.z.iloc[ibottom:(itop+1)].mean()
            if (
                (z_mean > (z_low+10.)) and \
                (theta_mean > (theta_low+0.2) ) and \
                (((theta_mean - theta_low)/(z_mean - z_low)) > 0.0001)):

                air_ap_tail = air_ap_tail.append(air_ap_tail_orig.iloc[ibottom:(itop+1)].mean(),ignore_index=True)
                ibottom = itop+1
                theta_low = air_ap_tail.theta.iloc[-1]
                z_low =     air_ap_tail.z.iloc[-1]
            # elif  (itop > len(air_ap_tail_orig)-10):
            #     air_ap_tail = air_ap_tail.append(air_ap_tail_orig.iloc[itop],ignore_index=True)





        air_ap = \
            pd.concat((air_ap_head,air_ap_tail)).reset_index().drop(['index'],axis=1)
        
        # we copy the pressure at ground level from balloon sounding. The
        # pressure at mixed-layer height will be determined internally by class
        #print(air_ap['PRES'].iloc[0])

        rho        = 1.2                   # density of air [kg m-3]
        g          = 9.81                  # gravity acceleration [m s-2]

        air_ap['p'].iloc[0] =dpars['Ps'] 
        air_ap['p'].iloc[1] =(dpars['Ps'] - rho * g * dpars['h'])
        air_ap['p'].iloc[2] =(dpars['Ps'] - rho * g * dpars['h'] -0.1)

        
        dpars['lat'] = dpars['Station latitude']
        dpars['latitude'] = dpars['lat']
        
        # this is set to zero because we use local (sun) time as input (as if we were in Greenwhich)
        dpars['lon'] = 0.
        # this is the real longitude that will be used to extract ground data
        dpars['longitude'] = dpars['Station longitude']
        
        dpars['ldatetime'] = dpars['datetime'] \
                            + \
                            dt.timedelta(minutes=int(dpars['longitude']/360.*24.*60.))
        dpars['doy'] = dpars['datetime'].timetuple().tm_yday
        dpars['SolarAltitude'] = \
                                GetAltitude(\
                                    dpars['latitude'],\
                                    dpars['longitude'],\
                                    dpars['datetime']\
                                )
        dpars['SolarAzimuth'] =  GetAzimuth(\
                                    dpars['latitude'],\
                                    dpars['longitude'],\
                                    dpars['datetime']\
                                )
        dpars['lSunrise'], dpars['lSunset'] \
        =  GetSunriseSunset(dpars['latitude'],
                                         0.,
                                         dpars['ldatetime'])
        #print(dpars['lSunrise'])
        dpars['lSunrise'] = dpars['lSunrise']
        dpars['lSunset'] = dpars['lSunset']
        # dpars['lSunrise'] = pytz.utc.localize(dpars['lSunrise'])
        # dpars['lSunset'] = pytz.utc.localize(dpars['lSunset'])
        # This is the nearest datetime when the sun is up (for class)
        dpars['ldatetime_daylight'] = \
                                np.min(\
                                    (np.max(\
                                        (dpars['ldatetime'],\
                                         dpars['lSunrise'])\
                                     ),\
                                     dpars['lSunset']\
                                    )\
                                )
        # apply the same time shift for UTC datetime
        dpars['datetime_daylight'] = dpars['datetime'] \
                                    +\
                                    (dpars['ldatetime_daylight']\
                                     -\
                                     dpars['ldatetime'])
        
        dpars['doy'] = dpars['datetime'].timetuple().tm_yday

        # We set the starting time to the local sun time, since the model 
        # thinks we are always at the meridian (lon=0). This way the solar
        # radiation is calculated correctly.
        dpars['tstart'] = dpars['ldatetime_daylight'].hour \
                         + \
                         dpars['ldatetime_daylight'].minute/60.\
                         + \
                         dpars['ldatetime_daylight'].second/3600.
        

        # convert numpy types to native python data types. This provides
        # cleaner data IO with yaml:
        for key,value in dpars.items():
            if type(value).__module__ == 'numpy':
                dpars[key] = dpars[key].item()

        # # we make a pars object that is similar to the destination object
        # pars = model_input()
        # for key,value in dpars.items():
        #     pars.__dict__[key] = value


        # we round the columns to a specified decimal, so that we get a clean
        # output format for yaml
        decimals = {'p':0,'HGHT':1,'t':2,'DWPT':2,'RELH':2,'MIXR':2,\
                   'DRCT':2 ,'SKNT':2,   'theta':4,   'THTE':2,  'THTV':2,\
                   'z':2, 'q':5, 'WSPD':2, 'u':4,       'v':4}
# 
        for column,decimal in decimals.items():
            air_balloon[column] = air_balloon[column].round(decimal)
            air_ap[column] = air_ap[column].round(decimal)

        # in order to avoid warnings: the ABL values should have the same
        # rounding as the values profile.
        dpars['h'] = round(dpars['h'],decimals['z'])
        dpars['theta'] = round(dpars['theta'],decimals['theta'])
        dpars['q'] = round(dpars['q'],decimals['q'])
        dpars['u'] = round(dpars['u'],decimals['u'])
        dpars['v'] = round(dpars['v'],decimals['v'])

        self.update(source='wyoming',\
                    # pars=pars,
                    pars=dpars,\
                    air_balloon=air_balloon,\
                    air_ap=air_ap)

        
    def get_global_input(self, globaldata,only_keys=None,exclude_keys=None):
    
        """
        Purpose: This sets copies the parameters from the global datasets into the self (or similar object) 
                 according to the position (lat lon) and the class datetime and timespan
                 globaldata should be a globaldata multifile object
        
        Input: 
            - globaldata: this is the library object
            - only_keys: only extract specified keys
            - exclude_keys: do not inherit specified keys
        """
        classdatetime      = np.datetime64(self.pars.datetime_daylight)
        classdatetime_stop = np.datetime64(self.pars.datetime_daylight \
                                           + \
                                           dt.timedelta(seconds=self.pars.runtime)\
                                          )


        # # list of variables that we get from global ground data
        # self.ground_keys = ['fW', 'fB', 'fH', 'fTC', 'alpha', 'z0m', 'z0h', 
        #                 'wsat', 'Tsoil', 'cc', 'T2', 'wg', 'w2', 'wfc', 
        #                 'wwilt', 'DSMW', 'tex_coarse_values', 'tex_medium_values', 'tex_fine_values', 'code_values', 
        #                 'texture', 'itex', 'isoil', 'BR',
        #                 'b', 'cveg',
        #                 'C1sat', 
        #                 'C2ref', 'p', 'a',
        #                 ] #globaldata.datasets.keys():

        # # these are the required class4gl 3d atmospheric input which is not provided by the soundings
        # self.atm_keys = ['advtheta_x','advtheta_y','advu_x','advu_y','advv_x','advv_y','advq_x','advq_y','w','p']


        if type(globaldata) is not data_global:
            raise TypeError("Wrong type of input library") 

        # by default, we get all dataset keys
        keys = list(globaldata.datasets.keys())

        #print('keys orig', keys)

        # # In case there is surface pressure, we also calculate the half-level
        # # and full-level pressure fields
        # if ('sp' in keys):
        #     keys.append('pfull')
        #     keys.append('phalf')

        # If specified, we only take the keys that are in only_keys
        if only_keys is not None:
            cycle_keys = list(keys)
            for key in cycle_keys:
                if key not in only_keys:
                    keys.remove(key)

        #print('keys 1', keys)
                
        # If specified, we take out keys that are in exclude keys
        if exclude_keys is not None:
            for key in keys:
                if key in exclude_keys:
                    keys.remove(key)

        # We add LAI manually, because it is not listed in the datasets and
        #they its retreival is hard coded below based on LAIpixel and cveg
        if ('LAIpixel' in keys) and ('cveg' in keys):
            keys.append('LAI')

        # we set everything to nan first in the pars section (non-profile parameters
        # without lev argument), so that we can check afterwards whether the
        # data is well-fetched or not.

        for key in keys:
            if not ((key in globaldata.datasets) and \
                (globaldata.datasets[key].page is not None) and \
                ('lev' in globaldata.datasets[key].page[key].dims)):
                self.update(source='globaldata',pars={key:np.nan})
            # # we do not check profile input for now. We assume it is
            # # available
            #else:
            #    self.update(source='globaldata',air_ac=pd.DataFrame({key:list([np.nan])}))

        #print('keys 2', keys)
        print(keys)

        for key in keys:
            # If we find it, then we obtain the variables
            #print('key 0', key)
            if ((key in globaldata.datasets) and \
                (globaldata.datasets[key].page is not None)):

                #print('key 1', key)
                # check first whether the dataset has a height coordinate (3d space)
                if 'lev' in globaldata.datasets[key].page[key].dims:

                    # first, we browse to the correct file that has the current time
                    if 'time' in list(globaldata.datasets[key].page[key].dims):
                        globaldata.datasets[key].browse_page(time=classdatetime)

                    
                    if (globaldata.datasets[key].page is not None):
                        # find longitude and latitude coordinates
                        ilats = (np.abs(globaldata.datasets[key].page.lat -
                                        self.pars.latitude) < 0.5)
                        ilons = (np.abs(globaldata.datasets[key].page.lon -
                                        self.pars.longitude) < 0.5)
                        
                        # if we have a time dimension, then we look up the required timesteps during the class simulation
                        if 'time' in list(globaldata.datasets[key].page[key].dims):

                            DIST = np.abs((globaldata.datasets[key].page['time'].values - classdatetime))
                            
                            idatetime = np.where((DIST) == np.min(DIST))[0][0]
                            #print('idatetime',idatetime,globaldata.datasets[key].variables['time'].values[idatetime],classdatetime)
                            if key not in ['t','u','v','q']:
                                if ((globaldata.datasets[key].page.variables['time'].values[idatetime] < classdatetime) ):
                                    idatetime += 1
                            
                            DIST = np.abs((globaldata.datasets[key].page['time'].values - classdatetime_stop))
                            idatetimeend = np.where((DIST) == np.min(DIST))[0][0]
                            #print('idatetimeend',idatetimeend,globaldata.datasets[key].variables['time'].values[idatetime],classdatetimeend)
                            if ((globaldata.datasets[key].page.variables['time'].values[idatetimeend] > classdatetime_stop)):
                                idatetimeend -= 1
                            idatetime = np.min((idatetime,idatetimeend))
                            #for gleam, we take the previous day values

                            # in case of soil temperature or profile temperature, we take the exact
                            # timing (which is the morning)
                            if key in ['t','u','v','q']:
                                idatetimeend = idatetime
                            
                            itimes = range(idatetime,idatetimeend+1)
                            #print(key,'itimes',itimes)


                            # In case we didn't find any correct time, we take the
                            # closest one.
                            if len(itimes) == 0:


                                classdatetimemean = \
                                    np.datetime64(self.pars.datetime_daylight + \
                                    dt.timedelta(seconds=int(self.pars.runtime/2.)
                                                ))

                                dstimes = globaldata.datasets[key].page.time
                                time = dstimes.sel(time=classdatetimemean,method='nearest')
                                itimes = (globaldata.datasets[key].page.time ==
                                          time)
                                
                        else:
                            # we don't have a time coordinate so it doesn't matter
                            # what itimes is
                            itimes = 0

                        #multiplication by 1 is a trick to remove the array()-type in case of zero dimensions (single value).

                        # over which dimensions we take a mean:
                        dims = globaldata.datasets[key].page[key].dims
                        namesmean = list(dims)
                        namesmean.remove('lev')
                        idxmean = [dims.index(namemean) for namemean in namesmean]
                        
                        value = \
                        globaldata.datasets[key].page[key].isel(time=itimes,
                                                                lat=ilats,lon=ilons).mean(axis=tuple(idxmean)).values * 1.

                        # Ideally, source should be equal to the datakey of globaldata.library 
                        # or globaldata.datasets (eg., DSMW, IGBP-DIS, ERA-INTERIM etc.) 
                        #  but therefore the globaldata class requires a revision to make this work
                        self.update(source='globaldata',air_ac=pd.DataFrame({key:list(value)})) 

                else:
                    # this procedure is for reading the ground fields (2d space). 
                    # Actually, the code should be simplified to a similar fasion as the 3d procedure above and tested again.
                    if 'time' in list(globaldata.datasets[key].page[key].dims):
    
                       # first, we browse to the correct file
                       #print(key)
                       globaldata.datasets[key].browse_page(time=classdatetime)
    
                    if globaldata.datasets[key].page is not None:
                        DIST = \
                        np.abs((globaldata.datasets[key].page.variables['lat'].values\
                                - self.pars.latitude))
                        ilat = np.where((DIST) == np.min(DIST))[0][0]
                        DIST = \
                        np.abs((globaldata.datasets[key].page.variables['lon'].values\
                                - self.pars.longitude))
                        ilon = np.where((DIST) == np.min(DIST))[0][0]
                        
                        DIST = \
                        np.abs((globaldata.datasets[key].page.variables['lat'].values\
                                - (self.pars.latitude + 0.5)))
                        ilatmax = np.where((DIST) == np.min(DIST))[0][0]
                        if globaldata.datasets[key].page.variables['lat'].values[ilatmax] < globaldata.datasets[key].page.variables['lat'].values[ilat]:
                            ilatmax = ilat
                        
                        DIST = \
                        np.abs((globaldata.datasets[key].page.variables['lon'].values\
                                - (self.pars.longitude  + 0.5)))
                        ilonmax = np.where((DIST) == np.min(DIST))[0][0]
                        if globaldata.datasets[key].page.variables['lon'].values[ilonmax] < globaldata.datasets[key].page.variables['lon'].values[ilon]:
                            ilonmax = ilon
                        
                        DIST = \
                        np.abs((globaldata.datasets[key].page.lat.values\
                                - (self.pars.latitude - 0.5)))
                        ilatmin = np.where((DIST) == np.min(DIST))[0][0]
                        if globaldata.datasets[key].page.variables['lat'].values[ilatmin] > globaldata.datasets[key].page.variables['lat'].values[ilat]:
                            ilatmin = ilat
                        DIST = \
                        np.abs((globaldata.datasets[key].page.lon.values\
                                - (self.pars.longitude  - 0.5)))
                        ilonmin = np.where((DIST) == np.min(DIST))[0][0]
                        if globaldata.datasets[key].page.variables['lon'].values[ilonmin] > globaldata.datasets[key].page.variables['lon'].values[ilon]:
                            ilonmin = ilon        
                        
                        # for the koeppen climate classification we just take nearest
                        print(key)
                        if key == 'KGC':
                            ilatrange = range(ilat,ilat+1)
                            ilonrange = range(ilon,ilon+1)
                        else:
                            if ilatmin < ilatmax:
                                ilatrange = range(ilatmin,ilatmax+1)
                            else:
                                ilatrange = range(ilatmax,ilatmin+1)
                                
                            if ilonmin < ilonmax:
                                ilonrange = range(ilonmin,ilonmax+1)
                            else:
                                ilonrange = range(ilonmax,ilonmin+1)     
                            
                        if 'time' in list(globaldata.datasets[key].page[key].dims):
                            DIST = np.abs((globaldata.datasets[key].page['time'].values - classdatetime))
                            
                            idatetime = np.where((DIST) == np.min(DIST))[0][0]
                            #print('idatetime',idatetime,globaldata.datasets[key].variables['time'].values[idatetime],classdatetime)
                            if key not in ['Tsoil','T2']:
                                if ((globaldata.datasets[key].page.variables['time'].values[idatetime] < classdatetime) ):
                                    idatetime += 1
                            
                            classdatetimeend = np.datetime64(\
                                                             self.pars.datetime_daylight +\
                                                             dt.timedelta(seconds=self.pars.runtime)\
                                                            ) 
                            DIST = np.abs((globaldata.datasets[key].page['time'].values - classdatetimeend))
                            idatetimeend = np.where((DIST) == np.min(DIST))[0][0]
                            #print('idatetimeend',idatetimeend,globaldata.datasets[key].variables['time'].values[idatetime],classdatetimeend)
                            if ((globaldata.datasets[key].page.variables['time'].values[idatetimeend] > classdatetimeend)):
                                idatetimeend -= 1
                            idatetime = np.min((idatetime,idatetimeend))
                            #for gleam, we take the previous day values
                            if key in ['wg', 'w2']:
                                idatetime = idatetime - 1
                                idatetimeend = idatetimeend - 1

                            # in case of soil temperature, we take the exact
                            # timing (which is the morning)
                            if key in ['Tsoil','T2']:
                                idatetimeend = idatetime
                            
                            idts = range(idatetime,idatetimeend+1)
                            
                            count = 0
                            self.__dict__[key] = 0.
                            value = 0.
                            for iilat in ilatrange:
                                for iilon in ilonrange:
                                    for iidts in idts:
                                        value += np.mean(globaldata.datasets[key].page[key].isel(time=iidts,lat=iilat,lon=iilon,drop=True).values)
                                        count += 1
                            value = value/count
                            self.update(source='globaldata',pars={key:value.item()})
                                
                        else:
                                
                            count = 0
                            value = 0.
                            for iilat in ilatrange:
                                for iilon in ilonrange:
                                    value += np.mean(globaldata.datasets[key].page[key].isel(lat=iilat,lon=iilon,drop=True).values)
                                    count += 1
                            value = value/count                        

                            self.update(source='globaldata',pars={key:value.item()})

        if ('LAIpixel' in keys) and ('cveg' in keys):
            self.logger.debug('also update LAI based on LAIpixel and cveg') 
            # I suppose LAI pixel is already determined in the previous
            # procedure. Anyway...
            key = 'LAIpixel'

            if globaldata.datasets[key].page is not None:
                # first, we browse to the correct file that has the current time
                if 'time' in list(globaldata.datasets[key].page[key].dims):
                    globaldata.datasets[key].browse_page(time=classdatetime)
            
                DIST = \
                np.abs((globaldata.datasets[key].page.lat.values\
                        - self.pars.latitude))
                ilat = np.where((DIST) == np.min(DIST))[0][0]
                DIST = \
                np.abs((globaldata.datasets[key].page.lon.values\
                        - self.pars.longitude))
                ilon = np.where((DIST) == np.min(DIST))[0][0]
                 
                
                DIST = \
                np.abs((globaldata.datasets[key].page.lat.values\
                        - (self.pars.latitude + 0.5)))
                ilatmax = np.where((DIST) == np.min(DIST))[0][0]
                if globaldata.datasets[key].page.variables['lat'].values[ilatmax] < globaldata.datasets[key].page.variables['lat'].values[ilat]:
                    ilatmax = ilat
                
                DIST = \
                np.abs((globaldata.datasets[key].page.lon.values \
                        - (self.pars.longitude  + 0.5)))
                ilonmax = np.where((DIST) == np.min(DIST))[0][0]
                if globaldata.datasets[key].page.variables['lon'].values[ilonmax] < globaldata.datasets[key].page.variables['lon'].values[ilon]:
                    ilonmax = ilon
                
                DIST = \
                np.abs((globaldata.datasets[key].page.lat.values\
                        - (self.pars.latitude - 0.5)))
                ilatmin = np.where((DIST) == np.min(DIST))[0][0]
                if globaldata.datasets[key].page.variables['lat'].values[ilatmin] > globaldata.datasets[key].page.variables['lat'].values[ilat]:
                    ilatmin = ilat
                DIST = \
                np.abs((globaldata.datasets[key].page.lon.values\
                        - (self.pars.longitude  - 0.5)))
                ilonmin = np.where((DIST) == np.min(DIST))[0][0]
                if globaldata.datasets[key].page.variables['lon'].values[ilonmin] > globaldata.datasets[key].page.variables['lon'].values[ilon]:
                    ilonmin = ilon        
                DIST = np.abs((globaldata.datasets[key].page['time'].values - classdatetime))
                idatetime = np.where((DIST) == np.min(DIST))[0][0]
                
                
                if ilatmin < ilatmax:
                    ilatrange = range(ilatmin,ilatmax+1)
                else:
                    ilatrange = range(ilatmax,ilatmin+1)
                    
                if ilonmin < ilonmax:
                    ilonrange = range(ilonmin,ilonmax+1)
                else:
                    ilonrange = range(ilonmax,ilonmin+1)           
                
                #tarray_res = np.zeros(shape=globaldata.datasets[key]['time'].shape)
                LAIpixel = 0.
                count = 0
                for iilat in [ilat]: #ilatrange
                    for iilon in [ilon]: #ilonrange
                        LAIpixel += globaldata.datasets[key].page[key].isel(time = idatetime,lat=iilat,lon=iilon,drop=True).values
                        
                                        
                        # if np.isnan(tarray[idatetime]):
                        #     print("interpolating GIMMS LAIpixel nan value")
                        #     
                        #     mask = np.isnan(tarray)
                        #     
                        #     #replace each nan value with a interpolated value
                        #     if np.where(mask)[0].shape[0] < 0.25*mask.shape[0]:
                        #         tarray[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), tarray[~mask])
                        #         
                        #     else:
                        #         print("Warning. Could not interpolate GIMMS LAIpixel nan value")
                    
                        #         tarray *= np.nan 
                        
                        count += 1
                        #tarray_res += tarray
                LAIpixel = LAIpixel/count
                
                count = 0
                #tarray = globaldata.keys[dataset][key].isel({'lat':[ilat],'lon':[ilon]}).mean(dim=['lat','lon']).values
  
                self.update(source='globaldata',pars={'LAIpixel':np.float(LAIpixel)}) 
                #print('LAIpixel:',self.__dict__['LAIpixel'])
                #print('cveg:',self.__dict__['cveg'])
                
                # finally, we rescale the LAI according to the vegetation
                # fraction
                value = 0. 
                if ((self.pars.cveg is not None) and (self.pars.cveg > 0.1)):
                   value =self.pars.LAIpixel/self.pars.cveg
                else:
                    # in case of small vegetation fraction, we take just a standard 
                    # LAI value. It doesn't have a big influence anyway for
                    # small vegetation
                    value = 2.
                #print('LAI:',self.__dict__['LAI'])
                self.update(source='globaldata',pars={'LAI':value}) 


        # in case we have 'sp', we also calculate the 3d pressure fields at
        # full level and half level
        if ('sp' in keys) and ('sp' in self.pars.__dict__):
            pdAB = pd.read_fwf('/user/data/gent/gvo000/gvo00090/EXT/scripts/ECMWF/ecmwf_coeffs_L60_wrf.txt',header=None,names=['A','B'],index_col=0)  

            phalf,pfull =calc_air_ac_pres_L60(self.pars.sp,pdAB.A.values,pdAB.B.values)


            # # # STRANGE, THIS DOESN'T GIVE REALISTIC VALUES, IT NEEDS TO BE
            # # # CHECKED AGAIN SINCE THERE IS SIMILAR STRATEGY USED FOR 
            # # # CALCULATING THE ADVECTION PROFILES
            # # hydrostatic thickness of each model layer
            delpdgrav = -(phalf[:-1] - phalf[1:])/grav
            # # dz = rhodz/(R * T / pfull)


            # # subsidence multiplied by density. We calculate the subsidence of
            # # the in class itself
            # wrho = np.zeros_like(phalf)
            # wrho[-1] = 0. 

            # for ihlev in range(0,wrho.shape[0]-1):
            #     # subsidence multiplied by density is the integral of
            #     # divergences multiplied by the layer thicknessies
            #     wrho[ihlev] = ((self.air_ac['divU_x'][ihlev:] + \
            #                     self.air_ac['divU_y'][ihlev:]) * \
            #                    delpdgrav[ihlev:]).sum()


            
            self.update(source='globaldata',\
                        air_ac=pd.DataFrame({'p':list(pfull)}))
            self.update(source='globaldata',\
                        air_ach=pd.DataFrame({'p':list(phalf)}))
            self.update(source='globaldata',\
                        air_ac=pd.DataFrame({'delpdgrav':list(delpdgrav)}))
            # self.update(source='globaldata',\
            #             air_ach=pd.DataFrame({'wrho':list(wrho)}))


    # def get_idx_in_dataset(self,
    #                        globaldata,
    #                        latspan = 0.5):
    #                        lonspan = 0.5):
    #     """ 
    #     purpose:
    #         get the xarray indices that are representative between the starting and
    #         stopping time of the class simulations

    #     input:
    #         self: definition of the class input
    #         globaldata: book of class4gl global dataset
    #         key: key variable in the global dataset
    #         latspan: the span of the lat coordinate
    #         lonspan: the span of the lon coordinate


    #     output:
    #         itimes: time coordinates during of the class simulatios
    #         lats: 
    #         lons:
    #         """

    #     # first, we browse to the correct file that has the current time
    #     if 'time' in list(globaldata.datasets[key].page[key].dims):
    #         globaldata.datasets[key].browse_page(time=classdatetime)
    #     
    #     if (globaldata.datasets[key].page is not None):
    #         # find longitude and latitude coordinates
    #         ilats = (np.abs(globaldata.datasets[key].page.lat -
    #                         self.pars.latitude) < latspan)
    #         # In case we didn't find any latitude in the allowed range, we take the closest one.
    #         if len(ilats) == 0:
    #             ilats = np.where(\
    #                      globaldata.datasets[key].page.lat.isin(
    #                       globaldata.datasets[key].page.lat.sel(lat=self.pars.latitude)\
    #                      ))[0]
    #         ilons = (np.abs(globaldata.datasets[key].page.lon -
    #                         self.pars.longitude) < lonspan)
    #         # In case we didn't find any longitude in the allowed range, we take the closest one.
    #         if len(ilon) == 0:
    #             ilon = np.where(\
    #                      globaldata.datasets[key].page.lon.isin(
    #                       globaldata.datasets[key].page.lon.sel(lon=self.pars.longitude)\
    #                      ))[0]
    #         
    #         # if we have a time dimension, then we look up the required timesteps during the class simulation
    #         if 'time' in list(globaldata.datasets[key].page[key].dims):

    #             DIST = np.abs((globaldata.datasets[key].page['time'].values - classdatetime))
    #             
    #             idatetime = np.where((DIST) == np.min(DIST))[0][0]
    #             #print('idatetime',idatetime,globaldata.datasets[key].variables['time'].values[idatetime],classdatetime)
    #             if key not in ['t','u','v','q']:
    #                 if ((globaldata.datasets[key].page.variables['time'].values[idatetime] < classdatetime) ):
    #                     idatetime += 1
    #             
    #             DIST = np.abs((globaldata.datasets[key].page['time'].values - classdatetime_stop))
    #             idatetimeend = np.where((DIST) == np.min(DIST))[0][0]
    #             #print('idatetimeend',idatetimeend,globaldata.datasets[key].variables['time'].values[idatetime],classdatetimeend)
    #             if ((globaldata.datasets[key].page.variables['time'].values[idatetimeend] > classdatetime_stop)):
    #                 idatetimeend -= 1
    #             idatetime = np.min((idatetime,idatetimeend))
    #             #for gleam, we take the previous day values

    #             # in case of soil temperature, we take the exact
    #             # timing (which is the morning)
    #             if key in ['t','u','v','q']:
    #                 idatetimeend = idatetime
    #             
    #             itimes = range(idatetime,idatetimeend+1)
    #             #print(key,'itimes',itimes)


    #             # In case we didn't find any correct time, we take the
    #             # closest one.
    #             if len(itimes) == 0:


    #                 classdatetimemean = \
    #                     np.datetime64(self.pars.datetime_daylight + \
    #                     dt.timedelta(seconds=int(self.pars.runtime/2.)
    #                                 ))

    #                 dstimes = globaldata.datasets[key].page.time
    #                 time = dstimes.sel(time=classdatetimemean,method='nearest')
    #                 itimes = (globaldata.datasets[key].page.time ==
    #                           time)
    #                 
    #         else:
    #             # we don't have a time coordinate so it doesn't matter
    #             # what itimes is
    #             itimes = 0

    #         #multiplication by 1 is a trick to remove the array()-type in case of zero dimensions (single value).
    #       return itimes,ilats,ilons


    def query_source(self,var):
        """ 
        purpose:
            this procedure returns the name of the data source for a certain
            variable
        
        input:
            var: this should be in the format "section:variable", eg.,
            "pars:h", or "air_ac:theta"

        """

        for source,vars_in_source in self.sources.items():
            if var in vars_in_source:
                return source

    def check_source(self,source,check_only_sections=None,ignore_keys=[]):
        """ this procedure checks whether data of a specified source is valid.

        INPUT:
            source: the data source we want to check
            check_only_sections: a string or list with sections to be checked
        OUTPUT:
            returns True or False
        """

        # we set source ok to false as soon as we find a invalid input
        source_ok = True

        # convert to a single-item list in case of a string
        check_only_sections_def = (([check_only_sections]) if \
                                   type(check_only_sections) is str else \
                                    check_only_sections)
                                  
        if source not in self.sources.keys():
            self.logger.info('Source '+source+' does not exist')
            source_ok = False

        for sectiondatakey in self.sources[source]:                             
            section,datakey = sectiondatakey.split(':')                         
            if ((check_only_sections_def is None) or \
                (section in check_only_sections_def)):                          
                checkdatakeys = []
                if type(self.__dict__[section]) is pd.DataFrame:
                    checkdata = self.__dict__[section]
                elif type(self.__dict__[section]) is model_input:
                    checkdata = self.__dict__[section].__dict__

                if (datakey not in checkdata):                              
                    # self.logger.info('Expected key '+datakey+\
                    #                  ' is not in parameter input')                        
                    source_ok = False                                           
                elif (datakey not in ignore_keys) and \
                     ((checkdata[datakey] is None) or \
                     (pd.isnull(checkdata[datakey]) is True)):                    
        
                    # self.logger.info('Key value of "'+datakey+\
                    #                  '" is invalid: ('+ \
                    # str(self.__dict__[section].__dict__[datakey])+')')         
                    source_ok = False
                    self.logger.warning(datakey+' is invalid: '+ str(checkdata[datakey]))

        return source_ok

    def check_source_globaldata(self):
        """ this procedure checks whether all global parameter data is
        available, according to the keys in the self.sources"""

        source_globaldata_ok = True

        #self.get_values_air_input()

        # and now we can get the surface values
        #class_settings = class4gl_input()
        #class_settings.set_air_input(input_atm)
        
        # we only allow non-polar stations
        if not (self.pars.lat <= 60.):
            source_globaldata_ok = False
            self.logger.warning('cveg  is invalid: ('+str(self.pars.cveg)+')')
        
        # check lat and lon
        if (pd.isnull(self.pars.lat)) or (pd.isnull(self.pars.lon)):
            source_globaldata_ok = False
            self.logger.warning('lat  is invalid: ('+str(self.pars.lat)+')')
            self.logger.warning('or lon  is invalid: ('+str(self.pars.lon)+')')
        else:
            # we only check the ground parameter data (pars section). The 
            # profile data (air_ap section) are supposed to be valid in any 
            # case.
            source_ok = self.check_source(source='globaldata',\
                                          check_only_sections=['air_ac',\
                                                               'air_ap',\
                                                               'pars'],
                                         ignore_keys=[])
            if not source_ok:
                source_globaldata_ok = False
                self.logger.warning('something was wrong with the profiles')
        
            # Additional check: we exclude desert-like
            if ((self.pars.cveg is None) or pd.isnull(self.pars.cveg)):
                source_globaldata_ok = False
                self.logger.warning('cveg  is invalid: ('+str(self.pars.cveg)+')')
            if ((self.pars.LAI is None) or pd.isnull(self.pars.LAI)):
                source_globaldata_ok = False
                self.logger.warning('LAI  is invalid: ('+str(self.pars.LAI)+')')
            elif self.pars.cveg < 0.02:
                self.logger.warning('cveg  is too low: ('+str(self.pars.cveg)+')')
                source_globaldata_ok = False

        return source_globaldata_ok

    def mixed_layer_fit(self,air_ap,source,mode):
        """ 
            Purpose: 
                make a profile fit and write it to the air_ap section of the
                class4gl_input object (self).
            Input:
                air_ap: input profile


        """


        # Raise an error in case the input stream is not the correct object
        # if type(wy_strm) is not wyoming:
        #    raise TypeError('Not a wyoming type input stream')

        # Let's tell the class_input object that it is a Wyoming fit type
        self.air_ap_type = source+'_fit'
        # ... and which mode of fitting we apply
        self.air_ap_mode = mode


        # Therefore, determine the sounding that are valid for 'any' column 
        # is_valid = ~np.isnan(air_ap).any(axis=1) & (air_ap.z >= 0)
        
        if len(~np.isnan(air_ap).any(axis=1) & (air_ap.z >= 0)) == 0.:
            self.logger.warning('Warning, not all profile input is valid!  Please check input fields!', air_ap)


        is_valid = (air_ap.z >= 0)
        # # this is an alternative pipe/numpy method
        # (~np.isnan(air_ap).any(axis=1) & (air_ap.z >= 0)).pipe(np.where)[0]
        valid_indices = air_ap.index[is_valid].values
        print(valid_indices)


        hvalues = {}
        if len(valid_indices) > 0:
            #calculated mixed-layer height considering the critical Richardson number of the virtual temperature profile
            hvalues['h_b'] ,hvalues['h_u'],hvalues['h_l']  = blh(air_ap.z,air_ap.thetav,np.sqrt(air_ap.u**2. + air_ap.v**2.))
            
            hvalues['h_b']  = np.max((hvalues['h_b'] ,10.))
            hvalues['h_u']  = np.max((hvalues['h_u'] ,10.)) #upper limit of mixed layer height
            hvalues['h_l']  = np.max((hvalues['h_l'] ,10.)) #low limit of mixed layer height
            hvalues['h_e']  = np.abs( hvalues['h_u']  - hvalues['h_l'] ) # error of mixed-layer height
            
            # the final mixed-layer height that will be used by class. We round it
            # to 1 decimal so that we get a clean yaml output format
            hvalues['h']  = np.round(hvalues['h_'+mode],1)
        else:
            hvalues['h_u']  =np.nan
            hvalues['h_l']  =np.nan
            hvalues['h_e']  =np.nan
            hvalues['h']    =np.nan

        self.update(source='fit_from_'+source,pars=hvalues)

        if np.isnan(self.pars.h ):
            self.pars.Ps  = np.nan

        mlvalues = {}
        if ~np.isnan(self.pars.h ):
            # determine mixed-layer properties (moisture, potential temperature...) from profile
            
            # ... and those of the mixed layer
            is_valid_below_h = is_valid & (air_ap.z < self.pars.h)
            valid_indices_below_h =  air_ap.index[is_valid_below_h].values
            if len(valid_indices) > 1:
                if len(valid_indices_below_h) >= 3.:
                    ml_mean = air_ap[is_valid_below_h].mean()
                else:
                    ml_mean = air_ap.iloc[valid_indices[0]:valid_indices[1]].mean()
            elif len(valid_indices) == 1:
                ml_mean = (air_ap.iloc[0:1]).mean()
            else:
                temp =  pd.DataFrame(air_ap)
                temp.iloc[0] = np.nan
                ml_mean = temp
                       
            mlvalues['theta'] = ml_mean.theta
            mlvalues['q']     = ml_mean.q
            mlvalues['u']     = ml_mean.u 
            mlvalues['v']     = ml_mean.v 
        else:
            mlvalues['theta']  = np.nan
            mlvalues['q']  = np.nan
            mlvalues['u']  = np.nan
            mlvalues['v']  = np.nan
            

        self.update(source='fit_from_'+source,pars=mlvalues)


        # First 3 data points of the mixed-layer fit. We create a empty head
        # first
        air_ap_head = air_ap[0:0] #pd.DataFrame(columns = air_ap.columns)
        # All other  data points above the mixed-layer fit
        air_ap_tail = air_ap[air_ap.z > self.pars.h ]
        
        #calculate mixed-layer jump ( this should be larger than 0.1)
        
        air_ap_head['z'] = pd.Series(np.array([2.,self.pars.h ,self.pars.h ]))
        #air_ap_head['HGHT'] = air_ap_head['z'] \
        #                        + \
        #                        np.round(dpars[ 'Station elevation'],1)
        
        # make a row object for defining the jump
        jump = air_ap_head.iloc[0] * np.nan
            
        if air_ap_tail.shape[0] > 1:

            # we originally used THTA, but that has another definition than the
            # variable theta that we need which should be the temperature that
            # one would have if brought to surface (NOT reference) pressure.
            for column in ['theta','q','u','v']:
               
               # initialize the profile head with the mixed-layer values
               air_ap_head[column] = ml_mean[column]
               # calculate jump values at mixed-layer height, which will be
               # added to the third datapoint of the profile head
               jump[column] = (air_ap_tail[column].iloc[1]\
                               -\
                               air_ap_tail[column].iloc[0])\
                              /\
                              (air_ap_tail.z.iloc[1]\
                               - air_ap_tail.z.iloc[0])\
                              *\
                              (self.pars.h - air_ap_tail.z.iloc[0])\
                              +\
                              air_ap_tail[column].iloc[0]\
                              -\
                              ml_mean[column] 
               if column == 'theta':
                  # for potential temperature, we need to set a lower limit to
                  # avoid the model to crash
                  jump.theta = np.max((0.1,jump.theta))
        
               air_ap_head[column][2] += jump[column]
        
        air_ap_head.WSPD = np.sqrt(air_ap_head.u**2 +air_ap_head.v**2)



        # make theta increase strong enough to avoid numerical
        # instability
        air_ap_tail_orig = pd.DataFrame(air_ap_tail)
        air_ap_tail = pd.DataFrame()
        #air_ap_tail = air_ap_tail.append(air_ap_tail_orig.iloc[0],ignore_index=True)
        #air_ap_tail = air_ap_tail.append(air_ap_tail_orig.iloc[0],ignore_index=True)
        theta_low = self.pars.theta 
        z_low =     self.pars.h 
        ibottom = 0
        for itop in range(0,len(air_ap_tail_orig)):
            theta_mean = air_ap_tail_orig.theta.iloc[ibottom:(itop+1)].mean()
            z_mean =     air_ap_tail_orig.z.iloc[ibottom:(itop+1)].mean()
            if (
                (z_mean > (z_low+10.)) and \
                (theta_mean > (theta_low+0.2) ) and \
                (((theta_mean - theta_low)/(z_mean - z_low)) > 0.0001)):

                air_ap_tail = air_ap_tail.append(air_ap_tail_orig.iloc[ibottom:(itop+1)].mean(),ignore_index=True)
                ibottom = itop+1
                theta_low = air_ap_tail.theta.iloc[-1]
                z_low =     air_ap_tail.z.iloc[-1]
            # elif  (itop > len(air_ap_tail_orig)-10):
            #     air_ap_tail = air_ap_tail.append(air_ap_tail_orig.iloc[itop],ignore_index=True)





        air_ap = \
            pd.concat((air_ap_head,air_ap_tail)).reset_index().drop(['index'],axis=1)
        
        # we copy the pressure at ground level from balloon sounding. The
        # pressure at mixed-layer height will be determined internally by class
        #print(air_ap['PRES'].iloc[0])

        rho        = 1.2                   # density of air [kg m-3]
        g          = 9.81                  # gravity acceleration [m s-2]

        air_ap['p'].iloc[0] =self.pars.Ps  
        air_ap['p'].iloc[1] =(self.pars.Ps - rho * g * self.pars.h )
        air_ap['p'].iloc[2] =(self.pars.Ps - rho * g * self.pars.h -0.1)

        self.update(source='fit_from_'+source,air_ap=air_ap)



class c4gli_iterator():
    """ this iterator allows to loop through an entire yaml file and load class4gl_input sequentially 
    
        for information/documentation on creating such iterator classes, see: https://stackoverflow.com/questions/19151/build-a-basic-python-iterator
    """
    def __init__(self,file):
        # take file as IO stream
        self.file = file
        self.yaml_generator = yaml.load_all(file)
        self.current_dict = {}
        self.current_class4gl_input = class4gl_input()
        separator = self.file.readline() # this is just dummy
        self.header = file.readline()
        if self.header != '# CLASS4GL record; format version: 0.1\n':
            raise NotImplementedError("Wrong format version: '"+self.header+"'")
    def __iter__(self):
        return self
    def __next__(self):
        self.current_dict = self.yaml_generator.__next__()
        self.current_class4gl_input.load_yaml_dict(self.current_dict)
        return self.current_class4gl_input



#get_cape and lift_parcel are adapted from the SkewT package
    
class gl_dia(object):
    def get_lifted_index(self,timestep=-1):
        self.LI = get_lifted_index(self.input.Ps,self.out.T2m[timestep],self.out.q[timestep],self.p_pro,self.theta_pro,endp=50000.)
    
#from SkewT
#def get_lcl(startp,startt,startdp,nsteps=101):
#    from numpy import interp
#    #--------------------------------------------------------------------
#    # Lift a parcel dry adiabatically from startp to LCL.
#    # Init temp is startt in K, Init dew point is stwrtdp,
#    # pressure levels are in Pa    
#    #--------------------------------------------------------------------
#
#    assert startdp<=startt
#
#    if startdp==startt:
#        return np.array([startp]),np.array([startt]),np.array([startdp]),
#
#    # Pres=linspace(startp,60000.,nsteps)
#    Pres=np.logspace(np.log10(startp),np.log10(60000.),nsteps)
#
#    # Lift the dry parcel
#    T_dry=(startt)*(Pres/startp)**(Rs_da/Cp_da) 
#    # Mixing ratio isopleth
#    starte=VaporPressure(startdp)
#    startw=MixRatio(starte,startp)
#    e=Pres*startw/(.622+startw)
#    T_iso=243.5/(17.67/np.log(e/6.112)-1.) + degCtoK
#
#    # Solve for the intersection of these lines (LCL).
#    # interp requires the x argument (argument 2)
#    # to be ascending in order!
#    P_lcl=interp(0.,T_iso-T_dry,Pres)
#    T_lcl=interp(P_lcl,Pres[::-1],T_dry[::-1])
#
#    # # presdry=linspace(startp,P_lcl)
#    # presdry=logspace(log10(startp),log10(P_lcl),nsteps)
#
#    # tempdry=interp(presdry,Pres[::-1],T_dry[::-1])
#    # tempiso=interp(presdry,Pres[::-1],T_iso[::-1])
#
#    return P_lcl,T_lcl






#from class
def get_lcl(startp,startt,startqv):
        # Find lifting condensation level iteratively
    lcl = 20.
    RHlcl = 0.5
    
    itmax = 30
    it = 0
    while(((RHlcl <= 0.9999) or (RHlcl >= 1.0001)) and it<itmax):
        lcl    += (1.-RHlcl)*1000.
        p_lcl        = startp - 1.2 * 9.81 * lcl
        T_lcl        = startt - 9.81/1005. * lcl
        RHlcl        = startqv / qsat(T_lcl, p_lcl)
        it          += 1

    if(it == itmax):

        print("LCL calculation not converged!!")
        print("RHlcl = %f, zlcl=%f, theta=%f, q=%f"%(RHlcl, lcl,startt,startqv))

    return p_lcl,T_lcl

    
def lift_moist(startp,startt,endp,nsteps=501):
    #--------------------------------------------------------------------
    # Lift a parcel moist adiabatically from startp to endp.
    # Init temp is startt in K, pressure levels are in Pa    
    #--------------------------------------------------------------------
    
    # returns the end temp at endp (K)
    # preswet=linspace(startp,ptop,nsteps)
    preswet=np.logspace(np.log10(startp),np.log10(endp),nsteps)
    temp=startt
    tempwet=np.zeros(preswet.shape);tempwet[0]=startt
    for ii in range(preswet.shape[0]-1):
        delp=preswet[ii]-preswet[ii+1]
        temp=temp+delp*GammaW(temp,(preswet[ii]-delp/2.))
        tempwet[ii+1]=temp

    return tempwet[-1]




def lift_parcel(startp,startt,startqv,endp):
    """Do a lifted parcel analysis on the sounding data"""

    #temp = TempK(theta,pres)
    #startdp = DewPoint(MixR2VaporPress(startqv,startp))+degCtoK

    # Get Sub-LCL traces
    #P_lcl,T_lcl=get_lcl(startp,startt,startdp)
    
    P_lcl,T_lcl = get_lcl(startp,startt,startqv)

    # Get moist ascent traces
    endtemp=lift_moist(P_lcl,T_lcl,endp)

    return endtemp


def calc_air_ac_pres_L60(sp,A,B):
    """ returns the half-level and full level-pressure
    coordinates of an ECMWF grid based on the the surface pressure
    """
    phalf = A+B*sp
    # full levels (like phalf, but with length of dimension flev that is
    # one less than lev dimension
    pfull = np.zeros((phalf.shape[0]-1,))
    for ilev in range(pfull.shape[0]):
        pfull[ilev] = 0.5*(phalf[ilev]+phalf[ilev+1]) 
    return phalf,pfull



def get_lifted_index(startp,startt,startqv,pres,theta,endp=50000.):
    endthetaenv=np.interp(endp,pres[::-1],theta[::-1])
    endtempenv = TempK(endthetaenv,endp)
    endtemp = lift_parcel(startp,startt,startqv,endp)
    #print("endthetaenv",endthetaenv)
    #print("endtempenv ",endtempenv)
    #print(endtemp)
    return endtempenv - endtemp

def blh(HAGL,THTV,WSPD,RiBc = 0.31,RiBce = 0.08):
    """ Calculate mixed-layer height from temperature and wind speed profile

        Input:
            HAGL: height coordinates [m]
            THTV: virtual potential temperature profile [K]
            WSPD: wind speed profile [m/s]
            RIBc: critical Richardson Number. 
                According to Zhang et al., 2014 (GMD), it should  equal to 0.24
                for strongly stable boundary layers, 0.31 for weakly stable
                boundary layers, and 0.39 for unstable boundary layers. By
                default, it is set to the average of the three cases and an
                error RiBce value that comprises all values.

        Output:
            BLH: best-guess mixed-layer height
            BLHu: upper limit of mixed-layer height
            BLHl: lower limit of mixed-layer height

    """
    
    #initialize error BLH
    BLHe = 0.
    eps = 2.#security limit
    iTHTV_0 = np.where(~np.isnan(THTV))[0]
    if len(iTHTV_0) > 0:
        iTHTV_0 = iTHTV_0[0]
        THTV_0 = THTV[iTHTV_0]
    else:
        THTV_0 = np.nan

    RiB = 9.81/THTV_0 * ( THTV - THTV_0) * HAGL / np.clip(WSPD,a_min=0.1,a_max=None)**2.

    
    
    #RiB = 9.81/THTV_0 * ( THTV[i-1] +  (HGHT[i] - HGHT[i-1])/ - THTV_0) * HAGL / WSPD**2
    #RiB - RiBc = 0
    
    #best guess of BLH
    
    #print("RiB: ",RiB)
    #print("RiBc: ",RiBc)
    
    
    
    BLHi = np.where(RiB > RiBc)[0]
    if len(BLHi ) > 0:
        BLHi = BLHi[0]
        #print("BLHi: ",BLHi)
        BLH = (HAGL[BLHi] - HAGL[BLHi-1])/(RiB[BLHi] -RiB[BLHi-1]) * (RiBc - RiB[BLHi-1]) + HAGL[BLHi-1]
        
        # possible error is calculated as the difference height levels used for the interpolation
        BLHu = np.max([BLH,HAGL[BLHi]-eps])
        BLHl = np.min([BLH,HAGL[BLHi-1]+eps])
        # calculate an alternative BLH based on another critical Richardson number (RiBce):
        BLHi =np.where(RiB > RiBce)[0]
        if len(BLHi ) > 0:    
            BLHi = BLHi[0]
                
            BLHa = (HAGL[BLHi] - HAGL[BLHi-1])/(RiB[BLHi] -RiB[BLHi-1]) * (RiBc - RiB[BLHi-1]) + HAGL[BLHi-1]
            BLHu = np.max([BLHu,HAGL[BLHi]-eps])
            BLHl = np.min([BLHl,HAGL[BLHi-1]+eps])
            
            BLHu = np.max([BLHu,BLH + abs(BLH-BLHa)])
            BLHl = np.min([BLHl,BLH - abs(BLH-BLHa)])
        
        else:
            BLH,BLHu,BLHl = np.nan, np.nan,np.nan

    else:
        BLH,BLHu,BLHl = np.nan, np.nan,np.nan
        
    return BLH,BLHu,BLHl

class class4gl(model):
    """ the extension of the 'class model' class """

    def dump(self,file,include_input=False,timeseries_only=None):
        """ this procedure dumps the class4gl object into a yaml file
            
            Input: 
                - self.__dict__ (internal): the dictionary from which we read 
                - timeseries_only: for the timeseries output, dump only
                                   specific output variables
            Output:
                - file: All the parameters in self.__init__() are written to
                the yaml file, including pars, air_ap, sources etc.
        """

        if include_input:
            self.input_c4gl.dump(file)

        file.write('---\n')
        index = file.tell()
        file.write('# CLASS4GL input; format version: 0.1\n')

        # write out the position of the current record
        yaml.dump({'index':index}, file, default_flow_style=False)


        # we only copy those variables that are also in the input
        # for pars, these variables are in the model object
        dictpars = {}
        for key in self.__dict__.keys():
            # we omit the profile data and the timeseries in the out section,
            # since they are in pandas format for convenience
            # We include it hereafter as separate sections in the yaml file
            if key not in ['air_ach',\
                           'air_ac',\
                           'air_ap',\
                           'out',\
                           'input',\
                           'logger',\
                           'input_c4gl']:
                dictpars[key] = self.__dict__[key]

                # convert numpy types to native python data types. This
                # provides cleaner data IO with yaml:
                if type(dictpars[key]).__module__ == 'numpy':
                    dictpars[key] = dictpars[key].item() 




        dictout = {}
        dictoutlast = {}
        if timeseries_only == None:
            outvars = self.__dict__['out'].__dict__.keys()
        else:
            outvars = timeseries_only
        for key in outvars:
            dictout[key] = self.__dict__['out'].__dict__[key]
            dictoutlast[key] = dictout[key][-1]

            if type(dictoutlast[key]).__module__ == 'numpy':
                dictoutlast[key] = dictoutlast[key].item() 
            # convert numpy types to native python data types. This
            # provides cleaner data IO with yaml:
            if type(dictout[key]).__module__ == 'numpy':
                dictout[key] = [ a.item() for a in \
                                 self.__dict__['out'].__dict__[key]]
            #dictout[key] = list(dictout[key] )

        yaml.dump({'pars' : {**dictoutlast,**dictpars}},file)

        if ('sw_ac' in self.input.__dict__.keys()) \
           and (self.input.__dict__['sw_ac']):
            yaml.dump({'air_ac' : \
                       self.__dict__['air_ac'].to_dict(orient='list')},file)
            #yaml.dump({'air_ach' : \
            #           self.__dict__['air_ach'].to_dict(orient='list')},file)
        if ('sw_ap' in self.input.__dict__.keys()) and (self.input.__dict__['sw_ap']):
            #print('hello',self.air_ap.to_dict(orient='list'))
            yaml.dump({'air_ap' : \
                       self.__dict__['air_ap'].to_dict(orient='list')},file)

        yaml.dump({'out' : dictout},file)

