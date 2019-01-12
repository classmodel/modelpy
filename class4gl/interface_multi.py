import pandas as pd
import numpy as np
import datetime as dt
import os
import xarray as xr
import sys
from contextlib import suppress
from time import sleep


# sys.path.insert(0, '/user/data/gent/gvo000/gvo00090/D2D/software/CLASS/class4gl/')

from class4gl import class4gl_input, data_global,class4gl,units
from interface_functions import *
# from data_soundings import wyoming
import yaml
import glob
import pandas as pd
import json
import io
import subprocess
import pytz
from scipy.stats import mstats

from matplotlib.colors import LinearSegmentedColormap
cdictpres = {'blue': (\
                   (0.,    0.,  0.),
                   (0.25,  0.25, 0.25),
                   (0.5,  .70, 0.70),
                   (0.75, 1.0, 1.0),
                   (1,     1.,  1.),
                   ),
       'green': (\
                   (0. ,   0., 0.0),
                   (0.25,  0.50, 0.50),
                   (0.5,  .70, 0.70),
                   (0.75,  0.50, 0.50),
                   (1  ,    0,  0.),
                   ),
       'red':  (\
                  (0 ,  1.0, 1.0),
                  (0.25 ,  1.0, 1.0),
                   (0.5,  .70, 0.70),
                  (0.75 , 0.25, 0.25),
                  (1,    0., 0.),
                  )}

statsviewcmap = LinearSegmentedColormap('statsviewcmap', cdictpres)


os.system('module load Ruby')

class c4gl_interface_soundings(object):
    def __init__(self,path_exp,path_obs=None,globaldata=None,refetch_records=False,refetch_stations=True,inputkeys = ['cveg','wg','w2','cc','sp','wwilt','Tsoil','T2','z0m','alpha','LAI',],obs_filter=False,tendencies_revised=False):
        """ creates an interactive interface for analysing class4gl experiments

        INPUT:
            path_exp : path of the experiment output
            path_obs : path of the observations 
            globaldata: global data that is being shown on the map
            obs_filtering: extra data filter considering observation tendencies
                           beyond what the model can capture
            refetch_stations: do we need to build the list of the stations again?
        OUTPUT:
            the procedure returns an interface object with interactive plots

        """
        
        # set the ground
        self.globaldata = globaldata

 
        self.obs_filter = obs_filter
        print(self.obs_filter)
        self.tendencies_revised = tendencies_revised
        self.path_exp = path_exp
        self.path_obs = path_obs
        self.exp_files = glob.glob(self.path_exp+'/?????.yaml')

        # # get the list of stations
        # stationsfile = self.path_exp+'/stations_list.csv'
        # if (os.path.isfile(stationsfile)) and (not refetch_stations):
        #     stations = pd.read_csv(stationsfile)
        # else:
        #     stations = get_stations(self.path_exp)
        #     stations.to_csv(stationsfile)

        # stations = stations.set_index('STNID')

        self.frames = {}

        self.frames['stats'] = {}
        self.frames['worldmap'] = {}
                
        self.frames['profiles'] = {}
        self.frames['profiles'] = {}
        self.frames['profiles']['DT'] = None
        self.frames['profiles']['STNID'] = None

        #self.frames['worldmap']['stationsfile'] = stationsfile
        self.frames['worldmap']['stations'] = stations(self.path_exp, \
                                                       suffix='ini',\
                                                       refetch_stations=refetch_stations)

        # Initially, the stats frame inherets the values/iterators of
        # worldmap
        for key in self.frames['worldmap'].keys():
            self.frames['stats'][key] = self.frames['worldmap'][key]

        # get its records and load it into the stats frame
        self.frames['stats']['records_all_stations_ini'] =\
                        get_records(self.frames['stats']['stations'].table,\
                                           self.path_exp,\
                                           subset='ini',\
                                           refetch_records=refetch_records
                                           )
        # get its records and load it into the stats frame
        self.frames['stats']['records_all_stations_mod'] =\
                        get_records(self.frames['stats']['stations'].table,\
                                           self.path_exp,\
                                           subset='mod',\
                                           refetch_records=refetch_records
                                           )

        if self.path_obs is not None:
            # get its records and load it into the stats frame
            self.frames['stats']['records_all_stations_obs_afternoon'] =\
                            get_records(self.frames['stats']['stations'].table,\
                                               self.path_obs,\
                                               subset='afternoon',\
                                               refetch_records=refetch_records
                                               )

        self.frames['stats']['records_all_stations_mod'].index = \
            self.frames['stats']['records_all_stations_ini'].index 

        
        if len(self.frames['stats']['records_all_stations_ini']) ==0:
            raise ValueError('no class records found. Aborting')

        self.frames['stats']['records_all_stations_ini']['dates'] = \
            self.frames['stats']['records_all_stations_ini']['ldatetime'].dt.date

        if self.path_obs is not None:
            self.frames['stats']['records_all_stations_obs_afternoon']['dates'] = \
                self.frames['stats']['records_all_stations_obs_afternoon']['ldatetime'].dt.date

            self.frames['stats']['records_all_stations_obs_afternoon'].set_index(['STNID','dates'],inplace=True)


            ini_index_dates = self.frames['stats']['records_all_stations_ini'].set_index(['STNID','dates']).index

            self.frames['stats']['records_all_stations_obs_afternoon'] = \
                self.frames['stats']['records_all_stations_obs_afternoon'].loc[ini_index_dates]

            self.frames['stats']['records_all_stations_obs_afternoon'].index = \
                self.frames['stats']['records_all_stations_ini'].index 

            self.frames['stats']['viewkeys'] = ['h','theta','q']
            print('Calculating table statistics')

            if self.tendencies_revised:
                self.frames['stats']['records_all_stations_mod_stats'] = \
                        tendencies_rev(self.frames['stats']['records_all_stations_mod'],\
                                           self.frames['stats']['records_all_stations_ini'],\
                                           self.frames['stats']['viewkeys']\
                                  )
                self.frames['stats']['records_all_stations_obs_afternoon_stats'] = \
                        tendencies_rev(self.frames['stats']['records_all_stations_obs_afternoon'],\
                                           self.frames['stats']['records_all_stations_ini'],\
                                           self.frames['stats']['viewkeys']\
                                  )

            else:
                self.frames['stats']['records_all_stations_mod_stats'] = \
                        tendencies(self.frames['stats']['records_all_stations_mod'],\
                                   self.frames['stats']['records_all_stations_obs_afternoon'],\
                                   self.frames['stats']['records_all_stations_ini'],\
                                   self.frames['stats']['viewkeys']\
                                  )
                self.frames['stats']['records_all_stations_obs_afternoon_stats'] = \
                        tendencies(self.frames['stats']['records_all_stations_obs_afternoon'],\
                                   self.frames['stats']['records_all_stations_obs_afternoon'],\
                                   self.frames['stats']['records_all_stations_ini'],\
                                   self.frames['stats']['viewkeys']\
                                  )

        self.frames['stats']['inputkeys'] = inputkeys
        
        # self.frames['stats']['inputkeys'] = \
        #     [ key for key in \
        #       self.globaldata.datasets.keys() \
        #       if key in \
        #       list(self.frames['stats']['records_all_stations_obs'].columns)]


        # get units from the class4gl units database
        self.units = dict(units)
        # for those that don't have a definition yet, we just ask a question
        # mark
        for var in self.frames['stats']['inputkeys']:
            self.units[var] = '?'

        self.frames['worldmap']['inputkeys'] = self.frames['stats']['inputkeys'] 
        self.frames['stats']['records_all_stations_ini_pct'] = \
                  pct(self.frames['stats']['records_all_stations_ini'], \
                      columns = self.frames['stats']['inputkeys'])

        #     pd.DataFrame(columns = self.frames['stats']['viewkeys'])
        # for ikey,key in enumerate(self.frames['stats']['viewkeys']):
        #     mod['

        # 
        # 
        # \
        #        self.frames['stats']['records_all_stations_mod'], \



        # self.frames['stats']['records_all_stations_mod_stats_stdrel'] = \
        #        stdrel(mod = self.frames['stats']['records_all_stations_mod_stats'], \
        #               obs = self.frames['stats']['records_all_stations_obs_afternoon_stats'], \
        #               columns = [ 'd'+key+'dt' for key in \
        #                           self.frames['stats']['viewkeys']], \
        #              )

        # self.frames['stats']['records_all_stations_obs_afternoon_stats_stdrel'] = \
        #        stdrel(mod = self.frames['stats']['records_all_stations_ini'], \
        #               obs = self.frames['stats']['records_all_stations_ini'], \
        #               columns = self.frames['stats']['viewkeys'], \
        #              )

        

        if self.path_obs is not None:
            print('filtering pathological data')
            indextype = self.frames['stats']['records_all_stations_mod_stats'].index.names
            # some observational sounding still seem problematic, which needs to be
            # investigated. In the meantime, we filter them

            print('hello',self.obs_filter)
            print ((self.path_obs is not None) and (self.obs_filter))
            if ((self.path_obs is not None) and (self.obs_filter)) is True:
                print('hallohallo')
            if ((self.path_obs is not None) and (self.obs_filter)) is True:
                print('exclude exceptional observations')
                print('exclude unrealistic model output -> should be investigated!')
                valid = (\
                         (self.frames['stats']['records_all_stations_obs_afternoon_stats'].dthetadt >  0.250) & 
                         #(self.frames['stats']['records_all_stations_mod_stats'].dthetadt >  0.25000) & 
                         #(self.frames['stats']['records_all_stations_mod_stats'].dthetadt <  1.8000) & 
                         (self.frames['stats']['records_all_stations_obs_afternoon_stats'].dthetadt <  1.8000) & 
                         #(self.frames['stats']['records_all_stations_mod_stats'].dhdt >  50.0000) & 
                         (self.frames['stats']['records_all_stations_obs_afternoon_stats'].dhdt >  40.0000) & 
                         #(self.frames['stats']['records_all_stations_mod_stats'].dhdt <  350.) & 
                         (self.frames['stats']['records_all_stations_obs_afternoon_stats'].dhdt <  400.) & 
                         (self.frames['stats']['records_all_stations_obs_afternoon_stats'].dqdt >  -.00055) & 
                         #(self.frames['stats']['records_all_stations_mod_stats'].dqdt >  -.00055) & 
                         (self.frames['stats']['records_all_stations_obs_afternoon_stats'].dqdt <  .0003) & 

                         # filter 'extreme' model output -> should be investigated!
                         (self.frames['stats']['records_all_stations_mod_stats'].dqdt <  .0006) & 
                         (self.frames['stats']['records_all_stations_mod_stats'].dqdt >  -.0006) & 
                         (self.frames['stats']['records_all_stations_mod_stats'].dthetadt >  .2) & 
                         (self.frames['stats']['records_all_stations_mod_stats'].dthetadt <  2.) & 
                         # (self.frames['stats']['records_all_stations_mod_stats'].dqdt <  .0003) & 
                         # (self.frames['stats']['records_all_stations_ini'].KGC != 'Cwb') & 
                         # (self.frames['stats']['records_all_stations_ini'].KGC != 'Dfc') & 
                         ~np.isnan(self.frames['stats']['records_all_stations_mod_stats'].dthetadt) & 
                         ~np.isnan(self.frames['stats']['records_all_stations_obs_afternoon_stats'].dthetadt))

                for key in self.frames['stats'].keys():
                    if (type(self.frames['stats'][key]) == pd.DataFrame) and \
                       (self.frames['stats'][key].index.names == indextype):
                        self.frames['stats'][key] = self.frames['stats'][key][valid]
                print("WARNING WARNING!: "+ str(len(valid) - np.sum(valid))+' soundings are filtered')

        self.frames['stats']['records_all_stations_index'] = self.frames['stats']['records_all_stations_mod'].index


        print("filtering stations from interface that have no records")
        for STNID,station in self.frames['worldmap']['stations'].table.iterrows():
            if ((self.frames['stats']['records_all_stations_index'].get_level_values('STNID')\
                    == STNID).sum() == 0):
                print("dropping", STNID)
                self.frames['worldmap']['stations'].table = \
                        self.frames['worldmap']['stations'].table.drop(STNID)
                    
        self.frames['worldmap']['stations_iterator'] = stations_iterator(self.frames['worldmap']['stations']) 
        
        # TO TEST: should be removed, since it's is also done just below
        self.frames['stats']['stations_iterator'] = \
            self.frames['worldmap']['stations_iterator'] 

        self.frames['worldmap']['inputkey'] = self.frames['worldmap']['inputkeys'][0]
        self.frames['worldmap']['inputkey'] = self.frames['worldmap']['inputkey']
        self.next_station()

        # self.goto_datetime_worldmap(
        #     self.frames['profiles']['current_record_obs'].datetime.to_pydatetime(),
        #     'after')
    def sel_station(self,STNID=None,rownumber=None):

        if (STNID is not None) and (rownumber is not None):
            raise ValueError('Please provide either STNID or rownumber, not both.')

        if (STNID is None) and (rownumber is None):
            raise ValueError('Please provide either STNID or rownumber.')
            
        if STNID is not None:
            self.frames['worldmap']['STNID'],\
            self.frames['worldmap']['current_station'] \
             = self.frames['worldmap']['stations_iterator'].set_STNID(STNID)
            print(
            self.frames['worldmap']['STNID'],\
            self.frames['worldmap']['current_station'] \
            )
            self.update_station()
        elif rownumber is not None:
            self.frames['worldmap']['STNID'],\
            self.frames['worldmap']['current_station'] \
             = STNID,station = self.frames['worldmap']['stations_iterator'].set_row(rownumber)
            self.update_station()



    def next_station(self,event=None,jump=1):
        with suppress(StopIteration):
            self.frames['worldmap']['STNID'],\
            self.frames['worldmap']['current_station'] \
                = self.frames['worldmap']['stations_iterator'].__next__(jump)
            # self.frames['worldmap']['stations_iterator'].close()
            # del(self.frames['worldmap']['stations_iterator'])
            # self.frames['worldmap']['stations_iterator'] = \
            #                 selfself.frames['worldmap']['stations'].iterrows()
            # self.frames['worldmap']['STNID'],\
            # self.frames['worldmap']['current_station'] \
            #     = self.frames['worldmap']['stations_iterator'].__next__()

        self.update_station()

    def prev_station(self,event=None):
        self.next_station(jump = -1,event=event)
    def update_station(self):
        for key in ['STNID','current_station','stations_iterator']: 
            self.frames['stats'][key] = self.frames['worldmap'][key] 



        # generate index of the current station
        self.frames['stats']['records_current_station_index'] = \
            (self.frames['stats']['records_all_stations_index'].get_level_values('STNID')\
             == \
             self.frames['stats']['current_station'].name)

        # create the value table of the records of the current station
        tab_suffixes = \
                ['_mod','_ini','_ini_pct']
        if self.path_obs is not None:
            tab_suffixes=tab_suffixes+['_obs_afternoon','_mod_stats','_obs_afternoon_stats']

        for tab_suffix in tab_suffixes:
            self.frames['stats']['records_current_station'+tab_suffix] = \
                self.frames['stats']['records_all_stations'+tab_suffix].iloc[self.frames['stats']['records_current_station_index']]

        # go to first record of current station
        self.frames['stats']['records_iterator'] = \
                        records_iterator(self.frames['stats']['records_current_station_mod'])
        (self.frames['stats']['STNID'] , \
        self.frames['stats']['current_record_chunk'] , \
        self.frames['stats']['current_record_index']) , \
        self.frames['stats']['current_record_mod'] = \
                        self.frames['stats']['records_iterator'].__next__()

        for key in self.frames['stats'].keys():
            self.frames['profiles'][key] = self.frames['stats'][key]

        STNID = self.frames['profiles']['STNID']
        chunk = self.frames['profiles']['current_record_chunk']
        if 'current_station_file_ini' in self.frames['profiles'].keys():
            self.frames['profiles']['current_station_file_ini'].close()
        self.frames['profiles']['current_station_file_ini'] = \
            open(self.path_exp+'/'+format(STNID,"05d")+'_'+str(chunk)+'_ini.yaml','r')

        if 'current_station_file_mod' in self.frames['profiles'].keys():
            self.frames['profiles']['current_station_file_mod'].close()
        self.frames['profiles']['current_station_file_mod'] = \
            open(self.path_exp+'/'+format(STNID,"05d")+'_'+str(chunk)+'_mod.yaml','r')
        if 'current_station_file_afternoon' in self.frames['profiles'].keys():
            self.frames['profiles']['current_station_file_afternoon'].close()
        if self.path_obs is not None:
            self.frames['profiles']['current_station_file_afternoon'] = \
                open(self.path_obs+'/'+format(STNID,"05d")+'_afternoon.yaml','r')

        # for the profiles we make a distinct record iterator, so that the
        # stats iterator can move independently
        self.frames['profiles']['records_iterator'] = \
                        records_iterator(self.frames['profiles']['records_current_station_mod'])
        (self.frames['profiles']['STNID'] , \
        self.frames['profiles']['current_record_chunk'] , \
        self.frames['profiles']['current_record_index']) , \
        self.frames['profiles']['current_record_mod'] = \
                        self.frames['profiles']['records_iterator'].__next__()


        # for the profiles we make a distinct record iterator, so that the
        # stats iterator can move independently

        self.update_record()

    def next_record(self,event=None,jump=1):
        
        old_chunk =  self.frames['profiles']['current_record_chunk']

        with suppress(StopIteration):
            (self.frames['profiles']['STNID'] , \
            self.frames['profiles']['current_record_chunk'] , \
            self.frames['profiles']['current_record_index']) , \
            self.frames['profiles']['current_record_mod'] = \
                      self.frames['profiles']['records_iterator'].__next__(jump)
        # except (StopIteration):
        #     self.frames['profiles']['records_iterator'].close()
        #     del( self.frames['profiles']['records_iterator'])
        #     self.frames['profiles']['records_iterator'] = \
        #                 self.frames['profiles']['records_current_station_mod'].iterrows()
        #     (self.frames['profiles']['STNID'] , \
        #     self.frames['profiles']['current_record_index']) , \
        #     self.frames['profiles']['current_record_mod'] = \
        #                     self.frames['profiles']['records_iterator'].__next__()

        for key in self.frames['profiles'].keys():
            self.frames['stats'][key] = self.frames['profiles'][key]

        # chunk file has changed! So we need to open it!
        if self.frames['profiles']['current_record_chunk'] != old_chunk:

            STNID = self.frames['profiles']['STNID']
            chunk = self.frames['profiles']['current_record_chunk']



            if 'current_station_file_ini' in self.frames['profiles'].keys():
                self.frames['profiles']['current_station_file_ini'].close()
            self.frames['profiles']['current_station_file_ini'] = \
                open(self.path_exp+'/'+format(STNID,"05d")+'_'+str(chunk)+'_ini.yaml','r')

            if 'current_station_file_mod' in self.frames['profiles'].keys():
                self.frames['profiles']['current_station_file_mod'].close()
            self.frames['profiles']['current_station_file_mod'] = \
                open(self.path_exp+'/'+format(STNID,"05d")+'_'+str(chunk)+'_mod.yaml','r')

            if self.path_obs is not None:
                if 'current_station_file_afternoon' in self.frames['profiles'].keys():
                    self.frames['profiles']['current_station_file_afternoon'].close()
                self.frames['profiles']['current_station_file_afternoon'] = \
                    open(self.path_obs+'/'+format(STNID,"05d")+'_afternoon.yaml','r')

        self.update_record()

    def prev_record(self,event=None):
        self.next_record(jump=-1,event=event)

    def update_record(self):
        self.frames['profiles']['current_record_ini'] =  \
            self.frames['profiles']['records_current_station_ini'].loc[\
                  (self.frames['profiles']['STNID'] , \
                  self.frames['profiles']['current_record_chunk'],\
                  self.frames['profiles']['current_record_index'])]
        if self.path_obs is not None:
            self.frames['profiles']['current_record_obs_afternoon'] =  \
                self.frames['profiles']['records_current_station_obs_afternoon'].loc[\
                      (self.frames['profiles']['STNID'] , \
                      self.frames['profiles']['current_record_chunk'] , \
                      self.frames['profiles']['current_record_index'])]

            self.frames['profiles']['current_record_mod_stats'] = \
                    self.frames['profiles']['records_all_stations_mod_stats'].loc[(\
                        self.frames['profiles']['STNID'], \
                        self.frames['profiles']['current_record_chunk'], \
                        self.frames['profiles']['current_record_index'])]
            self.frames['profiles']['current_record_obs_afternoon_stats'] = \
                    self.frames['profiles']['records_all_stations_obs_afternoon_stats'].loc[(\
                        self.frames['profiles']['STNID'],\
                        self.frames['profiles']['current_record_chunk'],\
                        self.frames['profiles']['current_record_index'])]
        self.frames['profiles']['current_record_ini_pct'] = \
                self.frames['profiles']['records_all_stations_ini_pct'].loc[(\
                    self.frames['profiles']['STNID'],\
                    self.frames['profiles']['current_record_chunk'],\
                    self.frames['profiles']['current_record_index'])]

        for key in self.frames['profiles'].keys():
            self.frames['stats'][key] = self.frames['profiles'][key]
        # frame
        # note that the current station, record is the same as the stats frame for initialization

        # select first 
        #self.frames['profiles']['current_record_index'], \
        #self.frames['profiles']['record_yaml_mod'] = \
        #   get_record_yaml(self.frames['profiles']['current_station']['filename'],\
        #                   self.frames['stats']['current_record_index'])
        self.frames['profiles']['record_yaml_mod'] = \
           get_record_yaml(
               self.frames['profiles']['current_station_file_mod'], \
               self.frames['profiles']['current_record_mod'].index_start,
               self.frames['profiles']['current_record_mod'].index_end,
               mode='mod')
                                
        record_ini = self.frames['profiles']['records_all_stations_ini'].loc[
                       (self.frames['stats']['STNID'] , \
                        self.frames['stats']['current_record_chunk'] , \
                        self.frames['stats']['current_record_index'])]

        self.frames['profiles']['record_yaml_ini'] = \
           get_record_yaml(
               self.frames['profiles']['current_station_file_ini'], \
               record_ini.index_start,
               record_ini.index_end,
                mode='ini')

        if self.path_obs is not None:
            record_afternoon = self.frames['profiles']['records_all_stations_obs_afternoon'].loc[
                           (self.frames['stats']['STNID'] , \
                            self.frames['stats']['current_record_chunk'] , \
                            self.frames['stats']['current_record_index'])]

            self.frames['profiles']['record_yaml_obs_afternoon'] = \
               get_record_yaml(
                   self.frames['profiles']['current_station_file_afternoon'], \
                   record_afternoon.index_start,
                   record_afternoon.index_end,
                    mode='ini')


        key = self.frames['worldmap']['inputkey']
        # only redraw the map if the current world map has a time
        # dimension
        if (self.globaldata is not None) and ('time' in self.globaldata.datasets[key].page[key].dims):
            self.goto_datetime_worldmap(
                self.frames['profiles']['current_record_ini'].datetime.to_pydatetime(),
                'after')
            if "fig" in self.__dict__.keys():
                self.refresh_plot_interface(only=['stats_lightupdate',
                                                  'worldmap',
                                                  'profiles'])
        else:
            if "fig" in self.__dict__.keys():
                self.refresh_plot_interface(only=['stats_lightupdate',
                                                  'worldmap_stations',
                                                  'profiles'])

    def abline(self,slope, intercept,axis):
        """Plot a line from slope and intercept"""
        #axis = plt.gca()
        x_vals = np.array(axis.get_xlim())
        y_vals = intercept + slope * x_vals
        axis.plot(x_vals, y_vals, 'k--')

    def plot(self):
        import pylab as pl
        from matplotlib.widgets import Button
        import matplotlib.pyplot as plt
        import matplotlib as mpl
        '''
        Definition of the axes for the sounding table stats
        '''
        
        fig = pl.figure(figsize=(14,9))
        axes = {} #axes
        btns = {} #buttons

        # frames, which sets attributes for a group of axes, buttens, 
        if self.path_obs is not None:

            for ikey,key in enumerate(list(self.frames['stats']['records_all_stations_mod_stats'].columns)):
                label = 'stats_'+str(key)
                axes[label] = fig.add_subplot(\
                                len(self.frames['stats']['viewkeys']),\
                                5,\
                                5*ikey+1,label=label)
                # Actually, the axes should be a part of the frame!
                #self.frames['stats']['axes'] = axes[

                # pointer to the axes' point data
                axes[label].data = {}

                # pointer to the axes' color fields
                axes[label].fields = {}


        fig.tight_layout()
        fig.subplots_adjust(top=0.95,bottom=0.15,left=0.05,right=0.99,hspace=0.26,wspace=0.08)

        label ='stats_colorbar'
        axes[label] = fig.add_axes([0.025,0.06,0.18,0.025])
        axes[label].fields = {}

        from matplotlib.colors import LinearSegmentedColormap
        cdictpres = {'blue': (\
                           (0.,    0.,  0.),
                           (0.25,  0.25, 0.25),
                           (0.5,  .70, 0.70),
                           (0.75, 1.0, 1.0),
                           (1,     1.,  1.),
                           ),
               'green': (\
                           (0. ,   0., 0.0),
                           (0.25,  0.50, 0.50),
                           (0.5,  .70, 0.70),
                           (0.75,  0.50, 0.50),
                           (1  ,    0,  0.),
                           ),
               'red':  (\
                          (0 ,  1.0, 1.0),
                          (0.25 ,  1.0, 1.0),
                           (0.5,  .70, 0.70),
                          (0.75 , 0.25, 0.25),
                          (1,    0., 0.),
                          )}
        
        self.statsviewcmap = LinearSegmentedColormap('statsviewcmap', cdictpres)


        label = 'times'
               
        axes[label] = fig.add_axes([0.30,0.87,0.30,0.10]) #[*left*, *bottom*, *width*,    *height*]
        # add pointers to the data of the axes
        axes[label].data = {}
        # add pointers to color fields (for maps and colorbars) in the axes
        axes[label].fields = {}


        label = 'worldmap'
               
        axes[label] = fig.add_axes([0.25,0.48,0.40,0.35]) #[*left*, *bottom*, *width*,    *height*]
        # add pointers to the data of the axes
        axes[label].data = {}
        # add pointers to color fields (for maps and colorbars) in the axes
        axes[label].fields = {}
        axes[label].lat = None
        axes[label].lon = None

        label = 'worldmap_colorbar'
        axes[label] = fig.add_axes([0.25,0.44,0.40,0.05])
        axes[label].fields = {}

        # we make a overlying axes for the animations on the map, so that we don't need to redraw the whole map over and over again
        label = 'worldmap_stations'
        axes[label] = fig.add_axes([0.25,0.48,0.40001,0.350001]) #[*left*, *bottom*, *width*,    *height*]
        axes[label].data = {}

        fig.canvas.mpl_connect('pick_event', self.on_pick)
        fig.canvas.callbacks.connect('motion_notify_event', self.on_plot_hover)


        """ buttons definitions """
        button_height = 0.055
        button_hspace = 0.005
        button_width  = 0.095
        button_wspace = 0.005
        buttons_upper = 0.28
        buttons_left = 0.25

        button_types = ['dataset','datetime','level','station','record']
        
        for ibutton_type,button_type in enumerate(button_types):
            label='bprev'+button_type
            axes[label] = fig.add_axes([
                buttons_left,\
                buttons_upper-ibutton_type*(button_height+button_hspace),\
                button_width,\
                button_height\
                                                     ])
            btns[label] = Button(axes[label], 'Previous '+button_type)
            btns[label].on_clicked(getattr(self, 'prev_'+button_type))

            label='bnext'+button_type
            axes[label] = fig.add_axes([
                buttons_left+button_width+button_wspace,\
                buttons_upper-ibutton_type*(button_height+button_hspace),\
                button_width,\
                button_height\
                                                     ])
            btns[label] = Button(axes[label], 'Next '+button_type)
            btns[label].on_clicked(getattr(self, 'next_'+button_type))

        
        # label = 'bprev_dataset'
        # axes[label] = fig.add_axes([0.25,0.28,0.10,0.075])
        # btns[label] = Button(axes[label], 'Previous dataset')
        # btns[label].on_clicked(self.prev_dataset)

        # label = 'bnext_dataset'
        # axes[label] = fig.add_axes([0.35,0.28,0.10,0.075])
        # btns[label] = Button(axes[label], 'Next dataset')
        # btns[label].on_clicked(self.next_dataset)

        # label = 'bprev_datetime'
        # axes[label] = fig.add_axes([0.25,0.20,0.10,0.075])
        # btns[label] = Button(axes[label], 'Previous datetime')
        # btns[label].on_clicked(self.prev_datetime)

        # label = 'bnext_datetime'
        # axes[label] = fig.add_axes([0.35,0.20,0.10,0.075])
        # btns[label] = Button(axes[label], 'Next datetime')
        # btns[label].on_clicked(self.next_datetime)


        # label = 'bprev_station'
        # axes[label] = fig.add_axes([0.25,0.12,0.10,0.075])
        # btns[label] = Button(axes[label], 'Previous station')
        # btns[label].on_clicked(self.prev_station)

        # label = 'bnext_station'
        # axes[label] = fig.add_axes([0.35,0.12,0.10,0.075])
        # btns[label] = Button(axes[label], 'Next station')
        # btns[label].on_clicked(self.next_station)

        # label = 'bprev_record'
        # axes[label] = fig.add_axes([0.25,0.04,0.10,0.075])
        # btns[label] = Button(axes[label], 'Previous record')
        # btns[label].on_clicked(self.prev_record)

        # label = 'bnext_record'
        # axes[label] = fig.add_axes([0.35,0.04,0.10,0.075])
        # btns[label] = Button(axes[label], 'Next record')
        # btns[label].on_clicked(self.next_record)


        # self.nstatsview = nstatsview
        # self.statsviewcmap = statsviewcmap
        self.fig = fig
        self.axes = axes
        self.btns = btns
        self.tbox = {}
        # self.hover_active = False

        #self.tbox['loading'] = fig.text(0.30,0.01, " ",fontsize=10, 
        #                                transform=plt.gcf().transFigure)

        self.tbox['datetime'] =  fig.text(0.70, 0.96, " ", fontsize=10,
                                          transform=plt.gcf().transFigure)

        label = 'air_ap:theta'
        self.axes[label] = fig.add_axes([0.70,0.44,0.12,0.50], label=label)

        label = 'air_ap:q'
        self.axes[label] = fig.add_axes([0.86,0.44,0.12,0.50], label=label)

        label = 'out:h'
        self.axes[label] = fig.add_axes([0.50,0.27,0.22,0.09], label=label)

        label = 'out:theta'
        self.axes[label] = fig.add_axes([0.50,0.17,0.22,0.09], label=label)

        label = 'out:q'
        self.axes[label] = fig.add_axes([0.50,0.07,0.22,0.09], label=label)

        label = 'SEB'
        self.axes[label] = fig.add_axes([0.77,0.07,0.22,0.30], label=label)


        self.hover_active = False
        self.fig = fig
        self.fig.show()
        self.fig.canvas.draw()
        self.refresh_plot_interface()


    # def scan_stations(self):
    #     blabla
        


    # def get_records(current_file):
    #     records = pd.DataFrame()

    #     # initial position
    #     next_record_found = False
    #     while(not next_record_found):
    #         next_record_found = (current_file.readline() == '---\n')
    #     next_tell = current_file.tell() 
    #     end_of_file = (currentline == '') # an empty line means we are at the end

    #     while not end_of_file:
    #         current_tell = next_tell
    #         next_record_found = False
    #         current_file.seek(current_tell)
    #         while ( (not next_record_found) and (not end_of_file)):
    #             current_line = current_file.readline()
    #             next_record_found = (currentline == '---\n')
    #             end_of_file = (currentline == '') # an empty line means we are at the end

    #         # we store the position of the next record
    #         next_tell = current_file.tell() 
    #         
    #         # we get the current record. Unfortunately we need to reset the
    #         # yaml record generator first.
    #         current_yamlgen.close()
    #         current_yamlgen = yaml.load_all(current_file)
    #         current_file.seek(current_tell)
    #         current_record_mod = current_yamlgen.__next__()
    #     current_yamlgen.close()

    #     return records

       #      next_record_found = False
       #      while(not record):
       #          next_record_found = (self.current_file.readline() == '---\n')
       #      self.current_tell0 = self.current_file.tell() 

       #  

       #  next_record_found = False
       #  while(not next_record_found):
       #      next_record_found = (self.current_file.readline() == '---\n')
       #  self.current_tell0 = self.current_file.tell() 

       #  next_record_found = False
       #  while(not next_record_found):
       #      next_record_found = (self.current_file.readline() == '---\n')
       #  self.current_tell1 = self.current_file.tell() 


       #  self.current_yamlgen.close()
       #  self.current_yamlgen = yaml.load_all(self.current_file)
       #  self.current_file.seek(self.current_tell0)
       #  self.r0 = self.current_yamlgen.__next__()

       #  self.current_file.seek(self.current_tell1)
       #  next_record_found = False
       #  while ( (not next_record_found) and (not end_of_file):
       #      current_line = self.current_file.readline()
       #      next_record_found = (currentline == '---\n')
       #      end_of_file = (currentline == '') # an empty line means we are at the end

       #  self.current_tell2 = self.current_file.tell() 


       #  self.current_yamlgen.close()
       #  self.current_yamlgen = yaml.load_all(self.current_file)
       #  self.current_file.seek(self.current_tell1)
       #  self.r1 = self.current_yamlgen.__next__()

       #  self.current_file.seek(self.current_tell2)
       #  next_record_found = False
       #  while(not next_record_found):
       #      next_record_found = (self.current_file.readline() == '---\n')
       #  self.current_tell3 = self.current_file.tell() 

       #  self.current_yamlgen.close()
       #  self.current_yamlgen = yaml.load_all(self.current_file)
       #  self.current_file.seek(self.current_tell2)
       #  self.r2 = self.current_yamlgen.__next__()

       #  # go to position of next record in file
       #  self.current_file.seek(self.current_tell3)
       #  next_record_found = False
       #  while(not next_record_found):
       #      next_record_found = (self.current_file.readline() == '---\n')
       #  self.current_tell4 = self.current_file.tell() 

       #  self.current_yamlgen.close()
       #  self.current_yamlgen = yaml.load_all(self.current_file)
       #  self.current_file.seek(self.current_tell3)
       #  self.r3 = self.current_yamlgen.__next__()
 
       #  #self.update_tablestats(SOUNDINGS_TABLESTATS)

    def goto_datetime_worldmap(self,DT,shift=None):
        DT = np.datetime64(DT) #self.globaldata.datasets[self.axes['worldmap'].focus['key']].variables['time'].values[self.axes['worldmap'].focus['iDT']]
        if self.globaldata is not None:
            if 'time' in self.globaldata.datasets[self.frames['worldmap']['inputkey']].page.variables[self.frames['worldmap']['inputkey']].dims:
                self.globaldata.datasets[self.frames['worldmap']['inputkey']].browse_page(time=DT)
                DIST = np.abs((self.globaldata.datasets[self.frames['worldmap']['inputkey']].page.variables['time'].values - DT))
                self.frames['worldmap']['iDT'] = np.where((DIST) == np.min(DIST))[0][0]
                if ((shift == 'after') and (self.globaldata.datasets[self.frames['worldmap']['inputkey']].page.variables['time'].values[self.frames['worldmap']['iDT']] < DT)):
                    self.frames['worldmap']['iDT'] += 1
                elif ((shift == 'before') and (self.globaldata.datasets[self.frames['worldmap']['inputkey']].page.variables['time'].values[self.frames['worldmap']['iDT']] > DT)):
                    self.frames['worldmap']['iDT'] -= 1 
                # for gleam, we take the values of the previous day
                if self.frames['worldmap']['inputkey'] in ['wg','w2']:
                    self.frames['worldmap']['iDT'] -= 2 
                self.frames['worldmap']['DT'] = self.globaldata.datasets[self.frames['worldmap']['inputkey']].page.variables['time'].values[self.frames['worldmap']['iDT']]
            #else:
            #    self.frames['worldmap'].pop('DT')

    def next_datetime(self,event=None):
        if 'time' in self.globaldata.datasets[self.frames['worldmap']['inputkey']].page.variables[self.frames['worldmap']['inputkey']].dims:
            # for now we don't go to different files, so we cannot go to
            # another file 
            self.frames['worldmap']['iDT'] = (self.frames['worldmap']['iDT'] + 1) % len(self.globaldata.datasets[self.frames['worldmap']['inputkey']].page.variables['time'].values)
            self.frames['worldmap']['DT'] = self.globaldata.datasets[self.frames['worldmap']['inputkey']].page.variables['time'].values[self.frames['worldmap']['iDT']]
            if "fig" in self.__dict__.keys():
                self.refresh_plot_interface(only='worldmap') 

    def prev_datetime(self,event=None):
        if 'time' in self.globaldata.datasets[self.frames['worldmap']['inputkey']].page.variables[self.frames['worldmap']['inputkey']].dims:
            # for now we don't go to different files, so we cannot go to
            # another file 
            self.frames['worldmap']['iDT'] = (self.frames['worldmap']['iDT'] - 1) % len(self.globaldata.datasets[self.frames['worldmap']['inputkey']].page.variables['time'].values)
            self.frames['worldmap']['DT'] = self.globaldata.datasets[self.frames['worldmap']['inputkey']].page.variables['time'].values[self.frames['worldmap']['iDT']]
            if "fig" in self.__dict__.keys():
                self.refresh_plot_interface(only='worldmap') 

    def next_dataset(self,event=None):
        ikey = self.frames['worldmap']['inputkeys'].index(self.frames['worldmap']['inputkey'])
        ikey = (ikey + 1) % len(self.frames['worldmap']['inputkeys'])
        self.sel_dataset(self.frames['worldmap']['inputkeys'][ikey])
    def prev_dataset(self,event=None):
        ikey = self.frames['worldmap']['inputkeys'].index(self.frames['worldmap']['inputkey'])
        ikey = (ikey - 1) % len(self.frames['worldmap']['inputkeys'])
        self.sel_dataset(self.frames['worldmap']['inputkeys'][ikey])

    def sel_dataset(self,inputkey):
        self.frames['worldmap']['inputkey'] = inputkey
        self.frames['stats']['inputkey'] = self.frames['worldmap']['inputkey'] # this is used for showing the percentiles per station in color.
        self.goto_datetime_worldmap(
            self.frames['profiles']['current_record_ini'].datetime.to_pydatetime(),
            'after')# get nearest datetime of the current dataset to the profile

        print('seldata0')
        if 'level' not in self.frames['worldmap'].keys():
            levels = self.globaldata.datasets[self.frames['worldmap']['inputkey']].page['lev']
            self.frames['worldmap']['level'] = np.max(levels)
            print('seldata1')

            minlev = np.min(levels)
            maxlev = np.max(levels)
            curlev = self.frames['worldmap']['level']
            curlev = np.max([curlev,np.min(levels)])
            curlev = np.min([curlev,np.max(levels)])
            print('seldata2')

            self.frames['worldmap']['level'] = curlev
            print('seldata3')


        print('seldata4')
        self.sel_level(self.frames['worldmap']['level'])



    def sel_level(self,level):

        if 'lev' not in list(self.globaldata.datasets[self.frames['worldmap']['inputkey']].page.dims):
            raise ValueError('lev dimension not in dataset '+self.frames['worldmap']['inputkey'])

        print('seldata5')


        if level > (np.max(self.globaldata.datasets[self.frames['worldmap']['inputkey']].page['lev'])):
            raise ValueError('Level '+str(level)+' exceed those of the current dataset: '+str(self.globaldata.datasets[frames['worldmap']['inputkey']].page['lev']))
        if level < (np.min(self.globaldata.datasets[self.frames['worldmap']['inputkey']].page['lev'])):
            raise ValueError('Level '+str(level)+' is lower than those of the current dataset: '+str(self.globaldata.datasets[frames['worldmap']['inputkey']].page['lev']))
        print('seldata6')
        self.frames['worldmap']['level'] = level

        print(level)
        if "fig" in self.__dict__.keys():
            self.refresh_plot_interface(only=['worldmap','stats_lightupdate','stats_colorbar']) 

        print('seldata7')

    def next_level(self,event=None,jump=1):
        if 'lev' not in list(self.globaldata.datasets[self.frames['worldmap']['inputkey']].page.dims.keys()):
            raise ValueError('lev dimension not in dataset'+self.frames['worldmap']['inputkey'])
        levels =  self.globaldata.datasets[self.frames['worldmap']['inputkey']].page['lev']
        level = self.frames['worldmap']['level']
        level =  ((level + jump - min(levels)) % (max(levels)-min(levels))) + min(levels)
        self.sel_level(level)

    def prev_level(self,event=None):
        self.next_level(jump=-1)

        #self.frames['worldmap']['level'] = level: 
       
    # def prev_station(self,event=None):
    #     self.istation = (self.istation - 1) % self.stations.shape[0]
    #     self.update_station()




    #def update_datetime(self):
    #    if 'time' in self.globaldata.datasets[self.worldmapfocus['key']].variables[self.worldmapfocus['key']].dims:
    #    #if 'time' in list(dict(self.globaldata.datasets[self.worldmapfocus['key']].variables[self.worldmapfocus['key']].dims).keys()):
    #        #self.worldmapfocus['DT'] = self.globaldata.datasets[self.worldmapfocus['key']].variables['time'].values[self.worldmapfocus['iDT']]
    #        print(self.worldmapfocus['DT'])
    #        self.refresh_plot_interface(only='worldmap')

    def refresh_plot_interface(self,only=None,statsnewdata=True,**args):

        #print('r1')
        for argkey in args.keys():
            self.__dict__[arg] = args[argkey]

        axes = self.axes
        tbox = self.tbox
        frames = self.frames
        fig = self.fig
 
        if self.globaldata is not None:
            if (only is None) or ('worldmap' in only):
                globaldata = self.globaldata
                print('hello0')
                if 'time' in globaldata.datasets[frames['worldmap']['inputkey']].page.variables[frames['worldmap']['inputkey']].dims:
                    globaldata.datasets[frames['worldmap']['inputkey']].browse_page(time=frames['worldmap']['DT'])
                    datasetxr = globaldata.datasets[frames['worldmap']['inputkey']].page.isel(time = frames['worldmap']['iDT'])
                else:
                    datasetxr = globaldata.datasets[frames['worldmap']['inputkey']].page
                if 'lev' in datasetxr.dims:
                    datasetxr = datasetxr.isel(lev=self.frames['worldmap']['level'])
                keystotranspose = ['lat','lon']
                for key in dict(datasetxr.dims).keys():
                    if key not in keystotranspose:
                        keystotranspose.append(key)

                datasetxr = datasetxr.transpose(*keystotranspose)
                datasetxr = datasetxr.sortby('lat',ascending=False)
                print('hello1')

                lonleft = datasetxr['lon'].where(datasetxr.lon > 180.,drop=True) 
                lonleft = lonleft - 360.
                print('hello2')
                lonright = datasetxr['lon'].where(datasetxr.lon <= 180.,drop=True) 
                label = 'worldmap'
                axes[label].clear()
                axes[label].lon = xr.concat([lonleft,lonright],'lon').values
                axes[label].lat = np.sort(globaldata.datasets[frames['worldmap']['inputkey']].page.variables['lat'].values)[::-1] #sortby('lat',ascending=False).values
                print('hello3')

            #if 'axmap' not in self.__dict__ :
            #    self.axmap = self.fig.add_axes([0.39,0.5,0.34,0.5])
            #else:

            #stations = self.stations


            # self.gmap = Basemap(projection='kav7', lat_0 = 0, lon_0 =0,
            #     resolution = 'l', 
            # area_thresh = 0.1,
            #     llcrnrlon=-180., llcrnrlat=-90.0,
            #     urcrnrlon=180., urcrnrlat=90.0,ax=self.axmap)
            # 
            # self.gmap.drawcoastlines(color='white',linewidth=0.3)
            # self.gmap.drawcountries(color='white',linewidth=0.3)
            # #self.gmap.fillcontinents(color = 'gray')
            # self.gmap.drawmapboundary(color='white',linewidth=0.3)
            # # self.gmap.drawmeridians(np.arange(-180, 180+45, 60.),labels=[1,1,0,1])
            # # self.gmap.drawparallels(np.arange(-90, 90, 30.),labels=[1,0,0,0])
            # self.gmap.drawmeridians(np.arange(-180, 180+45, 60.),color='white',linewidth=0.3,labels=[0,0,0,0])
            # self.gmap.drawparallels(np.arange(-90, 90, 30.),color='white',linewidth=0.3,labels=[0,0,0,0])
            # #self.ax5.shadedrelief()

           #if 'time' in list(dict(self.datasets[self.axes['worldmap'].focus['key']].variables[self.axes['worldmap'].focus['key']].dims).keys()):


                fieldleft =  datasetxr[frames['worldmap']['inputkey']].where(datasetxr.lon > 180.,drop=True) 
                fieldright = datasetxr[frames['worldmap']['inputkey']].where(datasetxr.lon <= 180.,drop=True) 
                print('hello4')

                field =xr.concat([fieldleft,fieldright],'lon') #.sortby('lat',ascending=False).values

                #np.concatenate([viewframe.datasets['cc']['cc'].page.isel(time=0).where(viewframe.datasets['cc'].lon > 180).values,viewframe.datasets['cc']['cc'].isel(time=0).where(viewframe.datasets['cc'].lon <= 180).values],axis=1)
                axes[label].axis('off')
                print('hello5')

                from matplotlib import cm
                axes[label].fields[label] = axes[label].imshow(field[:,:],interpolation='none',cmap = cm.viridis )
                
                print('hello6')
                title=frames['worldmap']['inputkey']
                if globaldata is not None: 
                    if 'time' in globaldata.datasets[frames['worldmap']['inputkey']].page.variables[frames['worldmap']['inputkey']].dims:
                        title = title+' ['+pd.to_datetime(frames['worldmap']['DT']).strftime("%Y/%m/%d %H:%M") +'UTC]'
                axes[label].set_title(title)
                print('hello7')

                label ='worldmap_colorbar'
                axes[label].clear()
                axes[label].fields[label] = fig.colorbar(axes['worldmap'].fields['worldmap'],cax=axes[label],orientation='horizontal',label=frames['worldmap']['inputkey']+' ['+self.units[frames['worldmap']['inputkey']]+']')


                # lons, lats = np.meshgrid(axes[label].lon,axes[label].lat)
                # x,y = self.gmap(lons,lats)
                # #self.cont_map = self.axmap.contourf(x,y,field.T,cmap=gmapcm)
                # self.cont_map = self.axmap.pcolormesh(x,y,field.T,cmap=gmapcm)

        if (self.path_obs is not None) and \
           (self.frames['worldmap']['inputkey'] in self.frames['stats']['records_all_stations_ini_pct'].keys()) and \
           (self.path_obs is not None) and \
           ((only is None) or ('stats' in only) or ('stats_lightupdate' in only)):

            statskeys_out = list(self.frames['stats']['records_all_stations_mod_stats'].columns)
            store_xlim = {}
            store_ylim = {}
            for ikey, key in enumerate(statskeys_out):
                if (only is not None) and ('stats_lightupdate' in only):
                    store_xlim[key] = axes['stats_'+key].get_xlim()
                    store_ylim[key] = axes['stats_'+key].get_ylim()
                self.axes['stats_'+key].clear()    

            label = 'times'
            self.axes[label].clear()

            key = 'dthetadt'
            x = self.frames['stats']['records_all_stations_ini']['datetime']
            #print(x)
            y = self.frames['stats']['records_all_stations_obs_afternoon_stats'][key]
            #print(y)
            z = self.frames['stats']['records_all_stations_ini_pct'][self.frames['worldmap']['inputkey'] ]
            #print(z)

            alpha_cloud_pixels = 1./(1.+1./(0.15 * 10000. / len(self.frames['stats']['records_all_stations_mod'])))
            self.axes[label].data[label] = self.axes[label].scatter(x.values,
                                                                    y.values,
                                                                    c=z.values,
                                                                    cmap=self.statsviewcmap,
                                                                    s=2,
                                                                    vmin=0.,
                                                                    vmax=1.,
                                                                    alpha=alpha_cloud_pixels)

            
            x = self.frames['stats']['records_current_station_ini']['datetime']
            y = self.frames['stats']['records_current_station_obs_afternoon_stats'][key]
            z = self.frames['stats']['records_current_station_ini_pct'][self.frames['worldmap']['inputkey'] ]
            self.axes[label].data[label+'_current_station_hover'] = self.axes[label].scatter(x.values,y.values,c=z.values,cmap=self.statsviewcmap,s=5,picker=5,vmin=0.,vmax=1.,edgecolor='k',linewidth=0.3)


            x = self.frames['profiles']['records_current_station_ini']['datetime']
            y = self.frames['profiles']['records_current_station_obs_afternoon_stats'][key]
            z = self.frames['profiles']['records_current_station_ini_pct'][self.frames['worldmap']['inputkey'] ]

            self.axes[label].data[label+'_current_station'] = self.axes[label].scatter(x.values,y.values,c=z.values,cmap=self.statsviewcmap,s=20,picker=20,vmin=0.,vmax=1.,edgecolor='k',linewidth=0.8)

            self.axes[label].set_xlim((dt.datetime(1981,1,1),dt.datetime(2018,1,1)))
            self.axes[label].set_ylabel(key+ ' ['+self.units[key]+']')

            for ikey, key in enumerate(statskeys_out):

                # show data of all stations
                x = self.frames['stats']['records_all_stations_obs_afternoon_stats'][key]
                y = self.frames['stats']['records_all_stations_mod_stats'][key]
                z = self.frames['stats']['records_all_stations_ini_pct'][self.frames['worldmap']['inputkey'] ]
                qvalmax = x.quantile(0.999)
                qvalmin = x.quantile(0.001)
                print('applying extra filter over extreme values for plotting stats')
                selx = (x >= qvalmin) & (x < qvalmax)
                sely = (x >= qvalmin) & (x < qvalmax)
                x = x[selx & sely]
                y = y[selx & sely]
                z = z[selx & sely]
                self.axes['stats_'+key].data['stats_'+key] = \
                       self.axes['stats_'+key].scatter(x,y, c=z,\
                                cmap=self.statsviewcmap,\
                                s=3,picker=3,label=key,vmin=0.,vmax=1.,alpha=alpha_cloud_pixels)

                if len(x) > 1:
                    fit = np.polyfit(x, y, deg=1)
                    self.axes['stats_'+key].data['stats_'+key+'_fit'] = \
                         self.axes['stats_'+key].plot(x, fit[0] * x + fit[1], color='k',alpha=0.4,lw=4)

                x = self.frames['stats']['records_current_station_obs_afternoon_stats'][key]
                y = self.frames['stats']['records_current_station_mod_stats'][key]
                z = self.frames['stats']['records_current_station_ini_pct'][self.frames['worldmap']['inputkey'] ]
                self.axes['stats_'+key].data['stats_'+key+'_current_station_hover'] = \
                       self.axes['stats_'+key].scatter(x.values,y.values, c=z.values,\
                                cmap=self.statsviewcmap,\
                                s=10,picker=10,label=key,vmin=0.,vmax=1.,edgecolor='k',linewidth=0.3)

                x = self.frames['profiles']['records_current_station_obs_afternoon_stats'][key]
                y = self.frames['profiles']['records_current_station_mod_stats'][key]
                z = self.frames['profiles']['records_current_station_ini_pct'][self.frames['worldmap']['inputkey'] ]
                self.axes['stats_'+key].data['stats_'+key+'_current_station'] = \
                       self.axes['stats_'+key].scatter(x.values,y.values, c=z.values,\
                                cmap=self.statsviewcmap,\
                                s=20,picker=20,label=key,vmin=0.,vmax=1.,edgecolor='k',linewidth=0.8)

                if len(x) > 1:
                    fit = np.polyfit(x, y, deg=1)
                    self.axes['stats_'+key].data['stats_'+key+'_fit'] = \
                         self.axes['stats_'+key].plot(x, fit[0] * x + fit[1], color='k',alpha=0.8,lw=3)

                x = self.frames['stats']['current_record_obs_afternoon_stats'][key]
                y = self.frames['stats']['current_record_mod_stats'][key]
                z = self.frames['stats']['current_record_ini_pct'][self.frames['worldmap']['inputkey'] ]

                text = 'EXT: '+ format(x,'2.4f')+ ', MOD: ' + format(y,'2.4f')
                self.axes['stats_'+key].data['stats_'+key+'_current_record'] = \
                    axes['stats_'+key].annotate(text, \
                                               xy=(x,y),\
                                               xytext=(0.05,0.05),\
                                               textcoords='axes fraction',\
                                               bbox=dict(boxstyle="round",fc=self.statsviewcmap(z),edgecolor='black'),\
                                               color='white',\
                                               arrowprops=dict(arrowstyle="->",linewidth=1.1,color='black'))
                # self.axes['stats_'+key].data[key+'_current_record'] = \
                #        self.axes['stats_'+key].scatter(x,y, c=z,\
                #                 cmap=self.statsviewcmap,\
                #                 s=30,picker=15,label=key,vmin=0.,vmax=1.,edgecolor='k',linewidth=1.1)

                # axes['stats_'+key].set_title('relative deviation per station of '+ key)
                self.axes['stats_'+key].set_title(key+ ' ['+self.units[key]+']')
                # # highlight data for curent station
                # self.frames['stats']['records_all_stations_mod_stats'].iloc[self.frames['stats']['records_all_stations_index'].get_level_values('STNID') == self.frames['stats']['current_station'].name]

                #text = 'EXT: '+format(seltablestatsstdrel_statannotate[key+'_ext'],'2.4f')+ ', MOD: '+format(seltablestatsstdrel_statannotate[key+'_mod'],'2.4f')

                if ikey == len(statskeys_out)-1:
                    self.axes['stats_'+key].set_xlabel('external')
                    #axes[label].set_xlabel('ext: '+ key+' ['+statsunits[ikey]+']')
                axes['stats_'+key].set_ylabel('model')


                if (only is not None) and ('stats_lightupdate' in only):
                    self.axes['stats_'+key].set_xlim(*store_xlim[key])
                    self.axes['stats_'+key].set_ylim(*store_ylim[key])
                else:
                    limlow = np.min((axes['stats_'+key].get_xlim()[0],axes['stats_'+key].get_ylim()[0]))
                    limhigh = np.max((axes['stats_'+key].get_xlim()[1],axes['stats_'+key].get_ylim()[1]))
                    self.axes['stats_'+key].set_xlim(limlow,limhigh)
                    self.axes['stats_'+key].set_ylim(limlow,limhigh)
                self.abline(1,0,axis=self.axes['stats_'+key])

        if (only is None) or ('stats_colorbar' in only):
            label ='stats_colorbar'
            axes[label].clear()
            import matplotlib as mpl
            norm = mpl.colors.Normalize(vmin=0.,vmax=1.)
            self.axes[label].fields[label] = \
             mpl.colorbar.ColorbarBase(self.axes[label],\
                        orientation='horizontal',\
                        label="percentile of "+self.frames['worldmap']['inputkey'],
                        alpha=1.,
                                cmap=self.statsviewcmap,\
                                       norm=norm
                         )

        #print('r1')
        if (only is None) or ('worldmap' in only) or ('worldmap_stations' in only):
            #print('r2')
            label = 'worldmap_stations'
            axes[label].clear()
            
            stations = self.frames['worldmap']['stations'].table
            globaldata = self.globaldata
            
            key = label

            #print('r3')
            if (stations is not None):
                xlist = []
                ylist = []
                #print('r4')
                for iSTN,STN in frames['worldmap']['stations'].table.iterrows():
            #        x,y =self.gmap(STN['longitude'],STN['latitude'])
            #        self.gmap.plot(x,y, 'mo' if (self.STNID == STN['ID']) else 'ro' , markersize=1)
                    x,y = len(axes['worldmap'].lon)*(STN['longitude']- axes['worldmap'].lon[0])/(axes['worldmap'].lon[-1] - axes['worldmap'].lon[0]) ,len(axes['worldmap'].lat)*(STN['latitude']- axes['worldmap'].lat[0])/(axes['worldmap'].lat[-1] - axes['worldmap'].lat[0])
                    xlist.append(x)
                    ylist.append(y)
                #picker is needed to make it clickable (pick_event)
                axes[label].data[label] = axes[label].scatter(xlist,ylist,
                                                              c='r', s=15,
                                                              picker = 15,
                                                              label=key,
                                                              edgecolor='k',
                                                              linewidth=0.8)

            # cb.set_label('Wilting point [kg kg-3]')
                #print('r5')

                
            #     xseries = []
            #     yseries = []
            #     for iSTN,STN in stations.iterrows():
            # #        x,y =self.gmap(STN['longitude'],STN['latitude'])
            # #        self.gmap.plot(x,y, 'mo' if (self.STNID == STN['ID']) else 'ro' , markersize=1)
            #         x,y = len(axes[label].lon)*(STN['longitude_ext']- axes[label].lon[0])/(axes[label].lon[-1] - axes[label].lon[0])  ,len(axes[label].lat)*(STN['latitude_ext']- axes[label].axes[label].lat[0])/(axes[label].lat[-1] - axes[label].axes[label].lat[0])
            #         xseries.append(x)                    
            #         yseries.append(y)
            #         
            #         
            #     axes[label].data[label] = axes[label].scatter(xseries,yseries, c='r' , s=15, edgecolor='none',label=key)
                    
                if ('current_station' in frames['worldmap']):
                    #print('r5')
                    STN = frames['stats']['current_station']
                    STNID = frames['stats']['STNID']
                    #print('r5')

                    x,y = len(axes['worldmap'].lon)* \
                            (STN['longitude']- axes['worldmap'].lon[0])/(axes['worldmap'].lon[-1] - axes['worldmap'].lon[0]),\
                          len(axes['worldmap'].lat)* \
                            (STN['latitude']- axes['worldmap'].lat[0])/(axes['worldmap'].lat[-1] - axes['worldmap'].lat[0])
                    #print('r6')
                    #VAL = self.seltablestats[(self.seltablestats['STNID'] \
                    #                          == \
                    #                          self.frames['worldmap']['STNID'])\
                    #                         & \
                    #                         (self.seltablestats['DT'] \
                    #                          == self.axes['statsview0].focus['DT']) \
                    #                        ][self.axes['worldmap'].focus['key']+'_ext'].iloc[0]
                    #print('r7')
                    text = 'STNID: '+ format(STNID,'10.0f') + \
                            ', LAT: '+format(STN['latitude'],'3.3f')+ \
                            ', LON: '+format(STN['longitude'],'3.3f')+ \
                            ', #SOUNDINGS: '+str(self.frames['stats']['records_current_station_mod'].shape[0]) \

                            #+', VAL: '+format(VAL,'.3e')

                    axes[label].scatter(x,y, c='r', s=30,\
                                        edgecolor='k',picker=30,label=key,linewidth=1.1)
                    #print('r8')
            
                    #colorrange = list(axes[label].fields['worldmap'].get_clim())
                    #colorstation = (VAL-colorrange[0])/(colorrange[1]-colorrange[0])
                    #colorstation = max((min((1.,colorstation)),0.))
                    colorstation =0.2
                    from matplotlib import cm
                    axes[label].annotate(text,
                                         xy=(x,y),
                                         xytext=(0.05,0.05),
                                         textcoords='axes fraction', 
                                         bbox=dict(boxstyle="round",
                                         fc = cm.viridis(colorstation),edgecolor='black'),
                                         arrowprops=dict(arrowstyle="->",
                                                         linewidth=1.1,color='black'),
                                         color='white' if colorstation < 0.5 else 'black')
                    #print('r9')

                    # #pos = sc.get_offsets()[ind["ind"][0]]
                    # 
                    # axes[label.data[label+'statannotate'].xy = (seltablestatsstdrel_statannotate[key+'_ext'],seltablestatsstdrel_statannotate[key+'_mod'])
                    # text = 'STN: '+str(int(axes['statsview0'].focus['STNID']))+', DT: '+str(axes['statsview0'].focus['DT'])+', EXT: '+str(seltablestatsstdrel_statannotate[key+'_ext'])+', MOD: '+str(seltablestatsstdrel_statannotate[key+'_mod'])
                    # axes[label].data[label+'statannotate'].set_text(text)
                    #axes[label].data[label+'statannotate'].get_bbox_patch().set_facecolor(statsviewcmap(seltablestatspct_statannotate[cmapkey]))
                    # axes[label].data[label+'statannotate'].get_bbox_patch().set_alpha(0.4)
            #print('r9')
            axes[label].axis('off')
            axes[label].set_xlim(0,(len(axes['worldmap'].lon)))
            axes[label].set_ylim((len(axes['worldmap'].lat),0))
            #print('r10')

        if (only is None) or ('profiles' in only): 
            #print('r11')

            # # self.istation = np.where(self.stations['ID'] == STNID)[0][0]
            # # self.update_station(goto_first_sounding=False)
            # isounding = np.where(pd.DatetimeIndex(self.df_soundings_eval_pairs.datetime) == self.profilefocus['DT'])[0][0]
            # #self.isounding = (self.isounding - 1) % self.df_soundings_eval_pairs.shape[0]
            # self.morning_sounding = self.df_soundings_eval_pairs.loc[isounding]
            # self.evening_sounding = self.df_soundings.loc[self.morning_sounding['eval0']]

            label = 'air_ap:theta'
            axes[label].clear()

            tbox['datetime'].set_text(\
                self.frames['profiles']['record_yaml_ini'].pars.datetime.strftime("%Y/%m/%d %H:%M"))
                # +\
                # ' -> '+ \
                # self.frames['profiles']['record_yaml_obs_afternoon'].pars.datetime.strftime("%Y/%m/%d %H:%M"))
            
            
            
            
            #+self.evening_sounding.datetime.strftime("%Y/%m/%d %H:%M")+ "UTC")
            # 
            #print('r12')

            # #axes[label].set_title(self.morning_sounding.ldatetime.strftime("local time:  %H:%M")+' -> '+self.evening_sounding.ldatetime.strftime("%H:%M"))
            # #axes[label].set_title(self.morning_sounding.datetime.strftime("%Y/%m/%d %H:%M") + ' -> '+self.evening_sounding.datetime.strftime("%Y/%m/%d %H:%M"))
            # 
            #print(self.frames['profiles']['record_yaml_ini'].pars.h)
            #print(self.frames['profiles']['record_yaml_obs_afternoon'].pars.h)
            #print(self.frames['profiles']['record_yaml_mod'].out['h'].values[-1])
            hmax = np.nanmax([self.frames['profiles']['record_yaml_ini'].pars.h,\
                           self.frames['profiles']['record_yaml_mod'].out.h[-1]])
            if self.path_obs is not None:
                hmax = np.nanmax([hmax,self.frames['profiles']['record_yaml_obs_afternoon'].pars.h])


                zidxmax = int(np.where((self.frames['profiles']['record_yaml_ini'].air_balloon.z.values
                                    < 2.*hmax))[0][-1])+2
                zidxmax = np.min((zidxmax,len(self.frames['profiles']['record_yaml_ini'].air_balloon.z.values)))
                zco = range(zidxmax)

                axes[label].plot(self.frames['profiles']['record_yaml_ini'].air_balloon.theta.values[zco], \
                                 self.frames['profiles']['record_yaml_ini'].air_balloon.z.values[zco],"b*", \
                                 label="obs "+\
                                 self.frames['profiles']['record_yaml_ini'].pars.ldatetime.strftime("%H:%M")\
                                 +'LT')
            #print('r14')
            zidxmax = int(np.where((self.frames['profiles']['record_yaml_ini'].air_ap.z.values
                                < 2.*hmax))[0][-1])+2
            zidxmax = np.min((zidxmax,len(self.frames['profiles']['record_yaml_ini'].air_ap.z.values)))
            zco = range(zidxmax)

            axes[label].plot(self.frames['profiles']['record_yaml_ini'].air_ap.theta.values[zco], \
                             self.frames['profiles']['record_yaml_ini'].air_ap.z.values[zco],"b:", \
                             label="fit "+\
                             self.frames['profiles']['record_yaml_ini'].pars.ldatetime.strftime("%H:%M")\
                             +'LT')


            #print('r15')
            if self.path_obs is not None:
                zidxmax = int(np.where((self.frames['profiles']['record_yaml_obs_afternoon'].air_balloon.z.values
                                    < 2.*hmax))[0][-1])+2
                zidxmax = np.min((zidxmax,len(self.frames['profiles']['record_yaml_obs_afternoon'].air_balloon.z.values)))
                zco = range(zidxmax)

                              
                axes[label].plot(self.frames['profiles']['record_yaml_obs_afternoon'].air_balloon.theta.values[zco], \
                                 self.frames['profiles']['record_yaml_obs_afternoon'].air_balloon.z.values[zco],"r*", \
                                 label="obs "+\
                                 self.frames['profiles']['record_yaml_obs_afternoon'].pars.ldatetime.strftime("%H:%M")\
                                 +'LT')

                #print('r16')

                zidxmax = int(np.where((self.frames['profiles']['record_yaml_obs_afternoon'].air_ap.z.values < 2.*hmax))[0][-1])+2
                zidxmax = np.min((zidxmax,len(self.frames['profiles']['record_yaml_obs_afternoon'].air_ap.z.values)))
                zco = range(zidxmax)

                axes[label].plot(self.frames['profiles']['record_yaml_obs_afternoon'].air_ap.theta.values[zco], \
                                 self.frames['profiles']['record_yaml_obs_afternoon'].air_ap.z.values[zco],"r:", \
                                 label="fit "+\
                                 self.frames['profiles']['record_yaml_obs_afternoon'].pars.ldatetime.strftime("%H:%M")\
                                 +'LT')

            #print('r17')
            #print(self.frames['profiles']['record_yaml_mod'].air_ap.z)
            #print(hmax)
            valid_mod = len(self.frames['profiles']['record_yaml_mod'].air_ap.z)>= 4
            if valid_mod:

                zidxmax = int(np.where((self.frames['profiles']['record_yaml_mod'].air_ap.z.values < 2.*hmax))[0][-1])+2
                zidxmax = np.min((zidxmax,len(self.frames['profiles']['record_yaml_mod'].air_ap.z.values)))
                zco = range(zidxmax)

                axes[label].plot(self.frames['profiles']['record_yaml_mod'].air_ap.theta.values[zco], \
                                 self.frames['profiles']['record_yaml_mod'].air_ap.z.values[zco],"r-", \
                                 label="mod "+\
                                 (self.frames['profiles']['record_yaml_ini'].pars.ldatetime
                                 +dt.timedelta(seconds=self.frames['profiles']['record_yaml_ini'].pars.runtime)).strftime("%H:%M")\
                                 +'LT')

            #print('r18')
            axes[label].legend(prop={'family':'monospace'},loc='upper left')
            axes[label].set_ylabel('height [m]')
            axes[label].set_xlabel('theta [K]')

            label = 'air_ap:q'
            axes[label].clear()

            tbox['datetime'].set_text(\
                 (self.frames['profiles']['record_yaml_ini'].pars.datetime_daylight+\
                  dt.timedelta(seconds=self.frames['profiles']['record_yaml_ini'].pars.runtime)).strftime("%Y/%m/%d %H:%M")
                )
            
            #+self.evening_sounding.datetime.strftime("%Y/%m/%d %H:%M")+ "UTC")
            # 

            #print('r19')
            # #axes[label].set_title(self.morning_sounding.ldatetime.strftime("local time:  %H:%M")+' -> '+self.evening_sounding.ldatetime.strftime("%H:%M"))
            # #axes[label].set_title(self.morning_sounding.datetime.strftime("%Y/%m/%d %H:%M") + ' -> '+self.evening_sounding.datetime.strftime("%Y/%m/%d %H:%M"))
            # 
            if valid_mod:
                hmax = np.nanmax([self.frames['profiles']['record_yaml_ini'].pars.h,\
                               self.frames['profiles']['record_yaml_mod'].out.h[-1]])
            else:
                hmax = self.frames['profiles']['record_yaml_ini'].pars.h

            if self.path_obs is not None:
                hmax = np.nanmax([hmax,self.frames['profiles']['record_yaml_obs_afternoon'].pars.h])
            # 
            #print('r20')

                zidxmax = int(np.where((self.frames['profiles']['record_yaml_ini'].air_balloon.z.values < 2.*hmax))[0][-1])+2
                zidxmax = np.min((zidxmax,len(self.frames['profiles']['record_yaml_ini'].air_balloon.z.values)))
                zco = range(zidxmax)

                axes[label].plot(self.frames['profiles']['record_yaml_ini'].air_balloon.q.values[zco], \
                                 self.frames['profiles']['record_yaml_ini'].air_balloon.z.values[zco],"b*", \
                                 label="obs "+\
                                 self.frames['profiles']['record_yaml_ini'].pars.ldatetime.strftime("%H:%M")\
                                 +'LT')
                #print('r21')


            zidxmax = int(np.where((self.frames['profiles']['record_yaml_ini'].air_ap.z.values < 2.*hmax))[0][-1])+2
            zidxmax = np.min((zidxmax,len(self.frames['profiles']['record_yaml_ini'].air_ap.z.values)))
            zco = range(zidxmax)

            axes[label].plot(self.frames['profiles']['record_yaml_ini'].air_ap.q.values[zco], \
                             self.frames['profiles']['record_yaml_ini'].air_ap.z.values[zco],"b:", \
                             label="fit "+\
                             self.frames['profiles']['record_yaml_ini'].pars.ldatetime.strftime("%H:%M")\
                             +'LT')

            if self.path_obs is not None:
                zidxmax = int(np.where((self.frames['profiles']['record_yaml_obs_afternoon'].air_balloon.z.values < 2.*hmax))[0][-1])+2
                zidxmax = np.min((zidxmax,len(self.frames['profiles']['record_yaml_obs_afternoon'].air_balloon.z.values)))
                zco = range(zidxmax)


                axes[label].plot(self.frames['profiles']['record_yaml_obs_afternoon'].air_balloon.q.values[zco], \
                                 self.frames['profiles']['record_yaml_obs_afternoon'].air_balloon.z.values[zco],"r*", \
                                 label="obs "+\
                                 self.frames['profiles']['record_yaml_obs_afternoon'].pars.ldatetime.strftime("%H:%M")\
                                 +'LT')

                zidxmax = int(np.where((self.frames['profiles']['record_yaml_obs_afternoon'].air_ap.z.values < 2.*hmax))[0][-1])+2
                zidxmax = np.min((zidxmax,len(self.frames['profiles']['record_yaml_obs_afternoon'].air_ap.z.values)))
                zco = range(zidxmax)

                #print('r23')
                axes[label].plot(self.frames['profiles']['record_yaml_obs_afternoon'].air_ap.q.values[zco], \
                                 self.frames['profiles']['record_yaml_obs_afternoon'].air_ap.z.values[zco],"r:", \
                                 label="fit "+\
                                 self.frames['profiles']['record_yaml_obs_afternoon'].pars.ldatetime.strftime("%H:%M")\
                                 +'LT')

            #print('r24')
            if valid_mod:
                zidxmax = int(np.where((self.frames['profiles']['record_yaml_mod'].air_ap.z.values < 2.*hmax))[0][-1])+2
                zidxmax = np.min((zidxmax,len(self.frames['profiles']['record_yaml_mod'].air_ap.z.values)))
                zco = range(zidxmax)
                axes[label].plot(self.frames['profiles']['record_yaml_mod'].air_ap.q.values[zco], \
                                 self.frames['profiles']['record_yaml_mod'].air_ap.z.values[zco],"r-", \
                                 label="fit ")#+\
                             #self.frames['profiles']['record_yaml_mod'].pars.ldatetime.strftime("%H:%M")\
                             #+'LT')
            #print('r25')
            #axes[label].legend()

            #axes[label].legend(prop={'family':'monospace'},loc='upper left')
            #axes[label].set_ylabel('height [m]')
            axes[label].set_xlabel('q [kg/kg]')

            # #axes[label].set_title(self.evening_sounding.datetime.strftime("%Y/%m/%d %H:%M"))
            # zco =  self.evening_sounding.obs.z_pro < 2.*hmax
            # axes[label].plot(self.evening_sounding.obs.theta_pro[zco], self.evening_sounding.obs.z_pro[zco],"r*",label="obs "+self.evening_sounding.ldatetime.strftime("%H:%M")+'LT')
            # 
            # zco =  self.evening_sounding.fit.z_pro < 2.*hmax
            # axes[label].plot(self.evening_sounding.fit.theta_pro[zco], self.evening_sounding.fit.z_pro[zco],"r:",label="fit "+self.evening_sounding.ldatetime.strftime("%H:%M")+'LT')
            # 
            # zco = self.morning_sounding.c4gl.z_pro < 2.*hmax
            # axes[label].plot(self.morning_sounding.c4gl.theta_pro[zco], self.morning_sounding.c4gl.z_pro[zco],"r-",label="mod "+self.evening_sounding.ldatetime.strftime("%H:%M")+'LT')

            # #pl.subplots_adjust(right=0.6)

            # label = 'q_pro'
            # axes[label].clear()

            # hmax = np.max([self.morning_sounding.c4gl.input.h,self.morning_sounding.c4gl.out.h[-1],self.evening_sounding.fit.h])
            # 
            # zco =  self.morning_sounding.obs.z_pro < 2.*hmax
            # axes[label].plot(self.morning_sounding.obs.q_pro[zco], self.morning_sounding.obs.z_pro[zco],"b*",label="obs")
            # 
            # zco =  self.morning_sounding.c4gl.input.z_pro < 2.*hmax
            # axes[label].plot(self.morning_sounding.c4gl.input.q_pro[zco], self.morning_sounding.c4gl.input.z_pro[zco ],"b:",label="fit")

            # #self.ax5.set_title(self.evening_sounding.ldatetime.strftime("local time: %H:%M"))
            # zco =  self.evening_sounding.obs.z_pro < 2.*hmax
            # axes[label].plot(self.evening_sounding.obs.q_pro[zco], self.evening_sounding.obs.z_pro[zco],"r*",label="obs")
            # 
            # zco =  self.evening_sounding.fit.z_pro < 2.*hmax
            # axes[label].plot(self.evening_sounding.fit.q_pro[zco], self.evening_sounding.fit.z_pro[zco],"r:",label="fit")
            # 
            # zco = self.morning_sounding.c4gl.z_pro < 2.*hmax
            # axes[label].plot(self.morning_sounding.c4gl.q_pro[zco], self.morning_sounding.c4gl.z_pro[zco],"r-",label="mod")
            # #pl.subplots_adjust(right=0.6)
            # axes[label].set_xlabel('specific humidity [kg/kg]')
 

            #print('r26')
            time = self.frames['profiles']['record_yaml_mod'].out.time
            for ilabel,label in enumerate(['h','theta','q']):
                axes["out:"+label].clear()
                axes["out:"+label].plot(time,self.frames['profiles']['record_yaml_mod'].out.__dict__[label],label=label)
                axes["out:"+label].set_ylabel(label)
                if ilabel == 2:
                    axes["out:"+label].set_xlabel('local sun time [h]')
                
            #print('r27')
            label = 'SEB'
            axes[label].clear()
            
            axes[label].plot(time,self.frames['profiles']['record_yaml_mod'].out.Swin - self.frames['profiles']['record_yaml_mod'].out.Swout,label='Sw')
            axes[label].plot(time,-self.frames['profiles']['record_yaml_mod'].out.H,label='H')
            axes[label].plot(time,self.frames['profiles']['record_yaml_mod'].out.Lwin - self.frames['profiles']['record_yaml_mod'].out.Lwout,label='Lw')
            axes[label].plot(time,-self.frames['profiles']['record_yaml_mod'].out.G,label='G')
            axes[label].plot(time,-self.frames['profiles']['record_yaml_mod'].out.LE,label='LE')
            axes[label].hlines(0.,*axes[label].get_xlim(),'k')
            axes[label].set_ylabel('energy flux [$\mathrm{W/m^2}$]')
            axes[label].set_xlabel('local sun time [$\mathrm{h}$]')
                
            #print('r28')
            
            axes[label].legend()
            
            #         for ax in self.fig_timeseries_axes:
#             ax.clear()
#         
#         self.fig_timeseries_axes[0].plot(self.morning_sounding.c4gl.out.h,label='h')
#         self.fig_timeseries_axes[1].plot(self.morning_sounding.c4gl.out.theta,label='theta')
#         self.fig_timeseries_axes[2].plot(self.morning_sounding.c4gl.out.q,label='q')
#         #print(self.morning_sounding.c4gl.out.Swin)
#         self.fig_timeseries_axes[3].plot(self.morning_sounding.c4gl.out.Swin - self.morning_sounding.c4gl.out.Swout,label='Sw')
#         self.fig_timeseries_axes[3].plot(-self.morning_sounding.c4gl.out.H,label='H')
#         self.fig_timeseries_axes[3].plot(self.morning_sounding.c4gl.out.Lwin - self.morning_sounding.c4gl.out.Lwout,label='Lw')
#         self.fig_timeseries_axes[3].plot(-self.morning_sounding.c4gl.out.G,label='G')
#         self.fig_timeseries_axes[3].plot(-self.morning_sounding.c4gl.out.LE,label='LE')
#         self.fig_timeseries_axes[3].hlines(0.,*self.fig_timeseries_axes[3].get_xlim(),'k')
#         self.fig_timeseries_axes[3].legend()
#         self.fig.canvas.draw()
            






        #self.ready()
        #print('r29')
        fig.canvas.draw()
        #fig.show()

        self.axes = axes
        self.tbox = tbox
        self.fig = fig

    def on_pick(self,event):
        #print("HELLO")
        # this makes clear that the dataset is loading (set_profile_focus takes a long time to load!)
        #self.axes['theta_pro'].clear()
        #self.axes['q_pro'].clear()
        

        # workaround because I cannot track the axes label here. I need it because the behaviour of this function should depend on which axes we are.
        # I can only track the label of the data points. So we make a definition that clarifies to which axes the select data points (having a 'key') belongs to
        keys_to_axes = {}
        for ikey,key in enumerate(self.frames['stats']['viewkeys']):
            keys_to_axes['d'+self.frames['stats']['viewkeys'][ikey]+'dt'] = 'stats_d'+key+'dt'

        keys_to_axes['worldmap_stations'] = 'worldmap_stations'
        keys_to_axes['worldmap'] = 'worldmap'
        
        axes = self.axes
        #nstatsview = self.nstatsview
        #statsviewcmap = self.statsviewcmap
        stations = self.frames['worldmap']['stations'].table


        #print("p1")
        current = event
        artist = event.artist
        
        selkey = artist.get_label()
        
        #print(keys_to_axes)
        
        label = keys_to_axes[selkey]
        #print("HELLO",selkey,label)

        # # Get to know in which axes we are
        # label = None
        # for axeskey in axes.keys():
        #     if event.inaxes == axes[axeskey]:
        #         label = axeskey
        #         

        # cont, pos = None, None
        
        xmouse, ymouse = event.mouseevent.xdata, event.mouseevent.ydata
        ind = event.ind
        # x, y = artist.get_xdata(), artist.get_ydata() # for some reason this doesnt work yet :/
        d = axes[label].collections[0]
        #d.set_offset_position('data')
        xy = d.get_offsets()
        x, y =  xy[:,0],xy[:,1]
        #axes[-1].plot(seltableoutput[key+'_obs']*3600.,seltableoutput[key+'_mod']*3600.,'ro', markersize=5, picker=5,label=key)

        #print("p2")
        if len(ind) > 0:
            #print("p3")
            pos = x[ind[0]], y[ind[0]]

            #if label[:-1] == 'statsview':
            #    #seltablestatsstdrel = self.seltablestatsstdrel
            #    #seltablestatspct = self.seltablestatspct

            #    #self.set_statsviewfocus('STNID' seltablestatsstdrel[(seltablestatsstdrel[selkey+'_ext'] == pos[0]) & (seltablestatsstdrel[selkey+'_mod'] == pos[1] )  ].STNID.iloc[0]
            #    #self.set_statsviewfocus('DT'] = seltablestatsstdrel[(seltablestatsstdrel[selkey+'_ext'] == pos[0]) & (seltablestatsstdrel[selkey+'_mod'] == pos[1] )  ].DT.iloc[0]
            #    
            #    self.axes['worldmap'].focus['STNID'] = self.axes['statsview0'].focus['STNID']
            #    self.set_profilefocus(STNID=self.axes['statsview0'].focus['STNID'],DT=self.axes['statsview0'].focus['DT'])
            #    self.goto_datetime_worldmap(self.profilefocus['DT'],'after')
            #    
            #    self.refresh_plot_interface(only=['statsviews_lightupdate','worldmap','profiles'],statsnewdata=False)
            #el
            if (label == 'worldmap') or (label == 'worldmap_stations'):
                self.hover_active = False
                if (self.frames['worldmap']['STNID'] !=
                    self.frames['profiles']['STNID']):
                # WE ALREADY HAVE the correct station from worldmap/stats because of the hovering!!
                # so we just need to perform update_station
                    self.update_station()
            elif (label[:5] == 'stats'):

                self.hover_active = False
                if (self.frames['stats']['STNID'] !=
                self.frames['profiles']['STNID']) or \
                   (self.frames['stats']['current_record_chunk'] != 
                    self.frames['profiles']['current_record_chunk']) or \
                   (self.frames['stats']['current_record_index'] != 
                    self.frames['profiles']['current_record_index']):



                    for key in ['STNID','current_station','stations_iterator']: 
                        self.frames['worldmap'][key] = self.frames['stats'][key] 

                    for key in self.frames['stats'].keys():
                        self.frames['profiles'][key] = self.frames['stats'][key]

                    STNID = self.frames['profiles']['STNID']
                    chunk = self.frames['profiles']['current_record_chunk']
                    if 'current_station_file_ini' in self.frames['profiles'].keys():
                        self.frames['profiles']['current_station_file_ini'].close()
                    self.frames['profiles']['current_station_file_ini'] = \
                        open(self.path_exp+'/'+format(STNID,"05d")+'_'+str(chunk)+'_ini.yaml','r')

                    if 'current_station_file_mod' in self.frames['profiles'].keys():
                        self.frames['profiles']['current_station_file_mod'].close()
                    self.frames['profiles']['current_station_file_mod'] = \
                        open(self.path_exp+'/'+format(STNID,"05d")+'_'+str(chunk)+'_mod.yaml','r')
                    if self.path_obs is not None:
                        if 'current_station_file_afternoon' in self.frames['profiles'].keys():
                            self.frames['profiles']['current_station_file_afternoon'].close()
                        self.frames['profiles']['current_station_file_afternoon'] = \
                            open(self.path_obs+'/'+format(STNID,"05d")+'_afternoon.yaml','r')

                    # go to hovered record of current station
                    self.frames['profiles']['records_iterator'] = \
                                    records_iterator(self.frames['profiles']['records_current_station_mod'])
                    # ... and go to the record of the profile window (last one that
                    # was picked by the user)
                    found = False
                    EOF = False
                    while (not found) and (not EOF):
                        try:
                            (STNID,chunk,index),record = self.frames['profiles']['records_iterator'].__next__()
                            #print("hello*")
                            #print(self.frames['profiles']['current_record_index'])
                            if (chunk == self.frames['profiles']['current_record_chunk']) and \
                               (index == self.frames['profiles']['current_record_index']) and \
                               (STNID == self.frames['profiles']['STNID']):
                                #print('found!')
                                found = True
                        except StopIteration:
                            EOF = True
                    if found:
                        self.frames['stats']['current_record_mod'] = record
                        self.frames['stats']['current_record_chunk'] = chunk
                        self.frames['stats']['current_record_index'] = index
                    # # for the profiles we make a distinct record iterator, so that the
                    # # stats iterator can move independently
                    # self.frames['profiles']['records_iterator'] = \
                    #                 records_iterator(self.frames['profiles']['records_current_station_mod'])
                    # (self.frames['profiles']['STNID'] , \
                    # self.frames['profiles']['current_record_index']) , \
                    # self.frames['profiles']['current_record_mod'] = \
                    #                 self.frames['profiles']['records_iterator'].__next__()


                    # for the profiles we make a distinct record iterator, so that the
                    # stats iterator can move independently

                    self.update_record()



    def on_plot_hover(self,event):
        axes = self.axes
        #print('h1')

        # Get to know in which axes we are
        label = None
        for axeskey in axes.keys():
            if event.inaxes == axes[axeskey]:
                label = axeskey
                
        #print('h2')

        cont, pos = None, None
        #print (label)
        
        if label is not None:
            if  ('data' in axes[label].__dict__.keys()) and \
                (label in axes[label].data.keys()) and \
                (axes[label].data[label] is not None):
                
                #print('h3')
                cont, ind =  axes[label].data[label].contains(event)
                selkey = axes[label].data[label].get_label()
                if len(ind["ind"]) > 0:
                    #print('h4')
                    pos = axes[label].data[label].get_offsets()[ind["ind"][0]]
                    #print('pos',pos,selkey)


                    #if label[:-1] == 'statsview':
                    #    seltablestatsstdrel = self.seltablestatsstdrel
                    #    seltablestatspct = self.seltablestatspct

                    #    self.set_statsviewfocus('STNID'] = seltablestatsstdrel[(seltablestatsstdrel[selkey+'_ext'] == pos[0]) & (seltablestatsstdrel[selkey+'_mod'] == pos[1] )  ].STNID.iloc[0]
                    #    self.set_statsviewfocus('DT'] = seltablestatsstdrel[(seltablestatsstdrel[selkey+'_ext'] == pos[0]) & (seltablestatsstdrel[selkey+'_mod'] == pos[1] )  ].DT.iloc[0]
                    #    self.axes['worldmap'].focus['STNID'] = self.axes['statsview0'].focus['STNID']
                    #    #self.goto_datetime_worldmap(self.axes['statsview0'].focus['DT'],'after')
                    #    self.hover_active = True
                    #    
                    #    self.refresh_plot_interface(only=['statsviews_lightupdate','worldmap_stations'])
                    #    
                    #el
                    #print(label[:5])
                    if (label[:5] == 'stats') or (label == 'times'):
                        # records_mod = self.frames['stats']['records_current_station_mod'][selkey]
                        # records_obs = self.frames['stats']['records_current_station_obs_afternoon'][selkey]
                        
                        if self.path_obs is not None:
                            if label[:5] == 'stats':
                                records_mod_stats = self.frames['stats']['records_all_stations_mod_stats']
                                records_obs_stats = self.frames['stats']['records_all_stations_obs_afternoon_stats']
                                (self.frames['stats']['STNID'] ,
                                 self.frames['stats']['current_record_chunk'], 
                                 self.frames['stats']['current_record_index']) = \
                                    records_mod_stats[(records_obs_stats[selkey] == pos[0]) & (records_mod_stats[selkey] == pos[1])].index[0]
                            # elif label[:5] == 'stats':
                            #     # records_mod_stats = self.frames['stats']['records_all_stations_mod_stats']
                            #     records_obs_stats = self.frames['stats']['records_all_stations_obs_afternoon_stats']
                            #     records_datetimes = self.frames['stats']['records_all_stations_ini']
                            #     (self.frames['stats']['STNID'] ,
                            #      self.frames['stats']['current_record_chunk'], 
                            #      self.frames['stats']['current_record_index']) = \
                            #         records_mod_stats[(records_obs_stats[selkey] == pos[0]) & (records_mod_stats[selkey] == pos[1])].index[0]


                        self.frames['stats']['stations_iterator'] = stations_iterator(self.frames['worldmap']['stations']) 
                        
                        # # TO TEST: should be removed, since it's is also done just below
                        # self.frames['stats']['stations_iterator'] = \
                        #     self.frames['worldmap']['stations_iterator'] 
                
                
                        # self.goto_datetime_worldmap(
                        #     self.frames['profiles']['current_record_obs'].datetime.to_pydatetime(),
                        #     'after')


                        # scrolling to the right station
                        STNID,station = self.frames['stats']['stations_iterator'].__next__()
                        EOF = False
                        found = False
                        while (not found and not EOF):
                            if (STNID == self.frames['stats']['STNID']):
                                   found = True 
                            if not found:
                                try:
                                    STNID,station = self.frames['stats']['stations_iterator'].__next__()
                                except (StopIteration):
                                    EOF = True
                        if found:
                        #    self.frames['stats']['STNID'] = STNID
                            self.frames['stats']['current_station'] =  station

                        #STNID = self.frames['profiles']['current_record_index'].iloc[0].name[0]
                        #index = self.frames['profiles']['current_record_index'].iloc[0].name[1]


                        # generate index of the current station
                        self.frames['stats']['records_current_station_index'] = \
                            (self.frames['stats']['records_all_stations_index'].get_level_values('STNID')\
                             == self.frames['stats']['STNID'])


                        tab_suffixes = \
                                ['_mod','_ini','_ini_pct']
                        if self.path_obs is not None:
                            tab_suffixes += \
                                ['_mod_stats','_obs_afternoon','_obs_afternoon_stats']
                            
                        for tab_suffix in tab_suffixes:
                            self.frames['stats']['records_current_station'+tab_suffix] = \
                                self.frames['stats']['records_all_stations'+tab_suffix].iloc[self.frames['stats']['records_current_station_index']]



                        # go to hovered record of current station
                        self.frames['stats']['records_iterator'] = \
                                        records_iterator(self.frames['stats']['records_current_station_mod'])


                        # ... and go to the record of the profile window (last one that
                        # was picked by the user)
                        found = False
                        EOF = False
                        while (not found) and (not EOF):
                            try:
                                (STNID,chunk,index),record = self.frames['stats']['records_iterator'].__next__()
                                #print("hello*")
                                #print(self.frames['profiles']['current_record_index'])
                                if (index == self.frames['stats']['current_record_index']) and \
                                   (chunk == self.frames['stats']['current_record_chunk']) and \
                                   (STNID == self.frames['stats']['STNID']):
                                    #print('found!')
                                    found = True
                            except StopIteration:
                                EOF = True
                        if found:
                            #print('h5')
                            self.frames['stats']['current_record_mod'] = record
                            self.frames['stats']['current_record_chunk'] = chunk
                            self.frames['stats']['current_record_index'] = index

                        #print(self.frames['stats']['STNID'],self.frames['stats']['current_record_index'])
                        tab_suffixes = \
                                ['_ini','_ini_pct']
                        if self.path_obs is not None:
                            tab_suffixes += \
                                ['_mod_stats','_obs_afternoon','_obs_afternoon_stats']
                        for tab_suffix in tab_suffixes:
                            #print(tab_suffix)
                            #print(self.frames['stats']['records_current_station'+tab_suffix])
                            self.frames['stats']['current_record'+tab_suffix] =  \
                                self.frames['stats']['records_current_station'+tab_suffix].loc[\
                                      (self.frames['stats']['STNID'] , \
                                       self.frames['stats']['current_record_chunk'] , \
                                       self.frames['stats']['current_record_index'])]


                        self.hover_active = True
                        self.refresh_plot_interface(only=['stats_lightupdate','worldmap_stations','profiles'])
                        # print('h13')
                        # if 'time' in self.globaldata.datasets[key].page[key].dims:
                        #     self.goto_datetime_worldmap(
                        #         self.frames['profiles']['current_record_ini'].datetime.to_pydatetime(),
                        #         'after')
                        #     if "fig" in self.__dict__.keys():
                        #         self.refresh_plot_interface(only=['stats_lightupdate',
                        #                                           'worldmap',
                        #                                           'profiles'])
                        # else:
                        #     if "fig" in self.__dict__.keys():
                        #         self.refresh_plot_interface(only=['stats_lightupdate',
                        #                                           'worldmap_stations',
                        #                                           'profiles'])



                    elif label in ['worldmap_stations','worldmap']:
                        #print('h5')

                        if (self.axes['worldmap'].lat is not None) and \
                           (self.axes['worldmap'].lon is not None):


                            #self.loading()
                            self.fig.canvas.draw()
                            self.fig.show()


                            # get position of 
                            latmap = round(pos[1]/len(self.axes['worldmap'].lat)*(self.axes['worldmap'].lat[-1] - \
                                                                 self.axes['worldmap'].lat[0]) + \
                                           self.axes['worldmap'].lat[0],4)
                            lonmap = round(pos[0]/len(self.axes['worldmap'].lon)*(self.axes['worldmap'].lon[-1] - \
                                                                 self.axes['worldmap'].lon[0]) + \
                                           self.axes['worldmap'].lon[0],4)
                        
                            stations = self.frames['worldmap']['stations'].table
                            #print('h7')
                        
                            #reset stations iterator:
                            # if 'stations_iterator' in self.frames['worldmap'].keys():
                            #     self.frames['worldmap']['stations_iterator'].close()
                            #     del(self.frames['worldmap']['stations_iterator'])
                            # if 'stations_iterator' in self.frames['stats'].keys():
                            #     self.frames['stats']['stations_iterator'].close()
                            #     del(self.frames['stats']['stations_iterator'])
                            self.frames['worldmap']['stations_iterator'] =\
                               stations_iterator(self.frames['worldmap']['stations'])
                            STNID,station = self.frames['worldmap']['stations_iterator'].__next__()
                            EOF = False
                            found = False
                            while (not found and not EOF):
                                #print('h8',station.latitude,latmap)
                                #print('h8',station.longitude,lonmap)
                                if (round(station.latitude,3) == round(latmap,3)) and \
                                    (round(station.longitude,3) == round(lonmap,3)):
                                       found = True 
                                if not found:
                                    try:
                                        STNID,station = self.frames['worldmap']['stations_iterator'].__next__()
                                    except (StopIteration):
                                        EOF = True
                            if found:
                                self.frames['worldmap']['STNID'] = STNID
                                self.frames['worldmap']['current_station'] = \
                                        station
                        
                            self.frames['stats']['stations_iterator'] = \
                                self.frames['worldmap']['stations_iterator'] 
                            #print('h8')
                            # inherit station position for the stats frame...
                            for key in self.frames['worldmap'].keys():
                                self.frames['stats'][key] = self.frames['worldmap'][key]
                                
                            ## fetch records of current station...
                            #self.frames['stats']['records_current_station_mod'] =\
                            #   get_records_mod(pd.DataFrame([self.frames['stats']['current_station']]),self.path_exp)

                            # ... and their indices
                            self.frames['stats']['records_current_station_index'] = \
                                    (self.frames['stats']['records_all_stations_index'].get_level_values('STNID')\
                                     == \
                                     self.frames['stats']['current_station'].name)

                            tab_suffixes = \
                                    ['_mod','_ini','_ini_pct']
                            if self.path_obs is not None:
                                tab_suffixes += \
                                    ['_mod_stats','_obs_afternoon','_obs_afternoon_stats']

                            for tab_suffix in tab_suffixes:
                                self.frames['stats']['records_current_station'+tab_suffix] = \
                                    self.frames['stats']['records_all_stations'+tab_suffix].iloc[self.frames['stats']['records_current_station_index']]


                            # ... create a record iterator ...
                            #self.frames['stats']['records_iterator'].close()
                            del(self.frames['stats']['records_iterator'])
                            self.frames['stats']['records_iterator'] = \
                                self.frames['stats']['records_current_station_mod'].iterrows()



                        
                            #print('h9')
                            # ... and go to to the first record of the current station
                            (self.frames['stats']['STNID'] , \
                             self.frames['stats']['current_record_chunk'] , \
                             self.frames['stats']['current_record_index']) , \
                            self.frames['stats']['current_record_mod'] = \
                                self.frames['stats']['records_iterator'].__next__()
                        
                            tab_suffixes = \
                                    ['_ini','_ini_pct']
                            if self.path_obs is not None:
                                tab_suffixes += \
                                    ['_mod_stats','_obs_afternoon','_obs_afternoon_stats']

                            for tab_suffix in tab_suffixes:
                                self.frames['stats']['current_record'+tab_suffix] =  \
                                    self.frames['stats']['records_current_station'+tab_suffix].loc[\
                                          (self.frames['stats']['STNID'] , \
                                           self.frames['stats']['current_record_chunk'] , \
                                           self.frames['stats']['current_record_index'])]

                            #print('h11')
                            
                            self.hover_active = True
                            self.refresh_plot_interface(only=['stats_lightupdate','worldmap_stations','profiles'])
                            #print('h13')

                        

            #if (stations is not None):
            #    for iSTN,STN in stations.iterrows():
            ##        x,y =self.gmap(STN['longitude'],STN['latitude'])
            ##        self.gmap.plot(x,y, 'mo' if (self.STNID == STN['ID']) else 'ro' , markersize=1)
            #        x,y = len(axes[label].lon)*(STN['longitude']- axes[label].lon[0])/(axes[label].lon[-1] - axes[label].lon[0])  ,len(lat)*(STN['latitude']- axes[label].lat[0])/(lat[-1] - axes[label].lat[0])
            #        axes['worldmap'].plot(x,y, 'mo' if (axes['worldmap'].focus['STNID'] == STN['ID']) else 'ro' , markersize=2)

        # self.fig.show()
 
        # we are hovering on nothing, so we are going back to the position of
        # the profile sounding
        if pos is None:
            if self.hover_active == True:
                #print('h1*')
                
                #self.loading()
                # to do: reset stations iterators

                # get station and record index from the current profile
                for key in ['STNID', 'current_station']:
                    self.frames['stats'][key] = self.frames['profiles'][key]

                self.frames['stats']['STNID'] = self.frames['profiles']['STNID']
                self.frames['stats']['current_station'] = \
                        self.frames['profiles']['current_station']
                #print('h3a*')
                self.frames['stats']['records_current_station_mod'] = \
                        self.frames['profiles']['records_current_station_mod']
                #print('h3b*')

                # the next lines recreate the records iterator. Probably it's
                # better to just copy the profile iterator and its position to
                # the worldmap/stats 

                # reset stations iterator...
                #self.frames['stats']['records_iterator'].close()
                del(self.frames['stats']['records_iterator'])
                self.frames['stats']['records_iterator'] = \
                    self.frames['stats']['records_current_station_mod'].iterrows()
                #print('h4*')

                # ... and go to the record of the profile window (last one that
                # was picked by the user)
                found = False
                EOF = False
                while (not found) and (not EOF):
                    try:
                        (STNID,chunk,index),record = self.frames['stats']['records_iterator'].__next__()
                        #print("hello*")
                        #print(self.frames['profiles']['current_record_index'])
                        #print(self.frames['profiles']['STNID'])
                        #print(STNID,index)
                        if (index == self.frames['profiles']['current_record_index']) and \
                            (chunk == self.frames['profiles']['current_record_chunk']) and \
                            (STNID == self.frames['profiles']['STNID']):
                            #print('found!')
                            found = True
                    except StopIteration:
                        EOF = True
                if found:
                    #print('h5*')
                    self.frames['stats']['current_record_mod'] = record
                    self.frames['stats']['current_record_chunk'] = chunk
                    self.frames['stats']['current_record_index'] = index

                #print('h6*')



                # # fetch records of current station...
                # self.frames['stats']['records_current_station_mod'] =\
                #    get_records_mod(pd.DataFrame([self.frames['stats']['current_station']]),self.path_exp)

                # ... and their indices
                self.frames['stats']['records_current_station_index'] = \
                        (self.frames['stats']['records_all_stations_index'].get_level_values('STNID')\
                         == \
                         self.frames['stats']['current_station'].name)

                
                tab_suffixes = \
                        ['_ini','_ini_pct']
                if self.path_obs is not None:
                    tab_suffixes += \
                        ['_mod_stats','_obs_afternoon','_obs_afternoon_stats']

                for tab_suffix in tab_suffixes:
                    self.frames['stats']['records_current_station'+tab_suffix] = \
                        self.frames['stats']['records_all_stations'+tab_suffix].iloc[self.frames['stats']['records_current_station_index']]


                for tab_suffix in tab_suffixes:
                    self.frames['stats']['current_record'+tab_suffix] =  \
                        self.frames['stats']['records_current_station'+tab_suffix].loc[\
                              (self.frames['stats']['STNID'] , \
                               self.frames['stats']['current_record_chunk'] , \
                               self.frames['stats']['current_record_index'])]


                # the next lines recreate the stations iterator. Probably it's
                # better to just copy the profile iterator and its position to
                # the worldmap/stats 
                #print('h7*')

                # reset the stations iterators
                for framekey in ['stats','worldmap']:
                    ##print(framekey)
                    if 'stations_iterator' in self.frames[framekey]:
                        #self.frames[framekey]['stations_iterator'].close()
                        del(self.frames[framekey]['stations_iterator'])

                self.frames['worldmap']['current_station'] = \
                        self.frames['profiles']['current_station']

                #recreate the stations iterator for the worldmap...
                self.frames['worldmap']['stations_iterator'] = stations_iterator(self.frames['worldmap']['stations']) 

                # ... and go the position of the profile
                #print('h8*')
                STNID,station = self.frames['worldmap']['stations_iterator'].__next__()
                EOF = False
                found = False
                while (not found and not EOF):
                    if STNID == self.frames['profiles']['STNID'] :
                        found = True 
                    if not found:
                        try:
                            STNID,station = self.frames['worldmap']['stations_iterator'].__next__()
                        except (StopIteration):
                            EOF = True
                if found:
                    self.frames['worldmap']['current_station'] = station
                    self.frames['worldmap']['STNID'] = STNID
                #print('h9*')
                self.frames['stats']['stations_iterator'] = \
                    self.frames['worldmap']['stations_iterator'] 

                # the stats window now inherits the current station from the
                # worldmap
                for key in ['STNID','current_station','stations_iterator']: 
                    self.frames['stats'][key] = self.frames['worldmap'][key] 
                #print('h10*')

                # # we now only need inherit station position and go to first record
                # for key in self.frames['worldmap'].keys():
                #     self.frames['stats'][key] = self.frames['worldmap'][key]

                # self.frames['stats']['records_current_station'] =\
                #     get_records(pd.DataFrame().append(self.frames['stats']['current_station']))

                # #print(self.frames['stats']['records_current_station'])
                # self.frames['stats']['records_iterator'] = \
                #                 self.frames['stats']['records_current_station'].iterrows()
                # (self.frames['stats']['STNID'] , \
                # self.frames['stats']['current_record_index']) , \
                # self.frames['stats']['current_record_mod'] = \
                #                 self.frames['stats']['records_iterator'].__next__()
                






                #self.set_statsviewfocus('STNID', self.profilefocus['STNID'])
                ##self.set_statsviewfocus('DT'], self.profilefocus['DT'])
                #self.axes['worldmap'].focus['STNID'] = self.profilefocus['STNID']
                ##self.goto_datetime_worldmap(self.profilefocus['DT'],'after')
                self.hover_active = False
                self.refresh_plot_interface(only=['stats_lightupdate','worldmap_stations'],statsnewdata=False)
    # def loading(self):
    #     self.tbox['loading'].set_text('Loading...')
    #     self.fig.canvas.draw()
    #     self.fig.show()
    #     sleep(0.1)
    # def ready(self):
    #     self.tbox['loading'].set_text('Ready')
    #     self.fig.canvas.draw()
    #     self.fig.show()



