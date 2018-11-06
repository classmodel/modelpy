#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thursday, March 29, 11:30 AM

@author: Hendrik Wouters

The dry-2-dry global radio sounding experiment.

usage:
    python setup_global.py <number>
    where <number> is an integer indicating the row index of the station list
    under args.path_output+'/'+fn_stations (see below)

this scripts should be called from the pbs script setup_global.pbs



dependencies:
    - pandas
    - class4gl
    - data_soundings


"""

""" import libraries """
import pandas as pd
import sys
#import copy as cp
import numpy as np
#from sklearn.metrics import mean_squared_error
import logging
import datetime as dt
import os
import math

import argparse

#if __name__ == '__main__':
parser = argparse.ArgumentParser()
#parser.add_argument('--timestamp')
parser.add_argument('--path_input')#,default='/user/data/gent/gvo000/gvo00090/D2D/data/SOUNDINGS/')
parser.add_argument('--path_output')#,default='/user/data/gent/gvo000/gvo00090/D2D/data/C4GL/')
# parser.add_argument('--first_YYYYMMDD',default="19810101")
# parser.add_argument('--last_YYYYMMDD',default="20180101")
parser.add_argument('--startyear',default="1981")
parser.add_argument('--first_station_row')
parser.add_argument('--last_station_row')
parser.add_argument('--c4gl_path_lib')#,default='/user/data/gent/gvo000/gvo00090/D2D/software/CLASS/class4gl/lib')
parser.add_argument('--station_id') # run a specific station id
# parser.add_argument('--error_handling',default='dump_on_success')
# parser.add_argument('--subset_output',default='morning') # this tells which yaml subset


# args.path_output = "/user/data/gent/gvo000/gvo00090/D2D/data/SOUNDINGS/GLOBAL/"
args = parser.parse_args()

sys.path.insert(0, args.c4gl_path_lib)
from interface_multi import stations,stations_iterator, records_iterator,get_record_yaml,get_records

fn_stations = args.path_input+'/igra-stations.txt'


#calculate the root mean square error

def rmse(y_actual,y_predicted,z_actual = None, z_predicted = None,filternan_actual = False):
    """ calculated root mean squared error 
        
    
        INPUT:
            y_actual: reference dataset
            y_predicted: predicting dataset
            z_actual: coordinate values of reference dataset
            z_predicted: coordinate values of the predicting dataset
            
            filternan_actual: throw away reference values that have nans
    """
    
    y_actual_temp = np.array(y_actual)
    y_predicted_temp = np.array(y_predicted)
    
    if z_actual is not None:
        z_actual_temp = np.array(z_actual)
    else: 
        z_actual_temp = None
        
    
    if filternan_actual:
        y_actual_temp = y_actual_temp[~np.isnan(y_actual_temp)]
        if z_actual_temp is not None:
            z_actual_temp = z_actual_temp[~np.isnan(y_actual_temp)]
    
    if ((z_actual_temp is not None) or (z_predicted is not None)):    
        if (z_actual_temp is None) or (z_predicted is None):
            raise ValueError('Input z_actual and z_predicted need \
                              to be specified simultaneously.')
        y_predicted_temp = np.interp(z_actual_temp,z_predicted, y_predicted)
    
    else:
        # this catches the situation that y_predicted is a single value (eg., 
        # which is the case for evaluating eg., mixed-layer estimates)
        y_predicted_temp = y_actual_temp*0. + y_predicted_temp
        
    rmse_temp = (y_actual_temp - y_predicted_temp)
    rmse_temp = np.mean(rmse_temp*rmse_temp)
    return np.sqrt(rmse_temp)


from class4gl import class4gl_input, data_global,class4gl
from data_soundings import wyoming
#from data_global import data_global

# iniitialize global data
globaldata = data_global()
# ...  and load initial data pages
globaldata.load_datasets(recalc=0)

# read the list of stations with valid ground data (list generated with
# get_valid_stations.py)
# args.path_input = "/user/data/gent/gvo000/gvo00090/EXT/data/SOUNDINGS/"

# df_stations = pd.read_fwf(fn_stations,names=['Country code',\
#                                                'ID',\
#                                                'Name',\
#                                                'latitude',\
#                                                'longitude',\
#                                                'height',\
#                                                'unknown',\
#                                                'startyear',\
#                                                'endyear'])
# 

# ===============================
print("getting a list of stations")
# ===============================
all_stations = stations(args.path_input,refetch_stations=False)
df_stations = all_stations.table
df_stations.columns

if args.station_id is not None:
    df_stations = df_stations.query('STNID == '+args.station_id)
else:
    if args.last_station_row is not None:
        df_stations = df_stations[:(int(args.last_station_row)+1)]
    if args.first_station_row is not None:
        df_stations = df_stations[int(args.first_station_row):]

STNlist = list(df_stations.iterrows())

os.system('mkdir -p '+args.path_output)
for iSTN,STN in STNlist:  
    one_run = False
# for iSTN,STN in STNlist[5:]:  
    
    fnout = args.path_output+"/"+format(STN.name,'05d')+"_morning.yaml"
    fnout_afternoon = args.path_output+"/"+format(STN.name,'05d')+"_afternoon.yaml"
    

    # c4glfiles = dict([(EXP,odirexperiments[EXP]+'/'+format(STN['ID'],'05d')+'.yaml') \
    #                   for EXP in experiments.keys()])
        
    with open(fnout,'w') as fileout, \
         open(fnout_afternoon,'w') as fileout_afternoon:
        wy_strm = wyoming(PATH=args.path_input, STNM=STN.name)
        wy_strm.set_STNM(int(STN.name))

        # we consider all soundings from 1981 onwards
        wy_strm.find_first(year=int(args.startyear))
        #wy_strm.find(dt.datetime(2004,10,19,6))
        
        c4gli = class4gl_input(debug_level=logging.INFO)
        c4gli_afternoon = class4gl_input(debug_level=logging.INFO)
        # so we continue as long as we can find a new sounding
                
        while wy_strm.current is not None:
            
            c4gli.clear()
            try: 
                c4gli.get_profile_wyoming(wy_strm)
                #print(STN['ID'],c4gli.pars.datetime)
                #c4gli.get_global_input(globaldata)

                print(c4gli.pars.STNID, c4gli.pars.ldatetime)

                logic = dict()
                logic['morning'] =  (c4gli.pars.ldatetime.hour <= 12.)

                # Sounding should have taken place after 3 hours before sunrise.
                # Note that the actual simulation only start at sunrise
                # (specified by ldatetime_daylight), so the ABL cooling af the time
                # before sunrise is ignored by the simulation.
                logic['daylight'] = \
                    ((c4gli.pars.ldatetime - 
                      c4gli.pars.lSunrise).total_seconds()/3600. >= -3.)
                
                logic['springsummer'] = (c4gli.pars.theta > 278.)
                
                # we take 3000 because previous analysis (ie., HUMPPA) has
                # focussed towards such altitude
                le3000 = (c4gli.air_balloon.z <= 3000.)
                logic['10measurements'] = (np.sum(le3000) >= 10) 

                leh = (c4gli.air_balloon.z <= c4gli.pars.h)

                logic['mlerrlow'] = (\
                        (len(np.where(leh)[0]) > 0) and \
                        # in cases where humidity is not defined, the mixed-layer
                        # values get corr
                        (not np.isnan(c4gli.pars.theta)) and \
                        (rmse(c4gli.air_balloon.theta[leh] , \
                              c4gli.pars.theta,filternan_actual=True) < 1.)\
                              )
    

                logic['mlherrlow'] = (c4gli.pars.h_e <= 150.)
                
                print('logic:', logic)
                # the result
                morning_ok = np.mean(list(logic.values()))
                print(morning_ok,c4gli.pars.ldatetime)

            except:
                morning_ok =False
                print('obtain morning not good')

            # the next sounding will be used either for an afternoon sounding
            # or for the morning sounding of the next day.
            wy_strm.find_next()
            # If the morning is ok, then we try to find a decent afternoon
            # sounding
            if morning_ok == 1.:
                print('MORNING OK!')
                # we get the current date
                current_date = dt.date(c4gli.pars.ldatetime.year, \
                                       c4gli.pars.ldatetime.month, \
                                       c4gli.pars.ldatetime.day)
                c4gli_afternoon.clear()
                print('AFTERNOON PROFILE CLEARED')
                try:
                    c4gli_afternoon.get_profile_wyoming(wy_strm)
                    print('AFTERNOON PROFILE OK')

                    if wy_strm.current is not None:
                        current_date_afternoon = \
                                   dt.date(c4gli_afternoon.pars.ldatetime.year, \
                                           c4gli_afternoon.pars.ldatetime.month, \
                                           c4gli_afternoon.pars.ldatetime.day)
                    else:
                        # a dummy date: this will be ignored anyway
                        current_date_afternoon = dt.date(1900,1,1)

                    # we will dump the latest afternoon sounding that fits the
                    # minimum criteria specified by logic_afternoon
                    print(current_date,current_date_afternoon)
                    c4gli_afternoon_for_dump = None
                    while ((current_date_afternoon == current_date) and \
                           (wy_strm.current is not None)):
                        logic_afternoon =dict()

                        logic_afternoon['afternoon'] = \
                            (c4gli_afternoon.pars.ldatetime.hour >= 12.)
                        # the sounding should have taken place before 1 hours
                        # before sunset. This is to minimize the chance that a
                        # stable boundary layer (yielding very low mixed layer
                        # heights) is formed which can not be represented by
                        # class.
                        logic_afternoon['daylight'] = \
                          ((c4gli_afternoon.pars.ldatetime - \
                            c4gli_afternoon.pars.lSunset \
                           ).total_seconds()/3600. <= 1.)


                        le3000_afternoon = \
                            (c4gli_afternoon.air_balloon.z <= 3000.)
                        logic_afternoon['5measurements'] = \
                            (np.sum(le3000_afternoon) >= 5) 

                        # we only store the last afternoon sounding that fits these
                        # minimum criteria

                        afternoon_ok = np.mean(list(logic_afternoon.values()))

                        print('logic_afternoon: ',logic_afternoon)
                        print(afternoon_ok,c4gli_afternoon.pars.ldatetime)
                        if afternoon_ok == 1.:
                            # # doesn't work :(
                            # c4gli_afternoon_for_dump = cp.deepcopy(c4gli_afternoon)
                            
                            # so we just create a new one from the same wyoming profile
                            c4gli_afternoon_for_dump = class4gl_input()
                            c4gli_afternoon_for_dump.get_profile_wyoming(wy_strm)

                        wy_strm.find_next()
                        c4gli_afternoon.clear()
                        c4gli_afternoon.get_profile_wyoming(wy_strm)

                        if wy_strm.current is not None:
                            current_date_afternoon = \
                                   dt.date(c4gli_afternoon.pars.ldatetime.year, \
                                           c4gli_afternoon.pars.ldatetime.month, \
                                           c4gli_afternoon.pars.ldatetime.day)
                        else:
                            # a dummy date: this will be ignored anyway
                            current_date_afternoon = dt.date(1900,1,1)

                        # Only in the case we have a good pair of soundings, we
                        # dump them to disk
                    if c4gli_afternoon_for_dump is not None:
                        c4gli.update(source='pairs',pars={'runtime' : \
                            int((c4gli_afternoon_for_dump.pars.datetime_daylight - 
                                 c4gli.pars.datetime_daylight).total_seconds())})
    
    
                        print('ALMOST...')
                        if c4gli.pars.runtime > 18000.: # more than 5 hours simulation
                                
        
                            c4gli.get_global_input(globaldata)
                            print('VERY CLOSE...')
                            if c4gli.check_source_globaldata() and \
                                (c4gli.check_source(source='wyoming',\
                                                   check_only_sections='pars')):
                                c4gli.dump(fileout)
                                
                                c4gli_afternoon_for_dump.dump(fileout_afternoon)
                                
                                
                                # for keyEXP,dictEXP in experiments.items():
                                #     
                                #     c4gli.update(source=keyEXP,pars = dictEXP)
                                #     c4gl = class4gl(c4gli)
                                #     # c4gl.run()
                                #     
                                #     c4gl.dump(c4glfiles[key])
                                
                                print('HIT!!!')
                                one_run = True
                except:
                    print('get profile failed')
                
    if one_run:
        #STN.name = STN.name
        all_records_morning = get_records(pd.DataFrame([STN]),\
                                      args.path_output,\
                                      subset='morning',
                                      refetch_records=True,
                                      )
        all_records_afternoon = get_records(pd.DataFrame([STN]),\
                                      args.path_output,\
                                      subset='afternoon',
                                      refetch_records=True,
                                      )
    else:
        os.system('rm '+fnout)
        os.system('rm '+fnout_afternoon)

    # for c4glfile in c4glfiles:
    #     c4glfile.close()            

