

# -*- coding: utf-8 -*-

import logging
import pandas as pd
import io
import os
import numpy as np
import datetime as dt
import sys
import pytz
import math
import importlib
spam_loader = importlib.find_loader('Pysolar')
found = spam_loader is not None
if found:
    import Pysolar
    import Pysolar.util.GetSunriseSunset
else:
    import pysolar as Pysolar
    GetSunriseSunset =  Pysolar.util.get_sunrise_sunset

import argparse

#if __name__ == '__main__':
parser = argparse.ArgumentParser()
#parser.add_argument('--timestamp')
parser.add_argument('--path_forcing')#,default='/user/data/gent/gvo000/gvo00090/D2D/data/SOUNDINGS/')
parser.add_argument('--path_experiments')#,default='/user/data/gent/gvo000/gvo00090/D2D/data/C4GL/')
parser.add_argument('--first_YYYYMMDD',default="19810101")
parser.add_argument('--last_YYYYMMDD',default="20180101")
parser.add_argument('--first_station_row')
parser.add_argument('--last_station_row')
parser.add_argument('--station_id') # run a specific station id
parser.add_argument('--latitude') # run a specific station id
parser.add_argument('--longitude') # run a specific station id
parser.add_argument('--error_handling',default='dump_on_success')
parser.add_argument('--subset_forcing',default='morning') # this tells which yaml subset
parser.add_argument('--subset_experiments',default='ini') # this tells which yaml subset
                                                      # to initialize with.
                                                      # Most common options are
                                                      # 'morning' and 'ini'.

# Tuntime is usually specified from the afternoon profile. You can also just
# specify the simulation length in seconds
parser.add_argument('--runtime',default='from_afternoon_profile')

parser.add_argument('--experiments')
parser.add_argument('--split_by',default=-1)# station soundings are split
                                            # up in chunks

#parser.add_argument('--station-chunk',default=0)
parser.add_argument('--c4gl_path_lib')#,default='/user/data/gent/gvo000/gvo00090/D2D/software/CLASS/class4gl/lib')
parser.add_argument('--global_chunk_number') # this is the batch number according to split-by in case of considering all stations
parser.add_argument('--station_chunk_number') # this is the batch number according to split-by in case of considering all stations
args = parser.parse_args()

sys.path.insert(0, args.c4gl_path_lib)
from class4gl import class4gl_input, data_global,class4gl
from interface_multi import stations,stations_iterator, records_iterator,get_record_yaml,get_records
from class4gl import blh,class4gl_input

EXP_DEFS  =\
{
  'ERA-INTERIM_NOAC':    {'sw_ac' : [],'sw_ap': True,'sw_lit': False},
  'ERA-INTERIM_ADV':{'sw_ac' : ['adv',],'sw_ap': True,'sw_lit': False},
  'ERA-INTERIM_W':  {'sw_ac' : ['w',],'sw_ap': True,'sw_lit': False},
  'ERA-INTERIM_AC': {'sw_ac' : ['adv','w'],'sw_ap': True,'sw_lit': False},
  'NOAC':    {'sw_ac' : [],'sw_ap': True,'sw_lit': False},
  'ADV':{'sw_ac' : ['adv',],'sw_ap': True,'sw_lit': False},
  'W':  {'sw_ac' : ['w',],'sw_ap': True,'sw_lit': False},
  'AC': {'sw_ac' : ['adv','w'],'sw_ap': True,'sw_lit': False},
  'GLOBAL_NOAC':    {'sw_ac' : [],'sw_ap': True,'sw_lit': False},
  'GLOBAL_ADV':{'sw_ac' : ['adv',],'sw_ap': True,'sw_lit': False},
  'GLOBAL_W':  {'sw_ac' : ['w',],'sw_ap': True,'sw_lit': False},
  'GLOBAL_AC': {'sw_ac' : ['adv','w'],'sw_ap': True,'sw_lit': False},
  'IOPS_NOAC':    {'sw_ac' : [],'sw_ap': True,'sw_lit': False},
  'IOPS_ADV':{'sw_ac' : ['adv',],'sw_ap': True,'sw_lit': False},
  'IOPS_W':  {'sw_ac' : ['w',],'sw_ap': True,'sw_lit': False},
  'IOPS_AC': {'sw_ac' : ['adv','w'],'sw_ap': True,'sw_lit': False},
}


# iniitialize global data
# ===============================
print("Initializing global data")
# ===============================
globaldata = data_global()
globaldata.sources = {**globaldata.sources,**{
        "ERAINT:t"     : "/user/data/gent/gvo000/gvo00090/EXT/data/ERA-INTERIM/by_var_nc/t_6hourly/t_*_6hourly.nc",
        "ERAINT:q"     : "/user/data/gent/gvo000/gvo00090/EXT/data/ERA-INTERIM/by_var_nc/q_6hourly/q_*_6hourly.nc",
        "ERAINT:u"     : "/user/data/gent/gvo000/gvo00090/EXT/data/ERA-INTERIM/by_var_nc/u_6hourly/u_*_6hourly.nc",
        "ERAINT:v"     : "/user/data/gent/gvo000/gvo00090/EXT/data/ERA-INTERIM/by_var_nc/v_6hourly/v_*_6hourly.nc",
        }}

# ...  and load initial data pages
globaldata.load_datasets(recalc=0)



# ===============================
print("getting a list of stations")
# ===============================
all_stations = stations(args.path_forcing,suffix=args.subset_forcing,refetch_stations=False)


# # ===============================
# print("Selecting station by ID")
# # ===============================
# stations_iter = stations_iterator(all_stations)
# STNID,run_station = stations_iter.set_STNID(STNID=int(args.station_id))
# all_stations_select = pd.DataFrame([run_station])
# print(run_station)


# ====================================
print('defining all_stations_select')
# ====================================

# these are all the stations that are supposed to run by the whole batch (all
# chunks). We narrow it down according to the station(s) specified.
if (args.latitude is not None) or (args.longitude is not None):
    print('custom coordinates not implemented yet, please ask developer.')
elif args.station_id is not None:
    print("Selecting station by ID")
    stations_iter = stations_iterator(all_stations)
    STNID,run_station = stations_iter.set_STNID(STNID=int(args.station_id))
    all_stations_select = pd.DataFrame([run_station])
#     print("making a custom station according to the coordinates")
# 
#     STNID = 43.23
else:
     print("Selecting stations from a row range in the table")
     all_stations_select = pd.DataFrame(all_stations.table)
     if args.last_station_row is not None:
         all_stations_select = all_station_select.iloc[:(int(args.last_station)+1)]
     if args.first_station_row is not None:
         all_stations_select = all_station_select.iloc[int(args.first_station):]

print("station numbers included in the whole batch "+\
      "(all chunks):",list(all_stations_select.index))

dtfirst = dt.datetime.strptime(args.first_YYYYMMDD,"%Y%m%d",)
dtlast = dt.datetime.strptime(args.last_YYYYMMDD,"%Y%m%d",)
# ===============================
print("Creating daily timeseries from", dtfirst," to ", dtlast)
# ===============================
DTS = [dtfirst + dt.timedelta(days=iday) for iday in \
       range(int((dtlast + dt.timedelta(days=1) -
                  dtfirst).total_seconds()/3600./24.))]

if args.split_by != -1:
    totalchunks = len(all_stations_select)*len(DTS)/int(args.split_by)
else:
    totalchunks = len(all_stations_select)


if args.global_chunk_number is not None:
    run_station_chunk = np.mod(int(args.global_chunk_number),len(DTS)/int(args.split_by))
else:
    if args.station_chunk_number is not None:
        run_station_chunk = int(args.station_chunk_number)
    else:
        if args.split_by != -1:
            raise ValueError("Chunks are defined by --split-by, but I don't know which chunk to run. Please provide --global_chunk_number or --station_chunk_number, or leave out --split-by.")
        run_station_chunk = 0
        print("stations that are processed.",list(run_stations.index))

DTS_chunk = DTS[(int(run_station_chunk)*int(args.split_by)):\
                 (int(run_station_chunk)+1)*int(args.split_by)]

# for the current implementation we only consider one station. Let's upgrade it
# later for more stations.
run_station_chunk = int(args.global_chunk_number)

# ===============================
print('start looping over chunk')
# ===============================

os.system('mkdir -p '+args.path_experiments)


fn_ini = args.path_experiments+'/'+format(run_station.name,'05d')+'_'+\
        str(int(run_station_chunk))+'_'+args.subset_experiments+'.yaml'
file_ini = open(fn_ini,'w');print('Writing to: ',fn_ini)

for iDT,DT in enumerate(DTS_chunk):
    print(iDT,DT)
    c4gli = class4gl_input(debug_level=logging.INFO)
    c4gli.update(source='STNID'+format(STNID,'05d'),\
                 pars=dict(latitude  = float(run_station.latitude), \
                           longitude = float(run_station.longitude),\
                           lat       = float(run_station.latitude), \
                           # Note the difference between longitude and lon. The
                           # lon variable should always be zero because we are
                           # always working in solar time for running CLASS
                           lon       = 0.,\
                           STNID     = int(STNID)))

    lSunrise, lSunset = GetSunriseSunset(c4gli.pars.latitude,0.,DT)

    #start simulation at sunrise and stop at one hour before sunset
    runtime = (lSunset - lSunrise).total_seconds() - 3600.*1.
    ldatetime = lSunrise
    datetime = ldatetime - dt.timedelta(hours=c4gli.pars.longitude/360.*24.)
    datetime_daylight = datetime
    c4gli.update(source='timeseries',   \
                 pars=dict(\
                           lSunrise = lSunrise, \
                           lSunset = lSunset, \
                           datetime = datetime, \
                           ldatetime = ldatetime, \
                           ldatetime_daylight = ldatetime, \
                           datetime_daylight = datetime, \
                           doy = datetime.timetuple().tm_yday,\
                           runtime = runtime,\
                          ))

    c4gli.get_global_input(globaldata)

    c4gli.update(source='era-interim',pars={'Ps' : c4gli.pars.sp})

    cp         = 1005.                 # specific heat of dry air [J kg-1 K-1]
    Rd         = 287.                  # gas constant for dry air [J kg-1 K-1]
    Rv         = 461.5                 # gas constant for moist air [J kg-1 K-1]
    R = (Rd*(1.-c4gli.air_ac.q) + Rv*c4gli.air_ac.q)
    rho = c4gli.air_ac.p/R/c4gli.air_ac.t
    dz = c4gli.air_ac.delpdgrav/rho
    z = [dz.iloc[-1]/2.]
    for idz in list(reversed(range(0,len(dz)-1,1))):
        z.append(z[-1]+ (dz[idz+1]+dz[idz])/2.)
    z = list(reversed(z))

    theta = c4gli.air_ac.t * \
               (c4gli.pars.sp/(c4gli.air_ac.p))**(R/cp)
    thetav   = theta*(1. + 0.61 * c4gli.air_ac.q)

    
    c4gli.update(source='era-interim',air_ac=pd.DataFrame({'z':list(z),
                                                           'theta':list(theta),
                                                           'thetav':list(thetav),
                                                          }))
    air_ap_input = c4gli.air_ac[::-1].reset_index().drop('index',axis=1)
    air_ap_mode = 'b'
    air_ap_input_source = c4gli.query_source('air_ac:theta')


    c4gli.mixed_layer_fit(air_ap=air_ap_input,
                         source=air_ap_input_source,
                         mode=air_ap_mode)

    if not c4gli.check_source_globaldata():
        print('Warning: some input sources appear invalid')

    c4gli.dump(file_ini)

file_ini.close()
all_records_morning = get_records(pd.DataFrame([run_station]),\
                              args.path_experiments,\
                              getchunk = int(run_station_chunk),\
                              subset=args.subset_experiments,
                              refetch_records=True,
                              )

