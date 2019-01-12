
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
parser.add_argument('--exec') # chunk simulation script
parser.add_argument('--pbs_string',default='')#,default='/user/data/gent/gvo000/gvo00090/D2D/data/SOUNDINGS/')
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
    
        "ERAINT:t"     : "/user/data/gent/gvo000/gvo00090/EXT/data/ERA-INTERIM/by_var_nc/t_6hourly/t_19830609-19830808_6hourly.nc",
        "ERAINT:q"     : "/user/data/gent/gvo000/gvo00090/EXT/data/ERA-INTERIM/by_var_nc/q_6hourly/q_19830609-19830808_6hourly.nc",
        "ERAINT:u"     : "/user/data/gent/gvo000/gvo00090/EXT/data/ERA-INTERIM/by_var_nc/u_6hourly/u_19830609-19830808_6hourly.nc",
        "ERAINT:v"     : "/user/data/gent/gvo000/gvo00090/EXT/data/ERA-INTERIM/by_var_nc/v_6hourly/v_19830609-19830808_6hourly.nc",
    
#        "ERAINT:q"     : "/user/data/gent/gvo000/gvo00090/EXT/data/ERA-INTERIM/by_var_nc/q_6hourly/q_19830209-19830410_6hourly.nc",
 #       "ERAINT:q"     : "/user/data/gent/gvo000/gvo00090/EXT/data/ERA-INTERIM/by_var_nc/q_6hourly/q*_6hourly.nc",
 #       "ERAINT:u"     : "/user/data/gent/gvo000/gvo00090/EXT/data/ERA-INTERIM/by_var_nc/u_6hourly/u*_6hourly.nc",
 #       "ERAINT:v"     : "/user/data/gent/gvo000/gvo00090/EXT/data/ERA-INTERIM/by_var_nc/v_6hourly/v*_6hourly.nc",
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
    totalchunks = len(all_stations_select)*math.ceil(len(DTS)/int(args.split_by))
else:
    totalchunks = len(all_stations_select)

print(totalchunks)

#if args.cleanup_experiments:
#    os.system("rm -R "+args.path_experiments+'/')

# C4GLJOB_timestamp="+dt.datetime.now().isoformat()+",
command = 'qsub '+args.pbs_string+' '+args.c4gl_path_lib+'/setup/batch_setup_era.pbs -t 0-'+\
            str(totalchunks-1)+" -v "
# propagate arguments towards the job script
lfirst = True
for argkey in args.__dict__.keys():
    if ((argkey not in ['experiments','pbs_string','cleanup_experiments']) and \
        # default values are specified in the simulation script, so
        # excluded here
        (args.__dict__[argkey] is not None)
       ):
        if lfirst:
            command +=' C4GLJOB_'+argkey+'='+args.__dict__[argkey]
        else:
            command +=',C4GLJOB_'+argkey+'='+args.__dict__[argkey]
        lfirst=False

print('Submitting array job: '+command)
os.system(command)
