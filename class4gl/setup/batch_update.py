# -*- coding: utf-8 -*-

"""
Usage:
python batch_update.py --exec $CLASS4GL/simulations/update_yaml_old.py
--path_experiments $VSC_DATA_VO/D2D/data/C4GL/GLOBAL_NOAC/ --path_input
$VSC_DATA_VO/D2D/data/C4GL/GLOBAL_NOAC_BACKUP_20180904/ --c4gl_path_lib
$CLASS4GL --split_by 50 --global_keys "KGC" --subset_input morning --experiments
"GLOBAL_NOAC"
"""


import pandas as pd
import io
import os
import numpy as np
import datetime as dt
import sys
import pytz
import math

import argparse

parser = argparse.ArgumentParser()
#if __name__ == '__main__':
parser.add_argument('--exec') # chunk simulation script
parser.add_argument('--first_station_row')
parser.add_argument('--last_station_row')
parser.add_argument('--pbs_string',default=' -l walltime=2:0:0')
parser.add_argument('--station_id') # run a specific station id
parser.add_argument('--error_handling')
parser.add_argument('--subset_input',default='morning') 
                                        # this tells which yaml subset
                                        # to initialize with.
                                        # Most common options are
                                        # 'morning' and 'ini'.
parser.add_argument('--subset_output',default='morning') 

# Tuntime is usually specified from the afternoon profile. You can also just
# specify the simulation length in seconds
parser.add_argument('--runtime')
# delete folders of experiments before running them
parser.add_argument('--split_by',default=50)# station soundings are split
                                            # up in chunks

parser.add_argument('--c4gl_path_lib')#,default='/user/data/gent/gvo000/gvo00090/D2D/software/CLASS/class4gl/lib')
parser.add_argument('--path_input') #,default='/user/data/gent/gvo000/gvo00090/D2D/data/SOUNDINGS/')
parser.add_argument('--path_output') #,default='/user/data/gent/gvo000/gvo00090/D2D/data/C4GL/')


#arguments only used for update_yaml.py
parser.add_argument('--path_dataset') 
parser.add_argument('--global_keys') 
parser.add_argument('--updates') 
args = parser.parse_args()

sys.path.insert(0, args.c4gl_path_lib)
from class4gl import class4gl_input, data_global,class4gl
from interface_multi import stations,stations_iterator, records_iterator,get_record_yaml,get_records
from class4gl import blh,class4gl_input


# this is a variant of global run in which the output of runs are still written
# out even when the run crashes.

# #only include the following timeseries in the model output
# timeseries_only = \
# ['Cm', 'Cs', 'G', 'H', 'L', 'LE', 'LEpot', 'LEref', 'LEsoil', 'LEveg', 'Lwin',
#  'Lwout', 'Q', 'RH_h', 'Rib', 'Swin', 'Swout', 'T2m', 'dq', 'dtheta',
#  'dthetav', 'du', 'dv', 'esat', 'gammaq', 'gammatheta', 'h', 'q', 'qsat',
#  'qsurf', 'ra', 'rs', 'theta', 'thetav', 'time', 'u', 'u2m', 'ustar', 'uw',
#  'v', 'v2m', 'vw', 'wq', 'wtheta', 'wthetae', 'wthetav', 'wthetae', 'zlcl']

print("getting all stations from --path_input")
# these are all the stations that are found in the input dataset
all_stations = stations(args.path_input,suffix=args.subset_input,refetch_stations=False)

print('defining all_stations_select')
# these are all the stations that are supposed to run by the whole batch (all
# chunks). We narrow it down according to the station(s) specified.
if args.station_id is not None:
    print("Selecting stations by --station_id")
    stations_iter = stations_iterator(all_stations)
    STNID,run_station = stations_iter.set_STNID(STNID=int(args.station_id))
    all_stations_select = pd.DataFrame([run_station])
else:
    print("Selecting stations from a row range in the table [--first_station_row,--last_station_row]")
    all_stations_select = pd.DataFrame(all_stations.table)
    if args.last_station_row is not None:
        all_stations_select = all_station_select.iloc[:(int(args.last_station)+1)]
    if args.first_station_row is not None:
        all_stations_select = all_station_select.iloc[int(args.first_station):]
print("station numbers included in the whole batch "+\
      "(all chunks):",list(all_stations_select.index))

print("getting all records of the whole batch")
all_records_morning_select = get_records(all_stations_select,\
                                         args.path_input,\
                                         subset=args.subset_input,\
                                         refetch_records=False,\
                                        )

print('splitting batch in --split_by='+args.split_by+' jobs.')
totalchunks = 0
for istation,current_station in all_stations_select.iterrows():
    records_morning_station_select = all_records_morning_select.query('STNID == '+str(current_station.name))
    chunks_current_station = len(records_morning_station_select.query('STNID == '+str(current_station.name)).chunk.unique())
    totalchunks +=chunks_current_station

print('total chunks (= size of array-job) per experiment: ' + str(totalchunks))

#if sys.argv[1] == 'qsub':
# with qsub



#C4GLJOB_timestamp="+dt.datetime.now().isoformat()+",
command = 'qsub '+args.pbs_string+' '+args.c4gl_path_lib+'/simulations/batch_simulations.pbs -t 0-'+\
            str(totalchunks-1)+" -v '"
# propagate arguments towards the job script
first = True
for argkey in args.__dict__.keys():
    if ((argkey not in ['pbs_string']) and \
        # default values are specified in the simulation script, so
        # excluded here
        (args.__dict__[argkey] is not None)
       ):
        if first:
            command +='C4GLJOB_'+argkey+'='+args.__dict__[argkey]
        else:
            command +=',C4GLJOB_'+argkey+'='+args.__dict__[argkey]
    first = False

command = command+"'"
print('Submitting array job: '+command)
os.system(command)


    #os.system(command)
# elif sys.argv[1] == 'wsub':
#     
#     # with wsub
#     STNlist = list(df_stations.iterrows())
#     NUMSTNS = len(STNlist)
#     PROCS = NUMSTNS 
#     BATCHSIZE = 1 #math.ceil(np.float(NUMSTNS)/np.float(PROCS))
# 
#     os.system('wsub -batch /user/data/gent/gvo000/gvo00090/D2D/scripts/C4GL/global_run.pbs -t 0-'+str(PROCS-1))

