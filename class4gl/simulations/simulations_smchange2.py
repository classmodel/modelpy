# -*- coding: utf-8 -*-

import pandas as pd
import io
import os
import numpy as np
import datetime as dt
import sys
import pytz
import math

import argparse

#if __name__ == '__main__':
parser = argparse.ArgumentParser()
#parser.add_argument('--timestamp')
parser.add_argument('--path_forcing')#,default='/user/data/gent/gvo000/gvo00090/D2D/data/SOUNDINGS/')
parser.add_argument('--path_experiments')#,default='/user/data/gent/gvo000/gvo00090/D2D/data/C4GL/')
parser.add_argument('--first_station_row')
parser.add_argument('--last_station_row')
parser.add_argument('--station_id') # run a specific station id
parser.add_argument('--error_handling',default='dump_on_success')
parser.add_argument('--subset_forcing',default='morning') # this tells which yaml subset
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

# this is a variant of global run in which the output of runs are still written
# out even when the run crashes.

# #only include the following timeseries in the model output
# timeseries_only = \
# ['Cm', 'Cs', 'G', 'H', 'L', 'LE', 'LEpot', 'LEref', 'LEsoil', 'LEveg', 'Lwin',
#  'Lwout', 'Q', 'RH_h', 'Rib', 'Swin', 'Swout', 'T2m', 'dq', 'dtheta',
#  'dthetav', 'du', 'dv', 'esat', 'gammaq', 'gammatheta', 'h', 'q', 'qsat',
#  'qsurf', 'ra', 'rs', 'theta', 'thetav', 'time', 'u', 'u2m', 'ustar', 'uw',
#  'v', 'v2m', 'vw', 'wq', 'wtheta', 'wthetae', 'wthetav', 'wthetae', 'zlcl']


EXP_DEFS  =\
{
  'SMa':    {'sw_ac' : [],'sw_ap': True,'sw_lit': False, 'wgchange': 0.05},
  'SMb':    {'sw_ac' : [],'sw_ap': True,'sw_lit': False, 'wgchange': 0.1},
  'SMc':    {'sw_ac' : [],'sw_ap': True,'sw_lit': False, 'wgchange': 0.15},
  'SMd':    {'sw_ac' : [],'sw_ap': True,'sw_lit': False, 'wgchange': 0.02},
  'SMe':    {'sw_ac' : [],'sw_ap': True,'sw_lit': False, 'wgchange': 0.25},
  'SMf':    {'sw_ac' : [],'sw_ap': True,'sw_lit': False, 'wgchange': 0.3},
  'SMg':    {'sw_ac' : [],'sw_ap': True,'sw_lit': False, 'wgchange': 0.35},
  'SMh':    {'sw_ac' : [],'sw_ap': True,'sw_lit': False, 'wgchange': 0.4},
  'SMi':    {'sw_ac' : [],'sw_ap': True,'sw_lit': False, 'wgchange': 0.45},
  'SMj':    {'sw_ac' : [],'sw_ap': True,'sw_lit': False, 'wgchange': 0.5},
  'SMk':    {'sw_ac' : [],'sw_ap': True,'sw_lit': False, 'wgchange': 0.55},
  'SMl':    {'sw_ac' : [],'sw_ap': True,'sw_lit': False, 'wgchange': 0.6},
  'SMm':    {'sw_ac' : [],'sw_ap': True,'sw_lit': False, 'wgchange': 0.65},
  'SMo':    {'sw_ac' : [],'sw_ap': True,'sw_lit': False, 'wgchange': 0.7},
  'SMp':    {'sw_ac' : [],'sw_ap': True,'sw_lit': False, 'wgchange': 0.75},
  'SMq':    {'sw_ac' : [],'sw_ap': True,'sw_lit': False, 'wgchange': 0.08},
  'SMr':    {'sw_ac' : [],'sw_ap': True,'sw_lit': False, 'wgchange': 0.85},
  'SMs':    {'sw_ac' : [],'sw_ap': True,'sw_lit': False, 'wgchange': 0.9},
  'SMt':    {'sw_ac' : [],'sw_ap': True,'sw_lit': False, 'wgchange': 0.95},
  'SMu':    {'sw_ac' : [],'sw_ap': True,'sw_lit': False, 'wgchange': 1.},
}

##for i in range(1,20):
##    EXP_DEFS['SM'+str(i)] = {'sw_ac' : [],'sw_ap': True,'sw_lit': False, 'wgchange':float(i)/20.}



# #SET = 'GLOBAL'
# SET = args.dataset


print("getting stations")
# these are all the stations that are found in the input dataset
all_stations = stations(args.path_forcing,suffix=args.subset_forcing,refetch_stations=False)

print('defining all_stations_select')
# these are all the stations that are supposed to run by the whole batch (all
# chunks). We narrow it down according to the station(s) specified.
if args.station_id is not None:
    print("Selecting station by ID")
    stations_iter = stations_iterator(all_stations)
    STNID,run_station = stations_iter.set_STNID(STNID=int(args.station_id))
    all_stations_select = pd.DataFrame([run_station])
else:
    print("Selecting stations from a row range in the table")
    all_stations_select = pd.DataFrame(all_stations.table)
    if args.last_station_row is not None:
        all_stations_select = all_station_select.iloc[:(int(args.last_station)+1)]
    if args.first_station_row is not None:
        all_stations_select = all_station_select.iloc[int(args.first_station):]
print("station numbers included in the whole batch "+\
      "(all chunks):",list(all_stations_select.index))

print("getting all records of the whole batch")
all_records_morning_select = get_records(all_stations_select,\
                                         args.path_forcing,\
                                         subset=args.subset_forcing,
                                         refetch_records=False,
                                         )

# only run a specific chunck from the selection
if args.global_chunk_number is not None:
    if args.station_chunk_number is not None:
        raise ValueError('You need to specify either global-chunk-number or station-chunk-number, not both.')


    if not (int(args.split_by) > 0) :
            raise ValueError("global_chunk_number is specified, but --split-by is not a strict positive number, so I don't know how to split the batch into chunks.")

    run_station_chunk = None
    print('determining the station and its chunk number according global_chunk_number ('+args.global_chunk_number+')')
    totalchunks = 0
    stations_iter = all_stations_select.iterrows()
    in_current_chunk = False
    try:
        while not in_current_chunk:
            istation,current_station = stations_iter.__next__()
            all_records_morning_station_select = all_records_morning_select.query('STNID == '+str(current_station.name))
            chunks_current_station = math.ceil(float(len(all_records_morning_station_select))/float(args.split_by))
            print('chunks_current_station',chunks_current_station)
            in_current_chunk = (int(args.global_chunk_number) < (totalchunks+chunks_current_station))
        
            if in_current_chunk:
                run_stations = pd.DataFrame([current_station])# run_stations.loc[(int(args.__dict__['last_station'])]
                run_station_chunk = int(args.global_chunk_number) - totalchunks 
        
            totalchunks +=chunks_current_station
        

    except StopIteration:
       raise ValueError("Could not determine station chunk number.  --global_chunk_number ("+args.global_chunk_number+") outside of range [0,"+ str(totalchunks)+'[')
    print("station = ",list(run_stations.index))
    print("station chunk number:",run_station_chunk)

# if no global chunk is specified, then run the whole station selection in one run, or
# a specific chunk for each selected station according to # args.station_chunk_number
else:
    run_stations = pd.DataFrame(all_stations_select)# run_stations.loc[(int(args.__dict__['last_station'])]
    if args.station_chunk_number is not None:
        run_station_chunk = int(args.station_chunk_number)
        print("station(s) that is processed.",list(run_stations.index))
        print("chunk number: ",run_station_chunk)
    else:
        if args.split_by != -1:
            raise ValueError("Chunks are defined by --split-by, but I don't know which chunk to run. Please provide --global_chunk_number or --station_chunk_number, or leave out --split-by.")
        run_station_chunk = 0
        print("stations that are processed.",list(run_stations.index))
        

#print(all_stations)
print('Fetching initial/forcing records')
records_morning = get_records(run_stations,\
                              args.path_forcing,\
                              subset=args.subset_forcing,
                              refetch_records=False,
                              )

# note that if runtime is an integer number, we don't need to get the afternoon
# profiles. 
if args.runtime == 'from_afternoon_profile':
    print('Fetching afternoon records for determining the simulation runtimes')
    records_afternoon = get_records(run_stations,\
                                    args.path_forcing,\
                                    subset='afternoon',
                                    refetch_records=False,
                                    )
    
    # print(records_morning.index)
    # print(records_afternoon.index)
    # align afternoon records with the noon records, and set same index
    print('hello')
    print(len(records_afternoon))
    print(len(records_morning))

    print("aligning morning and afternoon records")
    records_morning['dates'] = records_morning.ldatetime.dt.date
    records_afternoon['dates'] = records_afternoon.ldatetime.dt.date
    records_afternoon.set_index(['STNID','dates'],inplace=True)
    ini_index_dates = records_morning.set_index(['STNID','dates']).index
    records_afternoon = records_afternoon.loc[ini_index_dates]
    records_afternoon.index = records_morning.index

experiments = args.experiments.strip(' ').split(' ')
for expname in experiments:
    exp = EXP_DEFS[expname]
    path_exp = args.path_experiments+'/'+expname+'/'

    os.system('mkdir -p '+path_exp)
    for istation,current_station in run_stations.iterrows():
        records_morning_station = records_morning.query('STNID == '+str(current_station.name))
        if (int(args.split_by) * int(run_station_chunk)) >= (len(records_morning_station)):
            print("warning: outside of profile number range for station "+\
                  str(current_station)+". Skipping chunk number for this station.")
        else:
            file_morning = open(args.path_forcing+'/'+format(current_station.name,'05d')+'_morning.yaml')
            file_afternoon = open(args.path_forcing+'/'+format(current_station.name,'05d')+'_afternoon.yaml')
            fn_ini = path_exp+'/'+format(current_station.name,'05d')+'_'+\
                     str(int(run_station_chunk))+'_ini.yaml'
            fn_mod = path_exp+'/'+format(current_station.name,'05d')+'_'+\
                     str(int(run_station_chunk))+'_mod.yaml'
            file_ini = open(fn_ini,'w')
            file_mod = open(fn_mod,'w')

            #iexp = 0
            onerun = False
            print('starting station chunk number: '\
                  +str(run_station_chunk)+'(size: '+str(args.split_by)+' soundings)')

            records_morning_station_chunk = records_morning_station[(int(args.split_by)*run_station_chunk):(int(args.split_by)*(run_station_chunk+1))]

            isim = 0
            for (STNID,chunk,index),record_morning in records_morning_station_chunk.iterrows():
                    print('starting '+str(isim+1)+' out of '+\
                      str(len(records_morning_station_chunk) )+\
                      ' (station total: ',str(len(records_morning_station)),')')  
                
            
                    c4gli_morning = get_record_yaml(file_morning, 
                                                    record_morning.index_start, 
                                                    record_morning.index_end,
                                                    mode='ini')
                    
                    #print('c4gli_morning_ldatetime',c4gli_morning.pars.ldatetime)

                    sm = exp['wgchange']*(c4gli_morning.pars.wfc - c4gli_morning.pars.wwilt) + c4gli_morning.pars.wwilt
                    c4gli_morning.update(source=expname, pars={'wg': sm})
                    c4gli_morning.update(source=expname, pars={'w2': sm})
                    
                    
                    if args.runtime == 'from_afternoon_profile':
                        record_afternoon = records_afternoon.loc[(STNID,chunk,index)]
                        c4gli_afternoon = get_record_yaml(file_afternoon, 
                                                          record_afternoon.index_start, 
                                                          record_afternoon.index_end,
                                                        mode='ini')
                        runtime = int((c4gli_afternoon.pars.datetime_daylight - 
                                             c4gli_morning.pars.datetime_daylight).total_seconds())
                    else:
                        runtime = int(args.runtime)

            
                    c4gli_morning.update(source='pairs',pars={'runtime' : \
                                        runtime})
                    c4gli_morning.update(source=expname, pars=exp)

                    c4gl = class4gl(c4gli_morning)

                    if args.error_handling == 'dump_always':
                        try:
                            c4gl.run()
                            print('run succesfull')
                        except:
                            print('run not succesfull')
                        onerun = True

                        c4gli_morning.dump(file_ini)
                        
                        
                        c4gl.dump(file_mod,\
                                  include_input=False,\
                                  #timeseries_only=timeseries_only,\
                                 )
                        onerun = True
                    # in this case, only the file will dumped if the runs were
                    # successful
                    elif args.error_handling == 'dump_on_success':
                        try:
                            c4gl.run()
                            print('run succesfull')
                            c4gli_morning.dump(file_ini)
                            
                            
                            c4gl.dump(file_mod,\
                                      include_input=False,\
                                      #timeseries_only=timeseries_only,\
                                     )
                            onerun = True
                        except:
                            print('run not succesfull')
                    isim += 1


            file_ini.close()
            file_mod.close()
            file_morning.close()
            file_afternoon.close()
    
            if onerun:
                records_ini = get_records(pd.DataFrame([current_station]),\
                                                           path_exp,\
                                                           getchunk = int(run_station_chunk),\
                                                           subset='ini',
                                                           refetch_records=True,
                                                           )
                records_mod = get_records(pd.DataFrame([current_station]),\
                                                           path_exp,\
                                                           getchunk = int(run_station_chunk),\
                                                           subset='mod',\
                                                           refetch_records=True,\
                                                           )
            else:
                # remove empty files
                os.system('rm '+fn_ini)
                os.system('rm '+fn_mod)
    
    # # align afternoon records with initial records, and set same index
    # records_afternoon.index = records_afternoon.ldatetime.dt.date
    # records_afternoon = records_afternoon.loc[records_ini.ldatetime.dt.date]
    # records_afternoon.index = records_ini.index
    
    # stations_for_iter = stations(path_exp)
    # for STNID,station in stations_iterator(stations_for_iter):
    #     records_current_station_index = \
    #             (records_ini.index.get_level_values('STNID') == STNID)
    #     file_current_station_mod = STNID
    # 
    #     with \
    #     open(path_exp+'/'+format(STNID,"05d")+'_ini.yaml','r') as file_station_ini, \
    #     open(path_exp+'/'+format(STNID,"05d")+'_mod.yaml','r') as file_station_mod, \
    #     open(path_forcing+'/'+format(STNID,"05d")+'_afternoon.yaml','r') as file_station_afternoon:
    #         for (STNID,index),record_ini in records_iterator(records_ini):
    #             c4gli_ini = get_record_yaml(file_station_ini, 
    #                                         record_ini.index_start, 
    #                                         record_ini.index_end,
    #                                         mode='ini')
    #             #print('c4gli_in_ldatetime 3',c4gli_ini.pars.ldatetime)
    # 
    #             record_mod = records_mod.loc[(STNID,index)]
    #             c4gl_mod = get_record_yaml(file_station_mod, 
    #                                         record_mod.index_start, 
    #                                         record_mod.index_end,
    #                                         mode='mod')
    #             record_afternoon = records_afternoon.loc[(STNID,index)]
    #             c4gl_afternoon = get_record_yaml(file_station_afternoon, 
    #                                         record_afternoon.index_start, 
    #                                         record_afternoon.index_end,
    #                                         mode='ini')

