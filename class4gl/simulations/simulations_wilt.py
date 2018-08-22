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
parser.add_argument('--global-chunk') # this is the batch number according to split-by in case of considering all stations
parser.add_argument('--first-station-row')
parser.add_argument('--last-station-row')
parser.add_argument('--station-id') # run a specific station id
parser.add_argument('--path-experiments')#,default='/user/data/gent/gvo000/gvo00090/D2D/data/C4GL/')
parser.add_argument('--path-soundings')#,default='/user/data/gent/gvo000/gvo00090/D2D/data/SOUNDINGS/')
parser.add_argument('--error-handling',default='dump_on_success')
parser.add_argument('--experiments')
parser.add_argument('--split-by',default=-1)# station soundings are split
                                            # up in chunks

parser.add_argument('--station-chunk',default=0)
parser.add_argument('--c4gl-path-lib')#,default='/user/data/gent/gvo000/gvo00090/D2D/software/CLASS/class4gl/lib')
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
  'NOAC':    {'sw_ac' : [],'sw_ap': True,'sw_lit': False},
  'ADV':{'sw_ac' : ['adv',],'sw_ap': True,'sw_lit': False},
  'W':  {'sw_ac' : ['w',],'sw_ap': True,'sw_lit': False},
  'AC': {'sw_ac' : ['adv','w'],'sw_ap': True,'sw_lit': False},
  'WILT':    {'sw_ac' : [],'sw_ap': True,'sw_lit': False},
}


# #SET = 'GLOBAL'
# SET = args.dataset

# path_soundingsSET = args.path_soundings+'/'+SET+'/'

print("getting stations")
all_stations = stations(args.path_soundings,suffix='morning',refetch_stations=False)

if args.global_chunk is not None:
    
    all_records_morning = get_records(all_stations.table,\
                                  args.path_soundings,\
                                  subset='morning',
                                  refetch_records=False,
                                  )
    totalchunks = 0
    stations_iter = all_stations.table.iterrows()
    in_current_chunk = False
    while not in_current_chunk:
        istation,current_station = stations_iter.__next__()
        all_records_morning_station = all_records_morning.query('STNID == '+str(current_station.name))
        chunks_current_station = math.ceil(float(len(all_records_morning_station))/float(args.split_by))
        in_current_chunk = (int(args.global_chunk) < (totalchunks+chunks_current_station))

        if in_current_chunk:
            run_stations = pd.DataFrame([current_station])# run_stations.loc[(int(args.__dict__['last_station'])]
            run_station_chunk = int(args.global_chunk) - totalchunks 

        totalchunks +=chunks_current_station

else:
    if args.station_id is not None:
        print("Selecting station by ID")
        print(all_stations.table)
        stations_iter = stations_iterator(all_stations)
        STNID,run_station = stations_iter.set_STNID(STNID=int(args.station_id))
        run_stations = pd.DataFrame([run_station])
    else:
        print("Selecting stations from a row range in the table")
        run_stations = pd.DataFrame(all_stations.table)
        if args.last_station_row is not None:
            run_stations = run_stations.iloc[:(int(args.last_station)+1)]
        if args.first_station_row is not None:
            run_stations = run_stations.iloc[int(args.first_station):]
    run_station_chunk = args.station_chunk

#print(all_stations)
records_morning = get_records(run_stations,\
                              args.path_soundings,\
                              subset='morning',
                              refetch_records=False,
                              )
records_afternoon = get_records(run_stations,\
                                args.path_soundings,\
                                subset='afternoon',
                                refetch_records=False,
                                )

# print(records_morning.index)
# print(records_afternoon.index)
# align afternoon records with the noon records, and set same index
records_afternoon.index = records_afternoon.ldatetime.dt.date
records_afternoon = records_afternoon.loc[records_morning.ldatetime.dt.date]
records_afternoon.index = records_morning.index

experiments = args.experiments.split(';')
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
            file_morning = open(args.path_soundings+'/'+format(current_station.name,'05d')+'_morning.yaml')
            file_afternoon = open(args.path_soundings+'/'+format(current_station.name,'05d')+'_afternoon.yaml')
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
                    print('starting '+str(isim)+' out of '+\
                      str(len(records_morning_station_chunk) )+\
                      ' (station total: ',str(len(records_morning_station)),')')  
                
            
                    c4gli_morning = get_record_yaml(file_morning, 
                                                    record_morning.index_start, 
                                                    record_morning.index_end,
                                                    mode='ini')
                    
                    #print('c4gli_morning_ldatetime',c4gli_morning.pars.ldatetime)
                    
                    
                    record_afternoon = records_afternoon.loc[(STNID,chunk,index)]
                    c4gli_afternoon = get_record_yaml(file_afternoon, 
                                                      record_afternoon.index_start, 
                                                      record_afternoon.index_end,
                                                    mode='ini')
            
                    c4gli_morning.update(source='pairs',pars={'runtime' : \
                                        int((c4gli_afternoon.pars.datetime_daylight - 
                                             c4gli_morning.pars.datetime_daylight).total_seconds())})
                    c4gli_morning.update(source=expname, pars=exp)

                    c4gli_morning.update(source=expname, \
                                         pars={'wg':c4gli_morning.pars.wwilt,\
                                               'w2':c4gli_morning.pars.wwilt},
                                        )
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
    #     open(path_soundings+'/'+format(STNID,"05d")+'_afternoon.yaml','r') as file_station_afternoon:
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

