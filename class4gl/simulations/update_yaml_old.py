# -*- coding: utf-8 -*-

""" 
Purpose:
    update variables in class4gl yaml files, eg., when you need new categorical
    values in the table.


"""



import pandas as pd
import io
import os
import numpy as np
import datetime as dt
import sys
import pytz
import math
import dateutil.parser

import argparse


#if __name__ == '__main__':
parser = argparse.ArgumentParser()
parser.add_argument('--path_forcing')#,default='/user/data/gent/gvo000/gvo00090/D2D/data/C4GL/')
parser.add_argument('--first_station_row')
parser.add_argument('--last_station_row')
parser.add_argument('--path_experiments')
parser.add_argument('--experiments')
parser.add_argument('--station_id') # run a specific station id
parser.add_argument('--mode',default='ini') # this tells which yaml subset
parser.add_argument('--subset_forcing',default='morning') # this tells which yaml subset
                                                      # to update in the yaml
                                                      # dataset.
                                                      # Most common options are
                                                      # 'morning' and 'ini'.

parser.add_argument('--split_by',default=-1)# station soundings are split

#parser.add_argument('--station-chunk',default=0)
parser.add_argument('--c4gl_path_lib')#,default='/user/data/gent/gvo000/gvo00090/D2D/software/CLASS/class4gl/lib')
parser.add_argument('--global_chunk_number') # this is the batch number according to split-by in case of considering all stations
parser.add_argument('--station_chunk_number') # this is the batch number according to split-by in case of considering all stations
parser.add_argument('--global_keys') 
args = parser.parse_args()

sys.path.insert(0, args.c4gl_path_lib)
from class4gl import class4gl_input, data_global,class4gl
from interface_multi import stations,stations_iterator, records_iterator,get_record_yaml,get_records
from class4gl import blh,class4gl_input

# iniitialize global data
globaldata = data_global()
# ...  and load initial data pages
globaldata.load_datasets(recalc=0)


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


    # if not (int(args.split_by) > 0) :
    #         raise ValueError("global_chunk_number is specified, but --split-by is not a strict positive number, so I don't know how to split the batch into chunks.")

    run_station_chunk = None
    print('determining the station and its chunk number according global_chunk_number ('+args.global_chunk_number+')')
    totalchunks = 0
    stations_iter = all_stations_select.iterrows()
    in_current_chunk = False
    try:
        while not in_current_chunk:
            istation,current_station = stations_iter.__next__()
            all_records_morning_station_select = all_records_morning_select.query('STNID == '+str(current_station.name))
            #chunks_current_station = math.ceil(float(len(all_records_morning_station_select))/float(args.split_by))

            chunks_current_station = len(all_records_morning_station_select.query('STNID == '+str(current_station.name)).chunk.unique())
            print('chunks_current_station',chunks_current_station)

            in_current_chunk = (int(args.global_chunk_number) < (totalchunks+chunks_current_station))
        
            if in_current_chunk:
                run_stations = pd.DataFrame([current_station])# run_stations.loc[(int(args.__dict__['last_station'])]
                run_station_chunk =all_records_morning_station_select.query('STNID == '+str(current_station.name)).chunk.unique()[int(args.global_chunk_number) - totalchunks ]
        
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
print('Fetching current records')
records_forcing = get_records(run_stations,\
                              args.path_forcing,\
                              subset=args.subset_forcing,
                              refetch_records=False,
                              )

# if args.timestamp is None:
#     backupdir = args.path_forcing+'/'+dt.datetime.now().isoformat()+'/'
# else: 
#     backupdir = args.path_forcing+'/'+args.timestamp+'/'
# print('creating backup dir: '+backupdir)
# os.system('mkdir -p "'+backupdir+'"')


for EXP in args.experiments.strip().split(" "):
    os.system('mkdir -p '+args.path_experiments+'/'+EXP+'/')
    for istation,current_station in run_stations.iterrows():
        records_forcing_station = records_forcing.query('STNID == ' +\
                                                        str(current_station.name))
    
        records_forcing_station_chunk = records_forcing.query('STNID == ' +\
                                                        str(current_station.name)+\
                                                       '& chunk == '+str(run_station_chunk))
        print('lenrecords_forcing_station: ',len(records_forcing_station))
        print('split_by*run_station_chunk',int(args.split_by) * int(run_station_chunk))
        print('split_by*run_station_chunk+1',int(args.split_by) * int(run_station_chunk+1))
        
        # if (int(args.split_by) * int(run_station_chunk)) >= (len(records_forcing_station)):
        #     print("warning: outside of profile number range for station "+\
        #           str(current_station)+". Skipping chunk number for this station.")
        if len(records_forcing_station_chunk) == 0:
            print("warning: outside of profile number range for station "+\
                  str(current_station)+". Skipping chunk number for this station.")
        else:
            # normal case
            if ((int(args.split_by) > 0) or \
                (os.path.isfile(args.path_forcing+'/'+format(current_station.name,'05d')+'_'+\
                     str(run_station_chunk)+'_'+args.subset_forcing+'.yaml'))):
                fn_forcing = \
                        args.path_forcing+'/'+format(current_station.name,'05d')+'_'+\
                        str(run_station_chunk)+'_'+args.subset_forcing+'.yaml'
                file_forcing = \
                    open(fn_forcing,'r')
                fn_experiment = args.path_experiments+'/'+EXP+'/'+format(current_station.name,'05d')+'_'+\
                         str(run_station_chunk)+'_'+args.subset_forcing+'.yaml'
                file_experiment = \
                    open(fn_experiment,'w')
                fn_forcing_pkl = args.path_forcing+'/'+format(current_station.name,'05d')+'_'+\
                         str(run_station_chunk)+'_'+args.subset_forcing+'.pkl'
    
                # fn_backup = backupdir+format(current_station.name,'05d')+'_'+\
                #          str(run_station_chunk)+'_'+args.subset_forcing+'.yaml'
                # fn_backup_pkl = backupdir+format(current_station.name,'05d')+'_'+\
                #          str(run_station_chunk)+'_'+args.subset_forcing+'.pkl'
            else:
                print("\
    Warning. We are choosing chunk 0 without specifying it in filename.    \
     No-chunk naming will be removed in the future."\
                     )
    
                fn_forcing = \
                        args.path_forcing+'/'+format(current_station.name,'05d')+'_'+\
                        args.subset_forcing+'.yaml'
                file_forcing = \
                    open(fn_forcing,'r')
                fn_experiment = args.path_experiments+'/'+EXP+'/'+format(current_station.name,'05d')+'_'+\
                         str(run_station_chunk)+'_'+args.subset_forcing+'.yaml'
                file_experiment = \
                    open(fn_experiment,'w')
                fn_forcing_pkl = args.path_forcing+format(current_station.name,'05d')+'_'+\
                         str(run_station_chunk)+'_'+args.subset_forcing+'.pkl'
    
                # fn_backup = backupdir+format(current_station.name,'05d')+'_'+\
                #          str(run_station_chunk)+'_'+args.subset_forcing+'.yaml'
                # fn_backup_pkl = backupdir+format(current_station.name,'05d')+'_'+\
                #          args.subset_forcing+'.pkl'
    
            onerun = False
            print('starting station chunk number: '\
                  +str(run_station_chunk)+'(size: '+str(args.split_by)+' soundings)')
    
            #records_forcing_station_chunk = records_forcing_station[(int(args.split_by)*run_station_chunk):(int(args.split_by)*(run_station_chunk+1))]
    
            # records_forcing_station_chunk = records_forcing.query('STNID == ' +\
            #                                                 str(current_station.name)+\
            #                                                '& chunk == '+str(run_station_chunk))
            isim = 0
            for (STNID,chunk,index),record_forcing in records_forcing_station_chunk.iterrows():
                    print('starting '+str(isim+1)+' out of '+\
                      str(len(records_forcing_station_chunk) )+\
                      ' (station total: ',str(len(records_forcing_station)),')')  
                
                    c4gli_forcing = get_record_yaml(file_forcing, 
                                                    record_forcing.index_start, 
                                                    record_forcing.index_end,
                                                    mode=args.mode)
                    seltropo = (c4gli_forcing.air_ac.p > c4gli_forcing.air_ac.p.iloc[-1]+ 3000.*(- 1.2 * 9.81 ))
                    profile_tropo = c4gli_forcing.air_ac[seltropo]
                    mean_advt_tropo = np.mean(profile_tropo.advt_x +profile_tropo.advt_y )
                    c4gli_forcing.update(source='era-interim',pars={'advt_tropo':mean_advt_tropo})
                    
                    #print('c4gli_forcing_ldatetime',c4gli_forcing.pars.ldatetime)
                    
                    if args.global_keys is not None:
                        print(args.global_keys.strip(' ').split(' '))
                        c4gli_forcing.get_global_input(
                            globaldata, 
                            only_keys=args.global_keys.strip(' ').split(' ')
                        )
    
                    c4gli_forcing.dump(file_experiment)
                        
                        
                    onerun = True
                    isim += 1
    
    
            file_forcing.close()
            file_experiment.close()
    
            if onerun:
                # os.system('mv "'+fn_forcing+'" "'+fn_backup+'"')
                # if os.path.isfile(fn_forcing_pkl):
                #     os.system('mv "'+fn_forcing_pkl+'" "'+fn_backup_pkl+'"')
                # os.system('mv "'+fn_experiment+'" "'+fn_forcing+'"')
                # print('mv "'+fn_experiment+'" "'+fn_forcing+'"')
                records_forcing_current_cache = get_records(pd.DataFrame([current_station]),\
                                                           args.path_experiments+'/'+EXP+'/',\
                                                           getchunk = int(run_station_chunk),\
                                                           subset=args.subset_forcing,
                                                           refetch_records=True,
                                                           )
    
