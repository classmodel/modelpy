# -*- coding: utf-8 -*-

import pandas as pd
import io
import os
import numpy as np
import datetime as dt
import sys
import pytz
import math


arguments = []

#parser.add_argument('--timestamp')
arguments.append(dict(arg='--path_forcing',\
                    help='directory of forcing data to initialize and constrain the ABL model simulations'))
arguments.append(dict(arg='--path_output',
                    help='output directory in which the output as subdirectories are stored'))#,default='/user/data/gent/gvo000/gvo00090/D2D/data/C4GL/')
arguments.append(dict(arg='--first_station_row',\
                    help='starting row number of stations table'))
arguments.append(dict(arg='--last_station_row',\
                    help='ending row number of stations table'))
arguments.append(dict(arg='--global_vars',\
                    help="global vars to update"))
arguments.append(dict(arg='--station_id',\
                    help="process a specific station id"))
arguments.append(dict(arg='--error_handling',\
                    default='dump_on_success',\
                    help="type of error handling: either\n - 'dump_on_success' (default)\n - 'dump_always'"))
arguments.append(dict(arg='--diag_tropo',\
                    default=['advt','advq','advu','advv'],\
                    help="field to diagnose the mean in the troposphere (<= 3000m)"))
arguments.append(dict(arg='--subset_forcing',
                    default='ini', 
                    help="This indicates which yaml subset to initialize with.  Most common options are 'ini' (default) and 'morning'."))
# Tuntime is usually specified from the afternoon profile. You can also just
# specify the simulation length in seconds

arguments.append(dict(arg='--split_by',\
                    type=int,
                    help="the maxmimum number of soundings that are contained in each output file of a station. -1 means unlimited (default). In case of arrays experiments, this is usually overwritten by 50."))

#arguments.append(dict(arg='--station-chunk',default=0)
arguments.append(dict(arg='--c4gl_path_lib',help="the path of the CLASS4GL program"))#,default='/user/data/gent/gvo000/gvo00090/D2D/software/CLASS/class4gl/lib')
arguments.append(dict(arg='--global_chunk_number',help="this is the batch number of the expected series of experiments according to split_by"))
arguments.append(dict(arg='--station_chunk_number',help="this is the batch number according to split_by in case of considering one station"))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    #parser.add_argument('--timestamp')
    for argument in arguments:
        name = argument.pop('arg')
        parser.add_argument(name,**argument)

    args = parser.parse_args()
else:
    class Namespace:
        def __init__(self,**kwargs):
            self.__dict__.update(kwargs)

    args = Namespace()
    for argument in arguments:
        if 'default' in argument.keys():
            args.__dict__[argument['arg'].strip('-')] = argument['default']
        else:
            args.__dict__[argument['arg'].strip('-')] = None
    print(args.__dict__)
        

def execute(**kwargs):
    # note that with args, we actually mean the same as those specified with
    # the argparse module above
    
    # overwrite the args according to the kwargs when the procedure is called
    # as module function
    for key,value in kwargs.items():
        args.__dict__[key]  = value
    
    print("-- begin arguments --")
    for key,value in args.__dict__.items():
         print(key,': ',value)
    print("-- end arguments ----")
    
    # load specified class4gl library
    if args.c4gl_path_lib is not None:
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
    
    if (args.global_vars is not None):
        globaldata = data_global()
        globaldata.load_datasets(recalc=0)
    
    
    # ========================
    print("getting a list of stations")
    # ========================
    
    # these are all the stations that are found in the input dataset
    all_stations = stations(args.path_forcing,suffix=args.subset_forcing,refetch_stations=False)
    
    # ====================================
    print('defining all_stations_select')
    # ====================================
    
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
    
    print(all_stations_select)
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
    
        if (args.split_by is None) or (args.split_by <= 0):
                raise ValueError("global_chunk_number is specified, but --split_by is not a strict positive number, so I don't know how to split the batch into chunks.")
    
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
            if args.split_by is not None:
                raise ValueError("Chunks are defined by --split_by, but I don't know which chunk to run. Please provide --global_chunk_number or --station_chunk_number, or leave out --split_by.")
            run_station_chunk = 0
            print("stations that are processed.",list(run_stations.index))
            
    
    #print(all_stations)
    print('Fetching initial/forcing records')
    records_morning = get_records(run_stations,\
                                  args.path_forcing,\
                                  subset=args.subset_forcing,
                                  refetch_records=False,
                                  )
    if len(records_morning) == 0:
        raise IOError("No initialization records records found in "+\
                      args.path_forcing+' (subset: '+args_subset_forcing+')')
    
    # note that if runtime is an integer number, we don't need to get the afternoon
    # profiles. 
    

    
    path_output = args.path_output
    
    os.system('mkdir -p '+path_output)
    for istation,current_station in run_stations.iterrows():
        # records_morning_station = records_morning.query('STNID == '+str(current_station.name))
        records_morning_station = records_morning.loc[(current_station.name):(current_station.name)]

        fn_morning = args.path_forcing+'/'+format(current_station.name,'05d')+'_'+args.subset_forcing+'.yaml'
        if os.path.isfile(fn_morning):
            file_morning = open(fn_morning)
        else:
            fn_morning = \
                 args.path_forcing+'/'+format(current_station.name,'05d')+\
                 '_'+str(run_station_chunk)+'_'+args.subset_forcing+'.yaml'
            file_morning = open(fn_morning)
    
        # if args.runtime == 'from_profile_pair':
        #     file_afternoon = open(args.path_forcing+'/'+format(current_station.name,'05d')+'_end.yaml')
        fn_ini = path_output+'/'+format(current_station.name,'05d')+'_'+\
                 str(int(run_station_chunk))+'_ini.yaml'
        file_ini = open(fn_ini,'w')
    
        #iexp = 0
        onerun = False
        print('starting station chunk number: '\
              +str(run_station_chunk)+' (chunk size:',args.split_by,')')
    
        skip_chunk = False
        if 'chunk' in records_morning.index.names:
            records_morning_station_chunk = records_morning_station.loc[(current_station.name,run_station_chunk):(current_station.name,run_station_chunk)]
        else:
            start_record = run_station_chunk*args.split_by if run_station_chunk is not 0 else 0
            end_record = (run_station_chunk+1)*args.split_by if args.split_by is not None else None
            if start_record >= (len(records_morning_station)):
                print("warning: outside of profile number range for station "+\
                      str(current_station)+". Skipping chunk number for this station.")
                skip_chunk = True
                records_morning_station_chunk = None
            else:
                records_morning_station_chunk = records_morning_station.iloc[start_record:end_record] #  [(int(args.split_by)*run_station_chunk):(int(args.split_by)*(run_station_chunk+1))]

        if not skip_chunk:
    
            isim = 0
            for (STNID,chunk,index),record_morning in records_morning_station_chunk.iterrows():
                    print('starting '+str(isim+1)+' out of '+\
                      str(len(records_morning_station_chunk) )+\
                      ' (station total: ',str(len(records_morning_station)),')')  
                
            
                    c4gli_morning = get_record_yaml(file_morning, 
                                                    record_morning.index_start, 
                                                    record_morning.index_end,
                                                    mode='model_input')
                    
                    
    
            
                    if args.global_vars is not None:
                        c4gli_morning.get_global_input(globaldata,only_keys=args.global_vars.strip().split(','))
    
                    onerun = True
    
                    print("dumping to "+str(file_ini)+ ' ('+fn_ini+')') 
                    c4gli_morning.dump(file_ini)
                        
                        
                    isim += 1
    
    
            file_ini.close()
            file_morning.close()
    
            if onerun:
                records_ini = get_records(pd.DataFrame([current_station]),\
                                                           path_output,\
                                                           getchunk = int(run_station_chunk),\
                                                           subset='ini',
                                                           refetch_records=True,
                                                           )
            else:
                # remove empty files
                os.system('rm '+fn_ini)
    
    # # align afternoon records with initial records, and set same index
    # records_afternoon.index = records_afternoon.ldatetime.dt.date
    # records_afternoon = records_afternoon.loc[records_ini.ldatetime.dt.date]
    # records_afternoon.index = records_ini.index
    
    # stations_for_iter = stations(path_output)
    # for STNID,station in stations_iterator(stations_for_iter):
    #     records_current_station_index = \
    #             (records_ini.index.get_level_values('STNID') == STNID)
    #     file_current_station_end_mod = STNID
    # 
    #     with \
    #     open(path_output+'/'+format(STNID,"05d")+'_ini.yaml','r') as file_station_ini, \
    #     open(path_output+'/'+format(STNID,"05d")+'_end_mod.yaml','r') as file_station_end_mod, \
    #     open(path_forcing+'/'+format(STNID,"05d")+'_afternoon.yaml','r') as file_station_afternoon:
    #         for (STNID,index),record_ini in records_iterator(records_ini):
    #             c4gli_ini = get_record_yaml(file_station_ini, 
    #                                         record_ini.index_start, 
    #                                         record_ini.index_end,
    #                                         mode='ini')
    #             #print('c4gli_in_ldatetime 3',c4gli_ini.pars.ldatetime)
    # 
    #             record_end_mod = records_end_mod.loc[(STNID,index)]
    #             c4gl_end_mod = get_record_yaml(file_station_end_mod, 
    #                                         record_end_mod.index_start, 
    #                                         record_end_mod.index_end,
    #                                         mode='mod')
    #             record_afternoon = records_afternoon.loc[(STNID,index)]
    #             c4gl_afternoon = get_record_yaml(file_station_afternoon, 
    #                                         record_afternoon.index_start, 
    #                                         record_afternoon.index_end,
    #                                         mode='ini')


if __name__ == '__main__':
    #execute(**vars(args))
    execute()
