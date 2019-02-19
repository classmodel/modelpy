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

from time import sleep

parser = argparse.ArgumentParser()
#if __name__ == '__main__':
parser.add_argument('--exec') # chunk simulation script
parser.add_argument('--first_station_row',
                    help='starting row number of stations table')
parser.add_argument('--last_station_row')
parser.add_argument('--pbs_string',default=' -l walltime=2:0:0')
parser.add_argument('--station_id',
                    help="process a specific station id")
parser.add_argument('--error_handling')
parser.add_argument('--multi_processing_mode',default='pythonpool')
parser.add_argument('--cpu_count',type=int,default=2)
parser.add_argument('--subset_forcing',default='ini') 
                                        # this tells which yaml subset
                                        # to initialize with.
                                        # Most common options are
                                        # 'morning' and 'ini'.

# Tuntime is usually specified from the afternoon profile. You can also just
# specify the simulation length in seconds
parser.add_argument('--runtime',
                    help="set the runtime of the simulation in seconds, or get it from the daytime difference in the profile pairs 'from_profile_pair' (default)")
# delete folders of experiments before running them
parser.add_argument('--cleanup_output_directories',
                    default="False",
                    help="clean up output directories before executing the experiments")
parser.add_argument('--experiments', 
                    help="IDs of experiments, as a space-seperated list (default: 'BASE')")
parser.add_argument('--experiments_names', 
                    help="Alternative output names that are given to the experiments. By default, these are the same as --experiments") 
parser.add_argument('--split_by',
                    default=50,
                    type=int,
                    help="the maxmimum number of soundings that are contained in each output file of a station. -1 means unlimited. The default for array experiments is 50.")

parser.add_argument('--c4gl_path_lib',help="the path of the CLASS4GL program.")#,default='/user/data/gent/gvo000/gvo00090/D2D/software/CLASS/class4gl/lib')
parser.add_argument('--path_forcing',
                    help='directory of forcing data to initialize and constrain the ABL model simulations'
                   )
parser.add_argument('--path_experiments',
                    help='output directory in which the experiments as subdirectories are stored')#,default='/user/data/gent/gvo000/gvo00090/D2D/data/C4GL/')



#arguments only used for update_yaml.py
#parser.add_argument('--path_dataset') 
#parser.add_argument('--global_keys') 
batch_args = parser.parse_args()

if batch_args.c4gl_path_lib is not None:
    sys.path.insert(0, batch_args.c4gl_path_lib)
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



# #SET = 'GLOBAL'
# SET = batch_args.dataset

# path_forcingSET = batch_args.path_forcing+'/'+SET+'/'

print("getting all stations from "+batch_args.path_forcing)
# these are all the stations that are found in the input dataset
all_stations = stations(batch_args.path_forcing,suffix=batch_args.subset_forcing,refetch_stations=False)

print('defining all_stations_select')
# these are all the stations that are supposed to run by the whole batch (all
# chunks). We narrow it down according to the station(s) specified.
if batch_args.station_id is not None:
    print("Selecting stations by --station_id")
    stations_iter = stations_iterator(all_stations)
    STNID,run_station = stations_iter.set_STNID(STNID=int(batch_args.station_id))
    all_stations_select = pd.DataFrame([run_station])
else:
    print("Selecting stations from a row range in the table [--first_station_row,--last_station_row]")
    all_stations_select = pd.DataFrame(all_stations.table)
    if batch_args.last_station_row is not None:
        all_stations_select = all_station_select.iloc[:(int(batch_args.last_station)+1)]
    if batch_args.first_station_row is not None:
        all_stations_select = all_station_select.iloc[int(batch_args.first_station):]
print("station numbers included in the whole batch "+\
      "(all chunks):",list(all_stations_select.index))

print("getting all records of the whole batch")
all_records_morning_select = get_records(all_stations_select,\
                                         batch_args.path_forcing,\
                                         subset=batch_args.subset_forcing,\
                                         refetch_records=False,\
                                        )

print('splitting batch in --split_by='+str(batch_args.split_by)+' jobs.')
totalchunks = 0
for istation,current_station in all_stations_select.iterrows():
    records_morning_station_select = all_records_morning_select.query('STNID == '+str(current_station.name))
    chunks_current_station = math.ceil(float(len(records_morning_station_select))/float(batch_args.split_by))
    totalchunks +=chunks_current_station

print('total chunks of simulations (= size of array-job) per experiment: ' + str(totalchunks))

experiments = batch_args.experiments.strip(' ').split(' ')
if batch_args.experiments_names is not None:
    experiments_names = batch_args.experiments_names.strip(' ').split(' ')
    if len(experiments_names) != len(experiments):
        raise ValueError('Lenght of --experiments_names is different from --experiments')
else:
    experiments_names = experiments

odir_exists = False

cleanup = (batch_args.cleanup_output_directories == 'True')

if not cleanup:
    for expname in experiments_names:
        if os.path.exists(batch_args.path_experiments+'/'+expname):
            print("Output directory already exists: "+batch_args.path_experiments+'/'+expname+". ")
            odir_exists = True
if odir_exists:
    raise IOError("At least one of the output directories exists. Please use '--cleanup_output_directories True' to delete any output directory.")
else:
    for iexp,expname in enumerate(experiments_names):
        if cleanup:
            if os.path.exists(batch_args.path_experiments+'/'+expname):
                print("Warning! Output directory '"+batch_args.path_experiments+'/'+expname+"' exists! I'm removing it in 10 seconds!' Press ctrl-c to abort.")
                sleep(10)
                os.system("rm -R "+batch_args.path_experiments+'/'+expname)
        if batch_args.multi_processing_mode == 'qsub':
    
            # C4GLJOB_timestamp="+dt.datetime.now().isoformat()+",
            command = 'qsub '+batch_args.pbs_string+' '+batch_args.c4gl_path_lib+'/simulations/batch_simulations.pbs -t 0-'+\
                        str(totalchunks-1)+" -v C4GLJOB_experiments="+str(experiments[iexp])+",C4GLJOB_experiments_names="+str(expname)
            # propagate arguments towards the job script
            for argkey in batch_args.__dict__.keys():
                if ((argkey not in ['multi_processing_mode','cpu_count','experiments','experiments_names','pbs_string','cleanup_output_directories']) and \
                    # default values are specified in the simulation script, so
                    # excluded here
                    (batch_args.__dict__[argkey] is not None)
                   ):
                        command +=',C4GLJOB_'+argkey+'='+str(batch_args.__dict__[argkey])
    
            print('Submitting array job for experiment '+expname+': '+command)
            os.system(command)

        elif batch_args.multi_processing_mode == 'pythonpool':
            from multiprocessing import Pool                                       
            
            # # load moodule from absolute path
            # https://stackoverflow.com/questions/67631/how-to-import-a-module-given-the-full-path
            import importlib.util
            print(batch_args.exec)
            spec = importlib.util.spec_from_file_location("module.name", batch_args.exec)
            task_module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(task_module)
            print('hello')
            print(batch_args.exec)
            

            args_dict_current = {**batch_args.__dict__}

            # we avoid to pass Nones, so that defaults are taken from the child
            # script
            for removekey,nonevalue in batch_args.__dict__.items():
                if nonevalue is None:
                    args_dict_current.pop(removekey)

            # remove keys that are not relevant in the child script, so not
            # passed (or those that are redefined in the host script manually)
            for key in ['exec','multi_processing_mode','cpu_count','experiments','experiments_names','pbs_string','cleanup_output_directories']:
                if key in args_dict_current:
                    args_dict_current.pop(key)

            args_dict_current['experiments'] = experiments[iexp]
            args_dict_current['experiments_names'] = expname

            print(args_dict_current)
            all_tasks = []
            for ichunk in range(totalchunks):
                all_tasks.append({'global_chunk_number':str(ichunk),**args_dict_current}) 

            print(pd.DataFrame(all_tasks)) 
            def parallelize(analysis, filenames, processes):
                '''
                Call `analysis` for each file in the sequence `filenames`, using
                up to `processes` parallel processes. Wait for them all to complete
                and then return a list of results.
                '''
                return Pool(processes).map(analysis, filenames, chunksize = 1)
    
            def execute_kwargs(x):
                return task_module.execute(**x)

            parallelize(execute_kwargs,all_tasks,int(batch_args.cpu_count))

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

