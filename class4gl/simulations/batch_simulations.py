
import argparse

import pandas as pd
import os
import math
import numpy as np
import sys
import math
sys.path.insert(0, '/user/data/gent/gvo000/gvo00090/D2D/software/CLASS/class4gl/')
from class4gl import class4gl_input, data_global,class4gl
from interface_multi import stations,stations_iterator, records_iterator,get_record_yaml,get_records

odir = "/user/data/gent/gvo000/gvo00090/D2D/data/SOUNDINGS/GLOBAL/"
fn_stations = odir+'/igra-stations_sel.txt'
df_stations = pd.read_csv(fn_stations)

# if 'path-soundings' in args.__dict__.keys():
#     path_soundingsSET = args.__dict__['path-soundings']+'/'+SET+'/'
# else:



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--dataset')
    parser.add_argument('--path-soundings')
    #parser.add_argument('--exec',default='/user/data/gent/gvo000/gvo00090/D2D/scripts/C4GL/global_run.py')
    parser.add_argument('--exec')
    parser.add_argument('--experiments')#should be ';'-seperated list
    parser.add_argument('--split-by',default=-1)
    args = parser.parse_args()

experiments = args.experiments.split(';')
#SET = 'GLOBAL'
SET = args.dataset
print(args.experiments)

if 'path-soundings' in args.__dict__.keys():
    path_soundingsSET = args.__dict__['path-soundings']+'/'+SET+'/'
else:
    path_soundingsSET = '/user/data/gent/gvo000/gvo00090/D2D/data/SOUNDINGS/'+SET+'/'

all_stations = stations(path_soundingsSET,suffix='morning',refetch_stations=True).table
records_morning = get_records(all_stations,\
                              path_soundingsSET,\
                              subset='morning',
                              refetch_records=False,
                              )

for expname in experiments:
    #exp = EXP_DEFS[expname]
    path_exp = '/user/data/gent/gvo000/gvo00090/D2D/data/C4GL/'+SET+'_'+expname+'/'
    os.system('rm -R '+path_exp)

totalchunks = 0
for istation,current_station in all_stations.iterrows():
    records_morning_query = records_morning.query('STNID == '+str(current_station.name))
    chunks_current_station = math.ceil(float(len(records_morning_query))/float(args.split_by))
    totalchunks +=chunks_current_station

#if sys.argv[1] == 'qsub':
# with qsub
os.system('qsub /user/data/gent/gvo000/gvo00090/D2D/scripts/C4GL/global_run.pbs -t 0-'+str(totalchunks-1)+" -v dataset="+args.dataset+\
                                       ',split_by='+str(args.split_by)+\
                                       ',exec='+str(args.exec)+\
                                       ',experiments='+str(args.experiments))
# elif sys.argv[1] == 'wsub':
#     
#     # with wsub
#     STNlist = list(df_stations.iterrows())
#     NUMSTNS = len(STNlist)
#     PROCS = NUMSTNS 
#     BATCHSIZE = 1 #math.ceil(np.float(NUMSTNS)/np.float(PROCS))
# 
#     os.system('wsub -batch /user/data/gent/gvo000/gvo00090/D2D/scripts/C4GL/global_run.pbs -t 0-'+str(PROCS-1))

