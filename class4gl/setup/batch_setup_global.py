

''' 
Purpose: 
    launch array job to get sounding and other global forcing data in class4gl input format"
Usage:
    python start_setup_global.py

Author:
    Hendrik Wouters <hendrik.wouters@ugent.be>

'''

import pandas as pd
import os
import math
import numpy as np
import sys

odir = "/user/data/gent/gvo000/gvo00090/D2D/data/SOUNDINGS/GLOBAL/"
fn_stations = odir+'/igra-stations_sel.txt'
df_stations = pd.read_csv(fn_stations)

# if sys.argv[1] == 'qsub':
# with qsub
STNlist = list(df_stations.iterrows())
NUMSTNS = len(STNlist)
PROCS = len(STNlist) 
print(PROCS)
BATCHSIZE = math.ceil(np.float(NUMSTNS)/np.float(PROCS))
os.system('qsub /user/data/gent/gvo000/gvo00090/D2D/scripts/SOUNDINGS/setup_global.pbs -t 0-'+str(PROCS-1))
# elif sys.argv[1] == 'wsub':
#     
#     
#     # with wsub
#     STNlist = list(df_stations.iterrows())
#     NUMSTNS = len(STNlist)
#     PROCS = NUMSTNS 
#     BATCHSIZE = 1 #math.ceil(np.float(NUMSTNS)/np.float(PROCS))
# 
#     os.system('wsub -batch /user/data/gent/gvo000/gvo00090/D2D/scripts/SOUNDINGS/setup_global.pbs -t 0-'+str(PROCS-1))

