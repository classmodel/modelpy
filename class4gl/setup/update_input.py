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
parser.add_argument('--path_input')#,default='/user/data/gent/gvo000/gvo00090/D2D/data/SOUNDINGS/')
parser.add_argument('--path_output')#,default='/user/data/gent/gvo000/gvo00090/D2D/data/C4GL/')
parser.add_argument('--first_station_row')
parser.add_argument('--last_station_row')
parser.add_argument('--updates')
parser.add_argument('--global_vars')
parser.add_argument('--station_id') # run a specific station id
parser.add_argument('--error_handling',default='dump_on_success')
parser.add_argument('--diag_tropo',default=None)#['advt','advq','advu','advv'])
parser.add_argument('--subset_input',default='morning') # this tells which yaml subset
                                                      # to initialize with.
                                                      # Most common options are
                                                      # 'morning' and 'ini'.
parser.add_argument('--subset_output',default='morning')


# Tuntime is usually specified from the afternoon profile. You can also just
# specify the simulation length in seconds

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



# iniitialize global data
globaldata = data_global()
if (args.updates is not None) and ('era_profiles' in args.updates.strip().split(",")):
    globaldata.sources = {**globaldata.sources,**{
            "ERAINT:t"     : "/user/data/gent/gvo000/gvo00090/EXT/data/ERA-INTERIM/by_var_nc/t_6hourly/t_*_6hourly.nc",
            "ERAINT:q"     : "/user/data/gent/gvo000/gvo00090/EXT/data/ERA-INTERIM/by_var_nc/q_6hourly/q_*_6hourly.nc",
            "ERAINT:u"     : "/user/data/gent/gvo000/gvo00090/EXT/data/ERA-INTERIM/by_var_nc/u_6hourly/u_*_6hourly.nc",
            "ERAINT:v"     : "/user/data/gent/gvo000/gvo00090/EXT/data/ERA-INTERIM/by_var_nc/v_6hourly/v_*_6hourly.nc",
            }}

# ...  and load initial data pages
globaldata.load_datasets(recalc=0)

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
  'ERA_NOAC':    {'sw_ac' : [],'sw_ap': True,'sw_lit': False},
  'NOAC':    {'sw_ac' : [],'sw_ap': True,'sw_lit': False},
  'ADV':{'sw_ac' : ['adv',],'sw_ap': True,'sw_lit': False},
  'W':  {'sw_ac' : ['w',],'sw_ap': True,'sw_lit': False},
  'AC': {'sw_ac' : ['adv','w'],'sw_ap': True,'sw_lit': False},
  'GLOBAL_NOAC':    {'sw_ac' : [],'sw_ap': True,'sw_lit': False},
  'GLOBAL_ADV':{'sw_ac' : ['adv',],'sw_ap': True,'sw_lit': False},
  'GLOBAL_ADV_ERA_NEW':{'sw_ac' : ['adv',],'sw_ap': True,'sw_lit': False},
    'GLOBAL_ADV_SHR':{'sw_ac' : ['adv',],'sw_ap': True,'sw_lit': False,'sw_shr':True},
  'GLOBAL_W':  {'sw_ac' : ['w',],'sw_ap': True,'sw_lit': False},
  'GLOBAL_AC': {'sw_ac' : ['adv','w'],'sw_ap': True,'sw_lit': False},
  'IOPS_NOAC':    {'sw_ac' : [],'sw_ap': True,'sw_lit': False},
  'IOPS_ADV':{'sw_ac' : ['adv',],'sw_ap': True,'sw_lit': False},
  'IOPS_W':  {'sw_ac' : ['w',],'sw_ap': True,'sw_lit': False},
  'IOPS_AC': {'sw_ac' : ['adv','w'],'sw_ap': True,'sw_lit': False},
}

# ========================
print("getting a list of stations")
# ========================

# these are all the stations that are found in the input dataset
all_stations = stations(args.path_input,suffix=args.subset_input,refetch_stations=False)

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
all_records_input_select = get_records(all_stations_select,\
                                         args.path_input,\
                                         subset=args.subset_input,
                                         refetch_records=False,
                                         )

# only run a specific chunck from the selection
if args.global_chunk_number is not None:
    if args.station_chunk_number is not None:
        raise ValueError('You need to specify either global-chunk-number or station-chunk-number, not both.')


    if not (int(args.split_by) > 0) :
            raise ValueError("global_chunk_number is specified, but --split_by is not a strict positive number, so I don't know how to split the batch into chunks.")

    run_station_chunk = None
    print('determining the station and its chunk number according global_chunk_number ('+args.global_chunk_number+')')
    totalchunks = 0
    stations_iter = all_stations_select.iterrows()
    in_current_chunk = False
    try:
        while not in_current_chunk:
            istation,current_station = stations_iter.__next__()
            all_records_input_station_select = all_records_input_select.query('STNID == '+str(current_station.name))
            chunks_current_station = math.ceil(float(len(all_records_input_station_select))/float(args.split_by))
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
records_input = get_records(run_stations,\
                              args.path_input,\
                              subset=args.subset_input,
                              refetch_records=False,
                              )



os.system('mkdir -p '+args.path_output)
for istation,current_station in run_stations.iterrows():
    print(istation,current_station)
    records_input_station = records_input.query('STNID == '+str(current_station.name))
    if (int(args.split_by) * int(run_station_chunk)) >= (len(records_input_station)):
        print("warning: outside of profile number range for station "+\
              str(current_station)+". Skipping chunk number for this station.")
    else:
        fn_input = args.path_input+'/'+format(current_station.name,'05d')+'_'+args.subset_input+'.yaml'
        if os.path.isfile(fn_input):
            file_input = open(fn_input)
        else:
            fn_input = \
                 args.path_input+'/'+format(current_station.name,'05d')+\
                 '_'+str(run_station_chunk)+'_'+args.subset_input+'.yaml'
            file_input = open(fn_input)

        fn_output = args.path_output+'/'+format(current_station.name,'05d')+'_'+\
                 str(int(run_station_chunk))+'_'+args.subset_output+'.yaml'
        file_output = open(fn_output,'w')

        #iexp = 0
        onerun = False
        print('starting station chunk number: '\
              +str(run_station_chunk)+'(size: '+str(args.split_by)+' soundings)')

        records_input_station_chunk = records_input_station.iloc[((run_station_chunk)*int(args.split_by)):((run_station_chunk+1)*int(args.split_by))] #  [(int(args.split_by)*run_station_chunk):(int(args.split_by)*(run_station_chunk+1))]

        isim = 0
        for (STNID,chunk,index),record_input in records_input_station_chunk.iterrows():
                print('starting '+str(isim+1)+' out of '+\
                  str(len(records_input_station_chunk) )+\
                  ' (station total: ',str(len(records_input_station)),')')  
            
        
                c4gli_output = get_record_yaml(file_input, 
                                                record_input.index_start, 
                                                record_input.index_end,
                                                mode='ini')
                if args.global_vars is not None:
                    c4gli_output.get_global_input(globaldata,only_keys=args.global_vars.strip().split(','))

                if args.diag_tropo is not None:
                    print('add tropospheric parameters on advection and subsidence (for diagnosis)')
                    seltropo = (c4gli_output.air_ac.p > c4gli_output.air_ac.p.iloc[-1]+ 3000.*(- 1.2 * 9.81 ))
                    profile_tropo = c4gli_output.air_ac[seltropo]
                    for var in args.diag_tropo:#['t','q','u','v',]:
                        if var[:3] == 'adv':
                            mean_adv_tropo = np.mean(profile_tropo[var+'_x']+profile_tropo[var+'_y'] )
                            c4gli_output.update(source='era-interim',pars={var+'_tropo':mean_adv_tropo})
                        else:
                            print("warning: tropospheric variable "+var+" not recognized")


                if (args.updates is not None) and ('era_profiles' in args.updates.strip().split(",")):
                    c4gli_output.get_global_input(globaldata,only_keys=['t','u','v','q','sp'])

                    c4gli_output.update(source='era-interim',pars={'Ps' : c4gli_output.pars.sp})

                    cp         = 1005.                 # specific heat of dry air [J kg-1 K-1]
                    Rd         = 287.                  # gas constant for dry air [J kg-1 K-1]
                    Rv         = 461.5                 # gas constant for moist air [J kg-1 K-1]
                    R = (Rd*(1.-c4gli_output.air_ac.q) + Rv*c4gli_output.air_ac.q)
                    rho = c4gli_output.air_ac.p/R/c4gli_output.air_ac.t
                    dz = c4gli_output.air_ac.delpdgrav/rho
                    z = [dz.iloc[-1]/2.]
                    for idz in list(reversed(range(0,len(dz)-1,1))):
                        z.append(z[-1]+ (dz[idz+1]+dz[idz])/2.)
                    z = list(reversed(z))

                    theta = c4gli_output.air_ac.t * \
                               (c4gli_output.pars.sp/(c4gli_output.air_ac.p))**(R/cp)
                    thetav   = theta*(1. + 0.61 * c4gli_output.air_ac.q)

                    
                    c4gli_output.update(source='era-interim',air_ac=pd.DataFrame({'z':list(z),
                                                                           'theta':list(theta),
                                                                           'thetav':list(thetav),
                                                                          }))
                    air_ap_input = c4gli_output.air_ac[::-1].reset_index().drop('index',axis=1)
                    air_ap_mode = 'b'
                    air_ap_input_source = c4gli_output.query_source('air_ac:theta')


                    c4gli_output.mixed_layer_fit(air_ap=air_ap_input,
                                         source=air_ap_input_source,
                                         mode=air_ap_mode)


                onerun = True
                
                c4gli_output.dump(file_output)
                    
                    
        file_output.close()
        file_input.close()

        if onerun:
            records_output = get_records(pd.DataFrame([current_station]),\
                                                       args.path_output,\
                                                       getchunk = int(run_station_chunk),\
                                                       subset=args.subset_output,
                                                       refetch_records=True,
                                                       )
        else:
            # remove empty files
            os.system('rm '+fn_output)

# # align afternoon records with initial records, and set same index
# records_afternoon.index = records_afternoon.ldatetime.dt.date
# records_afternoon = records_afternoon.loc[records_output.ldatetime.dt.date]
# records_afternoon.index = records_output.index

# stations_for_iter = stations(path_exp)
# for STNID,station in stations_iterator(stations_for_iter):
#     records_current_station_index = \
#             (records_output.index.get_level_values('STNID') == STNID)
#     file_current_station_mod = STNID
# 
#     with \
#     open(path_exp+'/'+format(STNID,"05d")+'_output.yaml','r') as file_station_output, \
#     open(path_exp+'/'+format(STNID,"05d")+'_mod.yaml','r') as file_station_mod, \
#     open(path_input+'/'+format(STNID,"05d")+'_afternoon.yaml','r') as file_station_afternoon:
#         for (STNID,index),record_output in records_iterator(records_output):
#             c4gli_output = get_record_yaml(file_station_output, 
#                                         record_output.index_start, 
#                                         record_output.index_end,
#                                         mode='ini')
#             #print('c4gli_in_ldatetime 3',c4gli_output.pars.ldatetime)
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

