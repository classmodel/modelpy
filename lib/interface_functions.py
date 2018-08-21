import pandas as pd
import numpy as np
import datetime as dt
import os
import xarray as xr
import sys
from contextlib import suppress
from time import sleep


sys.path.insert(0, '/user/data/gent/gvo000/gvo00090/D2D/software/CLASS/class4gl/')
from class4gl import class4gl_input, data_global,class4gl,units
from interface_functions import *
#from data_soundings import wyoming
import yaml
import glob
import pandas as pd
import json
import io
import subprocess
import pytz
from scipy.stats import mstats

from matplotlib.colors import LinearSegmentedColormap

class records_iterator(object):
    def __init__(self,records):
            
        self.records = records
        self.ix = -1 
        
    def __iter__(self):
        return self

    def __next__(self,jump=1):
        self.ix = (self.ix+jump) 
        if self.ix >= len(self.records.index):
            raise StopIteration

        return self.records.index[self.ix], self.records.iloc[self.ix]
    def __prev__(self):
        return self.__next__(self,jump=-1)


#'_afternoon.yaml'
def get_record_yaml(yaml_file,index_start,index_end,mode='mod'):
    filename = yaml_file.name
    #filename = path_yaml+'/'+format(current_station.name,'05d')+suffix
    #yaml_file = open(filename)

    #print('going to next observation',filename)
    yaml_file.seek(index_start)

    buf =  yaml_file.read(index_end- index_start).replace('inf','9e19').replace('nan','9e19').replace('---','')

    filebuffer = open(filename+'.buffer.yaml.'+str(index_start),'w')
    filebuffer.write(buf)
    filebuffer.close()
    # print("HHHEEELOOOO",filename+'.buffer.yaml'+str(index_start))
    
    command = '/apps/gent/CO7/sandybridge/software/Ruby/2.4.2-foss-2017b/bin/ruby -rjson -ryaml -e "'+"puts YAML.load_file('"+filename+".buffer.yaml."+str(index_start)+"').to_json"+'" > '+filename+'.buffer.json.'+str(index_start)+' '

    #command = '/apps/gent/CO7/sandybridge/software/Ruby/2.4.2-foss-2017b/bin/ruby -rjson -ryaml -e "'+"puts YAML.load(ARGF.read()).to_json"+'"'
    print(command)
    os.system(command)
    jsonstream = open(filename+'.buffer.json.'+str(index_start))
    record_dict = json.load(jsonstream)
    jsonstream.close()
    os.system('rm '+filename+'.buffer.yaml.'+str(index_start))


    if mode =='mod':
        modelout = class4gl()
        modelout.load_yaml_dict(record_dict)
        os.system('rm '+filename+'.buffer.json.'+str(index_start))

        return modelout
    elif mode == 'ini':

 
        # datetimes are incorrectly converted to strings. We need to convert them
        # again to datetimes
        for key,value in record_dict['pars'].items():
            # we don't want the key with columns that have none values
            if value is not None: 
                if key in ['lSunrise','lSunset','datetime','ldatetime','ldatetime_daylight','datetime_daylight',]:#(type(value) == str):
               # elif (type(value) == str):
                    record_dict['pars'][key] = dt.datetime.strptime(value,"%Y-%m-%d %H:%M:%S %z")

            if (value == 0.9e19) or (value == '.9e19'):
                record_dict['pars'][key] = np.nan
        for key in record_dict.keys():
            #print(key)
            if key in ['air_ap','air_balloon',]:
                #NNprint('check')
                for datakey,datavalue in record_dict[key].items():
                    record_dict[key][datakey] = [ np.nan if (x =='.9e19') else x for x in record_dict[key][datakey]]

        #os.system('rm '+filename+'.buffer.json.'+str(index_start))

        c4gli = class4gl_input()
        print(c4gli.logger,'hello')
        c4gli.load_yaml_dict(record_dict)
        os.system('rm '+filename+'.buffer.json.'+str(index_start))
        return c4gli






        # self.frames['stats']['records_current_station_index'] = \
        #     (self.frames['stats']['records_all_stations_index'].get_level_values('STNID')\
        #      == \
        #      self.frames['stats']['current_station'].name)

        # # create the value table of the records of the current station
        # tab_suffixes = \
        #         ['_mod','_obs','_obs_afternoon','_mod_stats','_obs_afternoon_stats','_ini_pct']
        # for tab_suffix in tab_suffixes:
        #     self.frames['stats']['records_current_station'+tab_suffix] = \
        #         self.frames['stats']['records_all_stations'+tab_suffix].iloc[self.frames['stats']['records_current_station_index']]


# class records_selection(object):
#     def __init__

# class records(object):
#     def __init__(self,stations,path_obs,path_mod):
#         self.stations = stations
#         self.path_obs = path_obs
#         self.path_mod = path_mod
# 
#         self.ini =       self.get_records(self.path_mod,'ini')
#         self.mod =       self.get_records(self.path_mod,'mod')
#         #self.morning =   self.get_records(self.path_obs,'morning')
#         self.afternoon = self.get_records(self.path_obs,'afternoon')
# 
#         
#         self.afternoon.index = self.afternoon.ldatetime.dt.date
#         self.afternoon = self.afternoon.loc[records_ini.ldatetime.dt.date]
# 
#         self.index = self.ini.index
#         self.mod.index = self.index
#         self.afternoon.index = self.index
# 
# 
#         #self.records_iterator = records_current_station_mod.iterrows()




class stations(object):
    def __init__(self,path,suffix='ini',refetch_stations=False):

        self.path = path

        self.file = self.path+'/stations_list.csv'
        if (os.path.isfile(self.file)) and (not refetch_stations):
            self.table = pd.read_csv(self.file)
        else:
            self.table = self.get_stations(suffix=suffix)
            self.table.to_csv(self.file)
        
        self.table = self.table.set_index('STNID')
        #print(self.table)

    def get_stations(self,suffix):
        stations_list_files = glob.glob(self.path+'/?????_0_'+suffix+'.yaml')
        if len(stations_list_files) == 0:
            stations_list_files = glob.glob(self.path+'/?????_'+suffix+'.yaml')
        stations_list_files.sort()
        print(stations_list_files)
        if len(stations_list_files) == 0:
            raise ValueError('no stations found that match "'+self.path+'/?????[_0]_'+suffix+'.yaml'+'"')
        stations_list = []
        for stations_list_file in stations_list_files:
            thisfile = open(stations_list_file,'r')
            yamlgen = yaml.load_all(thisfile)
            try:
                first_record  = yamlgen.__next__()
            except:
                first_record = None
            if first_record is not None:
                stations_list.append({})
                for column in ['STNID','latitude','longitude']:
                    #print(first_record['pars'].keys())
                    stations_list[-1][column] = first_record['pars'][column]
                stations_list[-1]['filename'] = os.path.split(stations_list_file)[1]
            yamlgen.close()
            thisfile.close()
    
        print(stations_list)
        return pd.DataFrame(stations_list)

class stations_iterator(object):
    def __init__(self,stations):
        self.stations = stations
        self.ix = -1 
    def __iter__(self):
        return self
    def __next__(self,jump=1):
        self.ix = (self.ix+jump) 
        if ((self.ix >= len(self.stations.table.index)) or (self.ix < 0 )):
            raise StopIteration
        self.ix = np.mod(self.ix,len(self.stations.table)) 
        return self.stations.table.index[self.ix], self.stations.table.iloc[self.ix]
    def set_row(self,row):
        self.ix = row
        return self.stations.table.index[self.ix], self.stations.table.iloc[self.ix]
    def set_STNID(self,STNID):
        self.ix = np.where((self.stations.table.index == STNID))[0][0]
        print(self.ix)
        print( self.stations.table.index[self.ix], self.stations.table.iloc[self.ix])
        return self.stations.table.index[self.ix], self.stations.table.iloc[self.ix]

    def __prev__(self):
        return self.__next__(self,jump=-1)
    def close():
        del(self.ix)

class records_iterator(object):
    def __init__(self,records):
            
        self.records = records
        self.ix = -1 
        
    def __iter__(self):
        return self

    def __next__(self,jump=1):
        self.ix = (self.ix+jump) 
        if self.ix >= len(self.records.index):
            raise StopIteration
        self.ix = np.mod(self.ix,len(self.records))
        return self.records.index[self.ix], self.records.iloc[self.ix]
    def __prev__(self):
        return self.__next__(self,jump=-1)


# #'_afternoon.yaml'
# def get_record_yaml(yaml_file,index_start,index_end):
#     filename = yaml_file.name
#     #filename = path_yaml+'/'+format(current_station.name,'05d')+suffix
#     #yaml_file = open(filename)
# 
#     #print('going to next observation',filename)
#     yaml_file.seek(index_start)
# 
#     buf =  yaml_file.read(index_end- index_start).replace('inf','9e19').replace('nan','9e19').replace('---','')
# 
#     filebuffer = open(filename+'.buffer.yaml.'+str(index_start),'w')
#     filebuffer.write(buf)
#     filebuffer.close()
#     # print("HHHEEELOOOO",filename+'.buffer.yaml'+str(index_start))
#     
#     command = '/apps/gent/CO7/sandybridge/software/Ruby/2.4.2-foss-2017b/bin/ruby -rjson -ryaml -e "'+"puts YAML.load_file('"+filename+".buffer.yaml."+str(index_start)+"').to_json"+'" > '+filename+'.buffer.json.'+str(index_start)+' '
# 
#     #command = '/apps/gent/CO7/sandybridge/software/Ruby/2.4.2-foss-2017b/bin/ruby -rjson -ryaml -e "'+"puts YAML.load(ARGF.read()).to_json"+'"'
#     print(command)
#     os.system(command)
#     jsonstream = open(filename+'.buffer.json.'+str(index_start))
#     record_dict = json.load(jsonstream)
#     jsonstream.close()
#     os.system('rm '+filename+'.buffer.yaml.'+str(index_start))
#  
#     # datetimes are incorrectly converted to strings. We need to convert them
#     # again to datetimes
#     for key,value in record_dict['pars'].items():
#         # we don't want the key with columns that have none values
#         if value is not None: 
#             if key in ['lSunrise','lSunset','datetime','ldatetime','ldatetime_daylight','ldatetime_daylight','datetime_daylight','datetime_daylight']:#(type(value) == str):
#            # elif (type(value) == str):
#                 record_dict['pars'][key] = dt.datetime.strptime(value,"%Y-%m-%d %H:%M:%S %z")
#                 
#                 # Workaround. Unfortunately, Ruby puts it in local time of the computer. Turn it back to UTC (note that UTC means actually local time)!!!
#                 record_dict['pars'][key] = record_dict['pars'][key].astimezone(pytz.UTC)
# 
#         if (value == 0.9e19) or (value == '.9e19'):
#             record_dict['pars'][key] = np.nan
#     for key in record_dict.keys():
#         print(key)
#         if key in ['air_ap','air_balloon',]:
#             print('check')
#             for datakey,datavalue in record_dict[key].items():
#                 record_dict[key][datakey] = [ np.nan if (x =='.9e19') else x for x in record_dict[key][datakey]]
# 
#     #os.system('rm '+filename+'.buffer.json.'+str(index_start))
# 
#     c4gli = class4gl_input()
#     c4gli.load_yaml_dict(record_dict)
#     return c4gli






        # self.frames['stats']['records_current_station_index'] = \
        #     (self.frames['stats']['records_all_stations_index'].get_level_values('STNID')\
        #      == \
        #      self.frames['stats']['current_station'].name)

        # # create the value table of the records of the current station
        # tab_suffixes = \
        #         ['_mod','_obs','_obs_afternoon','_mod_stats','_obs_afternoon_stats','_ini_pct']
        # for tab_suffix in tab_suffixes:
        #     self.frames['stats']['records_current_station'+tab_suffix] = \
        #         self.frames['stats']['records_all_stations'+tab_suffix].iloc[self.frames['stats']['records_current_station_index']]


# class records_selection(object):
#     def __init__

# class records(object):
#     def __init__(self,stations,path_obs,path_mod):
#         self.stations = stations
#         self.path_obs = path_obs
#         self.path_mod = path_mod
# 
#         self.ini =       self.get_records(self.path_mod,'ini')
#         self.mod =       self.get_records(self.path_mod,'mod')
#         #self.morning =   self.get_records(self.path_obs,'morning')
#         self.afternoon = self.get_records(self.path_obs,'afternoon')
# 
#         
#         self.afternoon.index = self.afternoon.ldatetime.dt.date
#         self.afternoon = self.afternoon.loc[records_ini.ldatetime.dt.date]
# 
#         self.index = self.ini.index
#         self.mod.index = self.index
#         self.afternoon.index = self.index
# 
# 
#         #self.records_iterator = records_current_station_mod.iterrows()



def get_records(stations,path_yaml,getchunk='all',subset='morning',refetch_records=False):

    records = pd.DataFrame()
    for STNID,station in stations.iterrows():
        dictfnchunks = []
        if getchunk is 'all':

            # we try the old single-chunk filename format first (usually for
            # original profile pairs)
            fn = path_yaml+'/'+format(STNID,'05d')+'_'+subset+'.yaml'
            if os.path.isfile(fn):
                chunk = 0
                dictfnchunks.append(dict(fn=fn,chunk=chunk))

            # otherwise, we use the new multi-chunk filename format
            else:
                chunk = 0
                end_of_chunks = False
                while not end_of_chunks:
                    fn = path_yaml+'/'+format(STNID,'05d')+'_'+str(chunk)+'_'+subset+'.yaml'
                    if os.path.isfile(fn):
                        dictfnchunks.append(dict(fn=fn,chunk=chunk))
                    else:
                        end_of_chunks = True
                    chunk += 1

            # globyamlfilenames = path_yaml+'/'+format(STNID,'05d')+'*_'+subset+'.yaml'
            # yamlfilenames = glob.glob(globyamlfilenames)
            # yamlfilenames.sort()
        else:
            fn = path_yaml+'/'+format(STNID,'05d')+'_'+str(getchunk)+'_'+subset+'.yaml'
            dictfnchunks.append(dict(fn=fn,chunk=getchunk))
            
        if len(dictfnchunks) > 0:
            for dictfnchunk in dictfnchunks:
                yamlfilename = dictfnchunk['fn']
                chunk = dictfnchunk['chunk']
                print(chunk)

                #pklfilename = path_yaml+'/'+format(STNID,'05d')+'_'+subset+'.pkl'
                pklfilename = yamlfilename.replace('.yaml','.pkl')

                #print(yamlfilename+": "+str(os.path.getmtime(yamlfilename)))
                #print(pklfilename+": "+str(os.path.getmtime(pklfilename)))
                generate_pkl = False
                if not os.path.isfile(pklfilename): 
                    print('pkl file does not exist. I generate "'+\
                          pklfilename+'" from "'+yamlfilename+'"...')
                    generate_pkl = True
                elif not (os.path.getmtime(yamlfilename) <  \
                    os.path.getmtime(pklfilename)):
                    print('pkl file older than yaml file, so I regenerate "'+\
                          pklfilename+'" from "'+yamlfilename+'"...')
                    generate_pkl = True

                if refetch_records:
                    print('refetch_records flag is True. I regenerate "'+\
                          pklfilename+'" from "'+yamlfilename+'"...')
                    generate_pkl = True
                if not generate_pkl:
                    records = pd.concat([records,pd.read_pickle(pklfilename)])
                   # irecord = 0
                else:
                    with open(yamlfilename) as yaml_file:

                        dictout = {}

                        next_record_found = False
                        end_of_file = False
                        while (not next_record_found) and (not end_of_file):
                            linebuffer = yaml_file.readline()
                            next_record_found = (linebuffer == '---\n')
                            end_of_file = (linebuffer == '')
                        next_tell = yaml_file.tell()
                        
                        while not end_of_file:

                            print(' next record:',next_tell)
                            current_tell = next_tell
                            next_record_found = False
                            yaml_file.seek(current_tell)
                            filebuffer = open(yamlfilename+'.buffer.yaml.'+str(current_tell),'w')
                            linebuffer = ''
                            while ( (not next_record_found) and (not end_of_file)):
                                filebuffer.write(linebuffer.replace('inf','0').replace('nan','0'))
                                linebuffer = yaml_file.readline()
                                next_record_found = (linebuffer == '---\n')
                                end_of_file = (linebuffer == '')
                            filebuffer.close()
                            
                            next_tell = yaml_file.tell()
                            index_start = current_tell
                            index_end = next_tell

                            
                            #if ((irecord >= start) and (np.mod(irecord - start,2) == 0.) :
                            command = '/apps/gent/CO7/sandybridge/software/Ruby/2.4.2-foss-2017b/bin/ruby -rjson -ryaml -e "'+"puts YAML.load_file('"+yamlfilename+".buffer.yaml."+str(current_tell)+"').to_json"+'" > '+yamlfilename+'.buffer.json.'+str(current_tell)+' ' 
                            print(command)
                            
                            os.system(command)
                            #jsonoutput = subprocess.check_output(command,shell=True) 
                            #print(jsonoutput)
                            #jsonstream = io.StringIO(jsonoutput)
                            jsonstream = open(yamlfilename+'.buffer.json.'+str(current_tell))
                            record = json.load(jsonstream)
                            dictouttemp = {}
                            for key,value in record['pars'].items():
                                # we don't want the key with columns that have none values
                                if value is not None: 
                                   regular_numeric_types =[ type(x) for x in[0,False,0.0]]
                                   if (type(value) in regular_numeric_types):
                                        dictouttemp[key] = value
                                   elif key in ['lSunrise','lSunset','datetime','ldatetime','datetime_daylight','datetime_daylight','ldatetime_daylight','ldatetime_daylight']:#(type(value) == str):
                                       #print (key,value) # dictouttemp[key] = dt.datetime.strptime(value[:-6],"%Y-%m-%d %H:%M:%S")
                                       dictouttemp[key] = dt.datetime.strptime(value,"%Y-%m-%d %H:%M:%S %z")
                                       # Workaround. Unfortunately, Ruby puts it in local time of the computer. Turn it back to UTC (note that UTC means actually local time)!!!
                                       dictouttemp[key] = dictouttemp[key].astimezone(pytz.UTC)
                            recordindex = record['index']
                            dictouttemp['chunk'] = chunk
                            dictouttemp['index_start'] = index_start
                            dictouttemp['index_end'] = index_end
                            os.system('rm '+yamlfilename+'.buffer.json.'+str(current_tell))
                            for key,value in dictouttemp.items():
                                if key not in dictout.keys():
                                    dictout[key] = {}
                                dictout[key][(STNID,chunk,recordindex)] = dictouttemp[key]
                            print(' obs record registered')
                            jsonstream.close()
                            os.system('rm '+yamlfilename+'.buffer.yaml.'+str(current_tell))
                    records_station = pd.DataFrame.from_dict(dictout)
                    records_station.index.set_names(('STNID','chunk','index'),inplace=True)
                    print('writing table file ('+pklfilename+') for station '\
                          +str(STNID))
                    records_station.to_pickle(pklfilename)
                    # else:
                    #     os.system('rm '+pklfilename)
                    records = pd.concat([records,records_station])
    return records

def stdrel(mod,obs,columns):
    stdrel = pd.DataFrame(columns = columns)
    for column in columns:
        stdrel[column] = \
                (mod.groupby('STNID')[column].transform('mean') -
                 obs.groupby('STNID')[column].transform('mean')) /\
                obs.groupby('STNID')[column].transform('std') + \
                (mod[column] -
                 mod.groupby('STNID')[column].transform('mean')) /\
                obs.groupby('STNID')[column].transform('std') 
    return stdrel

def pct(obs,columns):
    pct = pd.DataFrame(columns=columns)
    for column in columns:
        #print(column)
        pct[column] = ""
        pct[column] = obs[column].rank(pct=True)
    return pct

def tendencies(mod_afternoon,obs_afternoon,obs_morning,keys):
    stats = pd.DataFrame()
    for key in keys: 
        stats['d'+key+'dt'] = ""
        stats['d'+key+'dt'] = (mod_afternoon[key] - obs_morning[key])/ \
                              (obs_afternoon.ldatetime - \
                               obs_morning.ldatetime).dt.seconds*3600.
    return stats
