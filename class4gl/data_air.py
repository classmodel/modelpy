import numpy as np

from bs4 import BeautifulSoup
import pandas as pd
import datetime as dt
#import pylab as pl
import io
import os
import calendar

import importlib
spam_loader = importlib.find_loader('Pysolar')
found = spam_loader is not None
if found:
    import Pysolar
else:
    import pysolar as Pysolar
import Pysolar.util as Pysolarutil





#from urllib import request
def blh(HAGL,THTV,WSPD,RiBc = 0.5,RiBce = 0.25):
    
    #initialize error BLH
    BLHe = 0.
    eps = 2.#security limit
    iTHTV_0 = np.where(~np.isnan(THTV))[0]
    if len(iTHTV_0) > 0:
        iTHTV_0 = iTHTV_0[0]
        THTV_0 = THTV[iTHTV_0]
    else:
        THTV_0 = np.nan
    RiB = 9.81/THTV_0 * ( THTV - THTV_0) * HAGL / WSPD**2.
    
    
    #RiB = 9.81/THTV_0 * ( THTV[i-1] +  (HGHT[i] - HGHT[i-1])/ - THTV_0) * HAGL / WSPD**2
    #RiB - RiBc = 0
    
    #best guess of BLH
    
    #print("RiB: ",RiB)
    #print("RiBc: ",RiBc)
    
    
    
    BLHi = np.where(RiB > RiBc)[0]
    if len(BLHi ) > 0:
        BLHi = BLHi[0]
        #print("BLHi: ",BLHi)
        BLH = (HAGL[BLHi] - HAGL[BLHi-1])/(RiB[BLHi] -RiB[BLHi-1]) * (RiBc - RiB[BLHi-1]) + HAGL[BLHi-1]
        
        # possible error is calculated as the difference height levels used for the interpolation
        BLHu = np.max([BLH,HAGL[BLHi]-eps])
        BLHd = np.min([BLH,HAGL[BLHi-1]+eps])
        # calculate an alternative BLH based on another critical Richardson number (RiBce):
        BLHi =np.where(RiB > RiBce)[0]
        if len(BLHi ) > 0:    
            BLHi = BLHi[0]
                
            BLHa = (HAGL[BLHi] - HAGL[BLHi-1])/(RiB[BLHi] -RiB[BLHi-1]) * (RiBc - RiB[BLHi-1]) + HAGL[BLHi-1]
            BLHu = np.max([BLHu,HAGL[BLHi]-eps])
            BLHd = np.min([BLHd,HAGL[BLHi-1]+eps])
            
            BLHu = np.max([BLHu,BLH + abs(BLH-BLHa)])
            BLHd = np.min([BLHd,BLH - abs(BLH-BLHa)])
        
        else:
            BLH,BLHu,BLHd = np.nan, np.nan,np.nan

    else:
        BLH,BLHu,BLHd = np.nan, np.nan,np.nan
        
    return BLH,BLHu,BLHd

def dtrange(STARTTIME,ENDTIME,TIMEJUMP=dt.timedelta(hours=24)):
    STEPS = int((ENDTIME - STARTTIME).total_seconds()/TIMEJUMP.total_seconds())
    return [STARTTIME + TIMEJUMP*i for i in range(0,STEPS)]


#from os import listdir
#from os.path import isfile #,join
import glob


class wyoming(object):
    def __init__(self):
       self.status = 'init'
       self.found = False
       self.DT = None
       self.current = None
       #self.mode = 'b'
       self.profile_type = 'wyoming'  
       self.MONTHS = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
       self.PATH = "/user/data/gent/gvo000/gvo00090/EXT/data/SOUNDINGS/"
         
    def set_STNM(self,STNM):
        self.__init__()
        self.STNM = STNM
        self.FILES = glob.glob(self.PATH+'/????/SOUNDINGS_????_'+format(STNM,'05d')+".html")
        self.FILES = [os.path.realpath(FILE) for FILE in self.FILES]
        self.current = None
        self.found = False
        self.FILES.sort()
        
    def find_first(self,year=None,get_atm=False):
        self.found = False    
                
        # check first file/year or specified year
        if year == None:
            self.iFN = 0
            self.FN = self.FILES[self.iFN]
        else:
            self.FN = os.path.realpath(self.PATH+'/'+str(year)+'/SOUNDINGS_'+str(year)+'_'+format(self.STNM,'05d')+".html")
            self.iFN = self.FILES.index(self.FN)
        self.sounding_series = BeautifulSoup(open(self.FN), "html.parser")
        self.current = self.sounding_series.find('h2')
        keepsearching = (self.current is None) #if we don't want later years, add here: "and (year is None)"
        
        # go through other files and find first sounding when year is not specified
        self.iFN=self.iFN+1
        while keepsearching:
            self.FN = self.FILES[self.iFN]
            self.sounding_series = BeautifulSoup(open(self.FN), "html.parser")
            self.current = self.sounding_series.find('h2')
            self.iFN=self.iFN+1
            keepsearching = (self.current is None) and (self.iFN < len(self.FILES))
        self.found = (self.current is not None)

        self.status = 'fetch'
        if self.found:
            self.DT = dt.datetime(int(self.current.text[-4:]),self.MONTHS.index(self.current.text[-8:-5])+1,int(self.current.text[-11:-9]),int(self.current.text[-15:-13]))
        
        if self.found and get_atm:
            self.get_values_air_input()
        
    
    def find(self,DT,get_atm=False):
        
        self.found = False
        keepsearching = True
        #print(DT)
        # we open a new file only when it's needed. Otherwise we just scroll to the right sounding.  
        if not ((self.current is not None) and (DT >= self.DT) and (self.DT.year == DT.year)):
            self.DT = DT
            self.FN = os.path.realpath(self.PATH+"/"+self.DT.strftime("%Y")+"/SOUNDINGS_"+self.DT.strftime("%Y")+"_"+format(self.STNM,'05d')+".html")
            self.iFN = self.FILES.index(self.FN)
            self.sounding_series = BeautifulSoup(open(self.FN), "html.parser")
            self.current = self.sounding_series.find('h2')
            
        keepsearching = (self.current is not None)
        while keepsearching:
            DTcurrent = dt.datetime(int(self.current.text[-4:]),self.MONTHS.index(self.current.text[-8:-5])+1,int(self.current.text[-11:-9]),int(self.current.text[-15:-13]))
            if DTcurrent == DT:
                self.found = True
                keepsearching = False
                if get_atm:
                    self.get_values_air_input()
                    self.DT = dt.datetime(int(self.current.text[-4:]),self.MONTHS.index(self.current.text[-8:-5])+1,int(self.current.text[-11:-9]),int(self.current.text[-15:-13]))
            elif DTcurrent > DT:
                keepsearching = False
                self.current = None
            else:
                self.current = self.current.find_next('h2')
                if self.current is None:
                    keepsearching = False
        self.found = (self.current is not None)
        self.status = 'fetch'

    def find_next(self,get_atm=False):
        self.found = False
        self.DT = None
        if self.current is None:
            self.find_first()
        else:                
            self.current = self.current.find_next('h2')
            self.found = (self.current is not None)
            keepsearching = ((self.current is None) and ((self.iFN+1) < len(self.FILES)))
            while keepsearching:
                self.iFN=self.iFN+1
                self.FN = self.FILES[self.iFN]
                self.sounding_series = BeautifulSoup(open(self.FN), "html.parser")
                self.current = self.sounding_series.find('h2')
                
                self.found = (self.current is not None)
                keepsearching = ((self.current is None) and (self.iFN < len(self.FILES)))
        if self.found:        
            self.DT = dt.datetime(int(self.current.text[-4:]),self.MONTHS.index(self.current.text[-8:-5])+1,int(self.current.text[-11:-9]),int(self.current.text[-15:-13]))
        if self.found and get_atm:
            self.get_values_air_input()
       


    def get_values_air_input(self,latitude=None,longitude=None):

        # for iDT,DT in enumerate(DTS):
        
            #websource = urllib.request.urlopen(webpage)
        #soup = BeautifulSoup(open(webpage), "html.parser")
        
       
        #workaround for ...last line has <pre> which results in stringlike first column
        string = self.current.find_next('pre').text
        string = string.split('\n')[:-1]
        string =  '\n'.join(string)
        columns = [ 'PRES', 'HGHT', 'TEMP', 'DWPT', 'RELH', 'MIXR', 'DRCT','SKNT' , 'THTA','THTE', 'THTV']             
        ONE_COLUMN = pd.read_fwf(io.StringIO(str(string)),widths=[7]*11,names=columns).iloc[5:-1]
        #ONE_COLUMN = pd.read_table(io.StringIO(str(string)),sep=r"\s*",skiprows=[0,1,3,4])
        
        #string =  soup.pre.next_sibling.next_sibling
        
        string = self.current.find_next('pre').find_next('pre').text

        PARAMS = pd.read_fwf(io.StringIO(str(string)),widths=[43,1,20],names=['descr','dummy','value']).iloc[1:-1].drop("dummy",1).set_index("descr").T
        #PARAMS.insert(0,'date',DT)

        PARAMS.insert(0,'datetime', dt.datetime.strptime(str(PARAMS['Observation time'][0]),"%y%m%d/%H%M"))
        PARAMS.insert(0,'STNID',int(PARAMS['Station number']))
        
        THTV = np.array(ONE_COLUMN.THTV,dtype='float')
        #THTA = np.array(ONE_COLUMN.THTA,dtype='float')
        HGHT = np.array(ONE_COLUMN.HGHT,dtype='float')
        HAGL = HGHT - np.float(PARAMS['Station elevation'])
        ONE_COLUMN.insert(0,'HAGL',HAGL)

        
        
        
        MIXR = np.array(ONE_COLUMN.MIXR,dtype='float')
        QABS = (MIXR/1000.)/(MIXR/1000.+1.)
        ONE_COLUMN.insert(0,'QABS',QABS)
        
        WSPD =0.51444 * np.array(ONE_COLUMN.SKNT,dtype='float')

        #mixed layer potential temperature
        #THTVM = np.float(VALUE['Mean mixed layer potential temperature'].iloc[0])

        #THTV_0 = THTA[np.where(~np.isnan(THTA))[0][0]]
        #RiBV = 9.81/THTV_0 * ( THTV - THTV_0) * HGHT / WSPD**2
        #RiBA = 9.81/THTA_0 * ( THTA - THTA_0) * HGHT / WSPD**2

        BLHV,BLHVu,BLHVd = blh(HAGL,THTV,WSPD)
        BLHV = np.max((BLHV,10.))
        BLHVu = np.max((BLHVu,10.))
        BLHVd = np.max((BLHVd,10.))
        #BLHA,BLHAu,BLHAd = BLH(HGHT,THTA,WSPD)

        #security values for mixed-layer jump values dthetav, dtheta and dq
        
        # fit new profiles taking the above-estimated mixed-layer height
        ONE_COLUMNNEW = []
        for BLH in [BLHV,BLHVu,BLHVd]:
            ONE_COLUMNNEW.append(pd.DataFrame())
            
            HAGLNEW = np.array([2.,BLH,BLH]+list(HAGL[HAGL > BLH]),dtype=np.float)
            ONE_COLUMNNEW[-1].insert(0,'HAGL',HAGLNEW)
            
            listHAGLNEW = list(HAGLNEW)
            for icol,col in enumerate(['THTA','THTV','QABS','SKNT','DRCT','PRES']):
                
                # get index of lowest valid observation. This seems to vary
                idxvalid = np.where((np.array(HAGL) >= 0) & (~pd.isnull(np.array(ONE_COLUMN[col],dtype=np.float) )))[0]
                if len(idxvalid) > 0:
                    #print('idxvalid',idxvalid)
                    if len(np.where(HAGL[idxvalid[0]:] < BLH)[0]) >= 3:
                        meanabl = np.nanmean(np.array(ONE_COLUMN[col][HAGL < BLH][(idxvalid[0]+1):],dtype=np.float))
                    else:
                        meanabl = np.nanmean(ONE_COLUMN[col][idxvalid[0]:(idxvalid[0]+1)],dtype=np.float)                    
                else:
                    meanabl = np.nanmean(ONE_COLUMN[col][0:1],dtype=np.float)
                    #print(col,meanabl)
               
                
                # if col == 'PRES':
                #     meanabl =  
            
                new_pro_h = list(np.array(ONE_COLUMN[col][HAGL > BLH],dtype=np.float))
                #THTVM = np.nanmean(THTV[HAGL <= BLH])
                #print("new_pro_h",new_pro_h)
                # calculate jump ath the top of the mixed layer
                if col in ['THTA','THTV',]:
                    #for moisture
                    #print('hello:',(new_pro_h[1] - new_pro_h[0])/(listHAGLNEW[4] - listHAGLNEW[3])*(BLH-listHAGLNEW[3]))
                    #print('hello:',new_pro_h[1] , new_pro_h[0],listHAGLNEW[4] , listHAGLNEW[3],BLH,listHAGLNEW[3])
                    if len(listHAGLNEW) > 4:
                        #print(type(new_pro_h[1]),type(new_pro_h[0]),type(listHAGLNEW[4]),type(listHAGLNEW[3]),type(BLH),type(meanabl))
                        dtheta_pre = (new_pro_h[1] - new_pro_h[0])/(listHAGLNEW[4] - listHAGLNEW[3])*(BLH-listHAGLNEW[3]) + new_pro_h[0] - meanabl 
                        dtheta = np.max((0.1,dtheta_pre))
                        #meanabl = meanabl - (dtheta - dtheta_pre)
                        #print('dtheta_pre',dtheta_pre)
                        #print('dtheta',dtheta)
                        #print('meanabl',meanabl)
                        #stop
                        
                    else:
                        dtheta = np.nan
                else:
                    if len(listHAGLNEW) > 4:
                        #for moisture (it can have both negative and positive slope)
                        dtheta = ((new_pro_h[1] - new_pro_h[0])/(listHAGLNEW[4] - listHAGLNEW[3])*(BLH-listHAGLNEW[3]) + new_pro_h[0] - meanabl ) 
                    else:
                        dtheta = np.nan
                #print('dtheta',dtheta)
                
                new_pro = np.array([meanabl,meanabl,meanabl+dtheta]+new_pro_h,dtype=np.float)
            
                
                ONE_COLUMNNEW[-1].insert(len(ONE_COLUMNNEW[-1].columns),col,new_pro)
                
            #QABSM = np.nanmean(QABS[HAGL <= BLH])
            #QABSNEW = np.array([QABSM,QABSM]+list(QABS[HAGL > BLH]))
            #ONE_COLUMNNEW.append(pd.DataFrame(zip(HAGLNEW,THTVNEW,QABSNEW),columns=('HAGL','THTV','QABS')))
            
        # we just make a copy of the fields, so that it can be read correctly by CLASS 
        for dataonecolumn in ONE_COLUMNNEW+[ONE_COLUMN]:
            dataonecolumn.insert(len(dataonecolumn.columns),'p_pro',np.array(dataonecolumn.PRES,dtype=np.float)*100.)
            dataonecolumn.insert(len(dataonecolumn.columns),'z_pro',np.array(dataonecolumn.HAGL,dtype=np.float))
            dataonecolumn.insert(len(dataonecolumn.columns),'theta_pro',np.array(dataonecolumn.THTA,dtype=np.float))
            dataonecolumn.insert(len(dataonecolumn.columns),'thetav_pro',np.array(dataonecolumn.THTV,dtype=np.float))
            dataonecolumn.insert(len(dataonecolumn.columns),'q_pro',np.array(dataonecolumn.QABS,dtype=np.float))
            
            angle_x = (90.-np.array(dataonecolumn.DRCT,dtype=np.float))/180.*np.pi # assuming that wind in direction of the south is 0 degrees.
            spd = 0.51444* np.array(dataonecolumn.SKNT,dtype=np.float)
        
            dataonecolumn.insert(len(dataonecolumn.columns),'u_pro',spd * np.sin(angle_x))
            dataonecolumn.insert(len(dataonecolumn.columns),'v_pro',spd * np.cos(angle_x))


        # assign fields adopted by CLASS
        if self.mode == 'o': #original 
            PARAMS.insert(0,'h',   np.float(BLHV))
        elif self.mode == 'b':
            PARAMS.insert(0,'h',   np.float(BLHV))
        elif self.mode == 'u':
            PARAMS.insert(0,'h',   BLHVu)
        elif self.mode == 'd':
            PARAMS.insert(0,'h',   BLHVd)
        else:
            PARAMS.insert(0,'h',   BLHV)
            

        try:
            PARAMS.insert(0,'lat', np.float(PARAMS['Station latitude'][0]))
            PARAMS.insert(0,'latitude', np.float(PARAMS['Station latitude'][0]))
        except:
            print("could not convert latitude coordinate")
            PARAMS.insert(0,'latitude', np.nan)
            PARAMS.insert(0,'lat', np.nan)
        try:
            PARAMS.insert(0,'longitude', np.float(PARAMS['Station longitude'][0]))
            # we set the actual input parameter value of lon to zero as we are working in local time (as if we were in Greenwhich) 
            PARAMS.insert(0,'lon', 0.)
        except:
            print("could not convert longitude coordinate")
            PARAMS.insert(0,'longitude', np.nan)
            PARAMS.insert(0,'lon', 0.)

        if latitude is not None:
            print('overwriting latitude with specified value')
            PARAMS['latitude'] = np.float(latitude)
            PARAMS['lat'] = np.float(latitude)
        if longitude is not None:
            print('overwriting longitude with specified value')
            PARAMS['longitude'] = np.float(longitude)
        try:
            #this is the local suntime datetime from which we calculate the hour of the day (assuming we would be in greenwhich hence taking lon=0)
            PARAMS['ldatetime'] = PARAMS.datetime.value + dt.timedelta(hours=PARAMS.longitude.value/360.*24.) 
            PARAMS['SolarAltitude'] = Pysolar.GetAltitude(PARAMS.lat.value,PARAMS.longitude.value,PARAMS.datetime.value)
            PARAMS['SolarAzimuth'] = Pysolar.GetAzimuth(PARAMS.lat.value,PARAMS.longitude.value,PARAMS.datetime.value)
            PARAMS['lSunrise'], PARAMS['lSunset'] = Pysolarutil.GetSunriseSunset(PARAMS.lat.value,0.,PARAMS.datetime.value,0.)
            # This is the nearest datetime when sun is up (for class)
            PARAMS['ldatetime_daylight'] = np.min(np.max(PARAMS['ldatetime'].value ,PARAMS['lSunrise'].value),PARAMS['lSunset'].value) 
            # apply the same time shift for UTC datetime
            PARAMS['datetime_daylight'] = PARAMS.datetime.value  + (PARAMS.ldatetime_daylight.value  - PARAMS.ldatetime.value)
            
        except:
            print("could not get local times for profile, perhaps because of wrong longitude or latitude in the profile description")
            PARAMS['ldatetime'] = dt.datetime(1900,1,1)
            PARAMS['SolarAltitude'] = np.nan #Pysolar.GetAltitude(PARAMS.lat.value,PARAMS.lon.value,PARAMS.datetime.value)
            PARAMS['SolarAzimuth'] = np.nan #Pysolar.GetAzimuth(PARAMS.lat.value,PARAMS.lon.value,PARAMS.datetime.value)
            PARAMS['lSunrise'], PARAMS['lSunset'] = dt.datetime(1900,1,1), dt.datetime(1900,1,1) #Pysolarutil.GetSunriseSunset(PARAMS.lat.value,0.,PARAMS.datetime.value,0.)
            PARAMS['ldatetime_daylight'] =PARAMS['ldatetime'].value
            PARAMS['datetime_daylight'] =PARAMS['datetime'].value

        

        PARAMS.insert(0,'day', PARAMS['ldatetime'][0].day)
        # as we are forcing lon equal to zero this is is expressed in local suntime
        PARAMS.insert(0,'tstart', PARAMS['ldatetime_daylight'][0].hour + PARAMS['ldatetime_daylight'][0].minute/60. + PARAMS['ldatetime_daylight'][0].second/3600.)

           
        ONE_COLUMNb = ONE_COLUMNNEW[0]
        ONE_COLUMNu = ONE_COLUMNNEW[1]
        ONE_COLUMNd = ONE_COLUMNNEW[2]
        

        THTVM = np.nanmean(THTV[HAGL <= BLHV])
        PARAMS.insert(len(PARAMS.columns),'THTVM',THTVM)
        
        QABSM = np.nanmean(QABS[HAGL <= BLHV])
        PARAMS.insert(len(PARAMS.columns),'QABSM',QABSM)
        
        PARAMS.insert(len(PARAMS.columns),'BLHV',BLHV)
        PARAMS.insert(len(PARAMS.columns),'BLHVu',BLHVu)
        PARAMS.insert(len(PARAMS.columns),'BLHVd',BLHVd)  

        BLHVe = abs(BLHV - BLHVu)
        BLHVe = max(BLHVe,abs(BLHV - BLHVd))

        #PARAMS.insert(0,'dq',0.)
        
        PARAMS.insert(len(PARAMS.columns),'BLHVe',BLHVe)  
        PARAMS.insert(0,'Ps',np.array(ONE_COLUMN.PRES,dtype='float')[0]*100.)
        #PARAMS.insert(len(PARAMS.columns),'STNM',STNM)
        #PARAMS.insert(len(PARAMS.columns),'PATH',webpage)
        
        if self.mode == 'o': #original 
            USE_ONECOLUMN = ONE_COLUMN
            BLCOLUMN = ONE_COLUMNb # this var is used for investigating whether the original profile is of sufficient quality to be used for analysis or class model input.
        elif self.mode == 'b': # best BLH
            USE_ONECOLUMN = ONE_COLUMNb
            BLCOLUMN = ONE_COLUMNb
        elif self.mode == 'u': # best BLH
            USE_ONECOLUMN = ONE_COLUMNu
            BLCOLUMN = ONE_COLUMNu
        elif self.mode == 'd': # best BLH
            USE_ONECOLUMN = ONE_COLUMNd
            BLCOLUMN = ONE_COLUMNd
        else:
            USE_ONECOLUMN = ONE_COLUMN
            BLCOLUMN = ONE_COLUMNb

        lt6000 = (BLCOLUMN['HAGL'] < 6000.)
        lt2500 = (BLCOLUMN['HAGL'] < 2500. + BLHV)
        # print(BLCOLUMN['HAGL'][lt6000])
        # print(BLCOLUMN['HAGL'][lt2500])
        # 
        # print(len(np.where(lt2500)[0]) > 9.) # distance between two points (lower than 2500m) should be smaller than 400 meters

        #print(BLCOLUMN['HAGL'][lt2500])
        PARAMS.insert(0,'OK',
                      ((BLHVe < 200.) and 
                       ( len(np.where(lt6000)[0]) > 5) and
                       (np.array(BLCOLUMN['HAGL'])[-1] >= 6000.) and # the last coordinate had a height higher than 5000.
                       (not len(np.where(pd.isnull(BLCOLUMN['THTA'][lt6000]))[0]) >0 ) and
                       (len(np.where(lt2500)[0]) > 10.) and # distance between two points (lower than 2500m) should be smaller than 400 meters
                       (not len(np.where(pd.isnull(BLCOLUMN['SKNT'][lt6000]))[0]) >0 ) and
                       (not len(np.where(pd.isnull(BLCOLUMN['DRCT'][lt6000]))[0]) >0 ) and
                       (not len(np.where(pd.isnull(BLCOLUMN['PRES'][lt6000]))[0]) >0 ) and
                       (not len(np.where(pd.isnull(BLCOLUMN['QABS'][lt6000]))[0]) >0 ) and
                       (not (len(np.where(np.array(BLCOLUMN['THTA'][lt6000])[2:] <= np.array(BLCOLUMN['THTA'][lt6000])[1:-1])[0]) >0) ) #absolute increasing
                      )
                     )

        PARAMS.insert(0,'theta',np.float(list(BLCOLUMN['THTA'])[1]))
        PARAMS.insert(0,'q',np.float(list(BLCOLUMN['QABS'])[1]))
        PARAMS.insert(0,'u',np.float(list(BLCOLUMN['u_pro'])[1]))  
        PARAMS.insert(0,'v',np.float(list(BLCOLUMN['v_pro'])[1]))
        PARAMS.insert(0,'dtheta',np.float(list(BLCOLUMN['THTA'])[2]-list(BLCOLUMN['THTA'])[1]))
        PARAMS.insert(0,'dq',np.float(list(BLCOLUMN['QABS'])[2]-list(BLCOLUMN['QABS'])[1]))
        PARAMS.insert(0,'du',np.float(list(BLCOLUMN['u_pro'])[2]-list(BLCOLUMN['u_pro'])[1]))
        PARAMS.insert(0,'dv',np.float(list(BLCOLUMN['v_pro'])[2]-list(BLCOLUMN['v_pro'])[1]))
        
        
        PARAMS = PARAMS.T

        
        self.PARAMS = PARAMS
        self.ONE_COLUMN = USE_ONECOLUMN
        # if self.mode == 'o': #original 
        #     self.ONE_COLUMN = ONE_COLUMN
        # elif self.mode == 'b': # best BLH
        #     self.ONE_COLUMN = ONE_COLUMNb
        # elif self.mode == 'u':# upper BLH
        #     self.ONE_COLUMN = ONE_COLUMNu
        # elif self.mode == 'd': # lower BLH
        #     self.ONE_COLUMN=ONE_COLUMNd
        # else:
        #     self.ONE_COLUMN = ONE_COLUMN

