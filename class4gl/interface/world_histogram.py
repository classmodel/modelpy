'''
import numpy as np
import pandas as pd
import sys
import matplotlib
matplotlib.use('TkAgg')

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--path_experiments')#,default='/user/data/gent/gvo000/gvo00090/D2D/data/C4GL/')
parser.add_argument('--path_forcing')#,default='/user/data/gent/gvo000/gvo00090/D2D/data/SOUNDINGS/')
parser.add_argument('--experiments')
parser.add_argument('--c4gl_path_lib')#,default='/user/data/gent/gvo000/gvo00090/D2D/software/CLASS/class4gl/lib')
parser.add_argument('--load_globaldata',default=False)
parser.add_argument('--make_figures',default=None)
parser.add_argument('--show_control_parameters',default=True)
parser.add_argument('--figure_filename',default=None)
parser.add_argument('--figure_filename_2',default=None)
parser.add_argument('--experiments_labels',default=None)
parser.add_argument('--obs_filter',default='True')
args = parser.parse_args()

print('Adding python library:',args.c4gl_path_lib)
sys.path.insert(0, args.c4gl_path_lib)
from interface_multi import c4gl_interface_soundings,get_record_yaml
from class4gl import class4gl_input, data_global,class4gl,units
#from sklearn.metrics import mean_squared_error
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn.apionly as sns
import pylab as pl
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import kde
from scipy.stats import pearsonr                                                
from taylorDiagram import TaylorDiagram
from matplotlib import ticker
import xarray as xr
# import importlib
# importlib.reload(mpl); importlib.reload(plt); importlib.reload(sns)


if args.experiments_labels is None:
    keylabels = args.experiments.strip().split(' ')
else:
    keylabels = args.experiments_labels.strip().split(';')


def rmse(y_actual,y_predicted,z_actual = None, z_predicted = None,filternan_actual = False):
    """ calculated root mean squared error 
        
    
        INPUT:
            y_actual: reference dataset
            y_predicted: predicting dataset
            z_actual: coordinate values of reference dataset
            z_predicted: coordinate values of the predicting dataset
            
            filternan_actual: throw away reference values that have nans
    """
    
    y_actual_temp = np.array(y_actual)
    y_predicted_temp = np.array(y_predicted)
    
    if z_actual is not None:
        z_actual_temp = np.array(z_actual)
    else: 
        z_actual_temp = None
        
    
    if filternan_actual:
        y_actual_temp = y_actual_temp[~np.isnan(y_actual_temp)]
        if z_actual_temp is not None:
            z_actual_temp = z_actual_temp[~np.isnan(y_actual_temp)]
    
    if ((z_actual_temp is not None) or (z_predicted is not None)):    
        if (z_actual_temp is None) or (z_predicted is None):
            raise ValueError('Input z_actual and z_predicted need \
                              to be specified simultaneously.')
        y_predicted_temp = np.interp(z_actual_temp,z_predicted, y_predicted)
    
    else:
        # this catches the situation that y_predicted is a single value (eg., 
        # which is the case for evaluating eg., mixed-layer estimates)
        y_predicted_temp = y_actual_temp*0. + y_predicted_temp
        
    rmse_temp = (y_actual_temp - y_predicted_temp)
    rmse_temp = np.mean(rmse_temp*rmse_temp)
    return np.sqrt(rmse_temp)





# EXPS  =\
# {
# 'GLOBAL_NOAC':    {'sw_ac' : [],'sw_ap': True,'sw_lit': False},
# #'GLOBAL_ADV':{'sw_ac' : ['adv',],'sw_ap': True,'sw_lit': False},
# #'GLOBAL_ITER_NOAC':    {'sw_ac' : [],'sw_ap': True,'sw_lit': False},
# #'GLOBAL_ITER_ADV':    {'sw_ac' : [],'sw_ap': True,'sw_lit': False},
# #'IOPS_ITER_ADV':{'sw_ac' : ['adv',],'sw_ap': True,'sw_lit': False},
# # 'IOPS_W':  {'sw_ac' : ['w',],'sw_ap': True,'sw_lit': False},
# # 'IOPS_AC': {'sw_ac' : ['adv','w'],'sw_ap': True,'sw_lit': False},
# }

if bool(args.load_globaldata):
    # iniitialize global data
    globaldata = data_global()
    # ...  and load initial data pages
    globaldata.load_datasets(recalc=0)
else:
    globaldata = None

c4gldata = {}
for key in args.experiments.strip(' ').split(' '):
    
    c4gldata[key] = c4gl_interface_soundings( \
                      args.path_experiments+'/'+key+'/',\
                      args.path_forcing+'/',\
                      globaldata,\
                      refetch_records=False,
                      obs_filter = (args.obs_filter == 'True')
                                            
                    )
'''

# kgccolors = {
#     'Dfa':['navy','white'],
#     'Cfb':['green','white']       ,
#     'BSk':['tan','black']      ,
#     'Csb':['lightgreen','black'] ,     
#     'Cfa':['darkgreen','white']  ,    
#     'BWh':['orange','black']      ,
#     'Aw' :['pink','black'],
#     'Dwc':['rebeccapurple','white'] ,    
#     'Dfb':['darkviolet','white']    , 
# }





import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
import seaborn as sns
# activate latex text rendering
sns.reset_orig()

data = [220,14.2,150,400,420,100,150,30,60,20,500,]
error = [10, 1, 20, 60, 10,10, 1, 20, 60, 10, 5, ]

ini = c4gldata['GLOBAL_ADV'].frames['stats']['records_all_stations_ini'].set_index(['STNID','dates'])
vals = c4gldata['GLOBAL_ADV'].frames['stats']['records_all_stations_mod']
stats = c4gldata['GLOBAL_ADV'].frames['stats']['records_all_stations_mod_stats']

ini_fc = c4gldata['GLOBAL_ADV_FC'].frames['stats']['records_all_stations_ini'].set_index(['STNID','dates'])
vals_fc = c4gldata['GLOBAL_ADV_FC'].frames['stats']['records_all_stations_mod']
stats_fc = c4gldata['GLOBAL_ADV_FC'].frames['stats']['records_all_stations_mod']

vals_obs = c4gldata['GLOBAL_ADV'].frames['stats']['records_all_stations_obs_afternoon']
stats_obs = c4gldata['GLOBAL_ADV'].frames['stats']['records_all_stations_obs_afternoon']

ini_common_index = ini.index.intersection(ini_fc.index)

vals.index = ini.index
vals = vals.loc[ini_common_index]

stats.index = ini.index
stats = stats.loc[ini_common_index]


vals_obs.index = ini.index
vals_obs = vals_obs.loc[ini_common_index]

stats_obs.index = ini.index
stats_obs = stats_obs.loc[ini_common_index]


ini = ini.loc[ini_common_index]


vals_fc.index = ini_fc.index
vals_fc = vals_fc.loc[ini_common_index]

stats_fc.index = ini_fc.index
stats_fc = stats_fc.loc[ini_common_index]

ini_fc = ini_fc.loc[ini_common_index]

dlat = 10
blat = np.arange(-55.,60.,dlat)[[2,3,4,7,8,9,10,11]]
lats = (blat[1:] + blat[:-1])/2.

fig = plt.figure(figsize=(10,6)) 
variables = ['h','theta',r'q']
labels = [r'$h$',r'$\theta$',r'$q$']
units = ['m','K','g/kg']
xlims = [(-10,3000),(280,310.),(2.5,17.5)]
#var = "theta"
for ivar,var in enumerate(variables):
    ax = fig.add_subplot(1,len(variables),ivar+1,)
    data = []
    data_025 = []
    data_075 = []
    
    data_fc = []
    data_fc_025 = []
    data_fc_075 = []
    
    data_obs = []
    data_obs_025 = []
    data_obs_075 = []
    
    data_ini = []
    data_ini_025 = []
    data_ini_075 = []
    
    lnts = []
    
    for ilat,lat in enumerate(lats):
        print(ilat,lat)
    # 
        query = 'latitude >= '+str(blat[ilat])+' and '+ 'latitude < '+str(blat[ilat+1])
        print(query)
        select = ini.query(query)
        if len(select) >= 7:
            lnts.append(len(select))
            #print(stats.iloc[select.index])
            print(stats.loc[select.index])
            data.append(vals.loc[select.index][var].mean()) 
            data_025.append(vals.loc[select.index][var].quantile(0.25)) 
            data_075.append(vals.loc[select.index][var].quantile(0.75)) 
            data_fc.append(vals_fc.loc[select.index][var].mean()) 
            data_fc_025.append(vals_fc.loc[select.index][var].quantile(0.25)) 
            data_fc_075.append(vals_fc.loc[select.index][var].quantile(0.75)) 
    
            data_obs.append(vals_obs.loc[select.index][var].mean()) 
            data_obs_025.append(vals_obs.loc[select.index][var].quantile(0.25)) 
            data_obs_075.append(vals_obs.loc[select.index][var].quantile(0.75)) 
    
            data_ini.append(ini.loc[select.index][var].mean()) 
            data_ini_025.append(ini.loc[select.index][var].quantile(0.25)) 
            data_ini_075.append(ini.loc[select.index][var].quantile(0.75)) 
        else:
            lnts.append(0)
            data.append(np.nan)
            data_025.append(np.nan)
            data_075.append(np.nan)
            data_fc.append(np.nan)
            data_fc_025.append(np.nan)
            data_fc_075.append(np.nan)
    
            data_obs.append(np.nan)
            data_obs_025.append(np.nan)
            data_obs_075.append(np.nan)
    
    
            data_ini.append(np.nan)
            data_ini_025.append(np.nan)
            data_ini_075.append(np.nan)
    
    data = np.array(data)
    data_025 = np.array(data_025)
    data_075 = np.array(data_075)
    
    data_fc = np.array(data_fc)
    data_fc_025 = np.array(data_fc_025)
    data_fc_075 = np.array(data_fc_075)
    
    
    data_obs = np.array(data_obs)
    data_obs_025 = np.array(data_obs_025)
    data_obs_075 = np.array(data_obs_075)
    
    data_ini = np.array(data_ini)
    data_ini_025 = np.array(data_ini_025)
    data_ini_075 = np.array(data_ini_075)

    if var == 'q':
        data = data*1000.
        data_025 = data_025*1000.
        data_075 = data_075*1000.

        data_fc = data_fc*1000.
        data_fc_025 = data_fc_025*1000.
        data_fc_075 = data_fc_075*1000.
    
        data_obs = data_obs*1000.
        data_obs_025 = data_obs_025*1000.
        data_obs_075 = data_obs_075*1000.

        data_ini = data_ini*1000.
        data_ini_025 = data_ini_025*1000.
        data_ini_075 = data_ini_075*1000.
    

    

    data_left = np.zeros_like(data)
    data_right = np.zeros_like(data)
    select = (data >= data_ini)
    data_left[select] = data[select]
    data_right[select] = data_ini[select]
    data_left[~select] = data_ini[~select]
    data_right[~select] = data[~select]

    data_diff_left = np.zeros_like(data)
    data_diff_right = np.zeros_like(data)
    select = (data_fc >= data)
    data_diff_left[select] = data[select]
    data_diff_right[select] = data_fc[select]
    data_diff_left[~select] = data_fc[~select]
    data_diff_right[~select] = data[~select]


    
    
    #bar = ax.barh(lats, data,blat[1:] - blat[:-1] , align="center", xerr=error)
    
    
    erb = ax.errorbar( data_ini, lats+2.*dlat/8.,xerr=[data_ini-data_ini_025,data_ini_075-data_ini], fmt='s', color='darkgrey',mfc='white', ms=3, mew=1)
    erb = ax.errorbar( data_obs, lats+2.*dlat/8.,xerr=[data_obs-data_obs_025,data_obs_075-data_obs], fmt='s', color='darkgrey',mfc='black', ms=3, mew=1)
    erb = ax.errorbar( data,    lats-1.*dlat/8,xerr=[data-data_025,data_075-data], fmt='s', color='black',mfc='black', ms=3, mew=1)
    er2 = ax.errorbar( data_fc, lats-2*dlat/8.,xerr=[data_fc-data_fc_025,data_fc_075-data_fc], fmt='s', color='blue',mfc='blue', ms=3, mew=1)
    
    ba2 = ax.barh(lats, data_diff_right - data_diff_left  ,(blat[1:] -
                                                            blat[:-1])*0.85 ,
                  align="center", left=data_diff_left,color='red',
                  edgecolor='lightgrey',linewidth=2.)
    bar = ax.barh(lats, data_right - data_left ,(blat[1:] - blat[:-1])*0.85 ,
                  align="center", left=data_left, color='none',edgecolor='black',linewidth=2.)
    # ax.set_yticks(blat)
    # labels = [ w.get_text() for w in ax.get_yticklabels()]
    # ax.set_yticklabels(labels)
    # ax.set_yticks(blat)
    ax.set_yticks([["test",-30.]])
    if ivar == 0:
        plt.legend([r"morning observations",\
                    r"afternoon observations",\
                    r"afternoon control",\
                    r"afternoon $\it{wet}$ "])
    
    #er2 = ax.errorbar( data_fc, lats,xerr=[data_fc-data_fc_025,data_fc_075-data_fc], fmt='s', mfc='blue', ms=3, mew=1)
    
    
    
    # plot = ax.plot(x, data)
    ax.set_yticks(lats)
    ax.set_xlim(xlims[ivar])
    ax.set_ylim((-65.,65.))
    #ax.set_yticklabels(('wt', 'N23PP', 'N23PP/PEKN', 'PEKN', 'N23PP/PEKN/L28F'))
    #ax.set_title(r"Everything in the document can use m$\alpha$th language", y=1.05)
    ax.set_title(labels[ivar],fontsize=17.)
    ax.set_xlabel("["+units[ivar]+"]", labelpad=10,fontsize=15.)
    if ivar == 0:
        ax.set_ylabel("Latitude [°]",labelpad=10,fontsize=15.)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    ax.tick_params(labelsize=15.)
fig.tight_layout()
fig.subplots_adjust(left=0.10,bottom=0.13,right=0.98,top=0.94,wspace=0.24,hspace=0.20)
fig.savefig('test.png',dpi=200)
fig.show()


