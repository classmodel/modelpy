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
import numpy as np

data = [220,14.2,150,400,420]
error = [10, 1, 20, 60, 10]
x = [i + .5 for i in range(5)]

fig, ax = plt.subplots()
bar = ax.barh(x, data,0.1 +np.arange(len(data))*0.9/len(data), align="center", yerr=error)
plot = ax.plot(x, data)
ax.set_xticks(x)
ax.set_xticklabels(('wt', 'N23PP', 'N23PP/PEKN', 'PEKN', 'N23PP/PEKN/L28F'))
ax.set_title(r"Everything in the document can use m$\alpha$th language",
             y=1.05)
ax.set_ylabel(r"Rate (s$^{-1}$)", labelpad=10)
ax.set_xlabel("Mutant",labelpad=10)
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
plt.savefig('test.png')
plt.show()


