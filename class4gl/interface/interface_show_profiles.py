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
parser.add_argument('--load_globaldata',default=False) # load the data needed for the interface
parser.add_argument('--make_figures',default=None)
parser.add_argument('--figure_filename',default=None)
parser.add_argument('--tendencies_revised',default=False)
parser.add_argument('--obs_filter',default="True")
args = parser.parse_args()

print('Adding python library:',args.c4gl_path_lib)
sys.path.insert(0, args.c4gl_path_lib)
from interface_multi import c4gl_interface_soundings,get_record_yaml
from class4gl import class4gl_input, data_global,class4gl,units
#from sklearn.metrics import mean_squared_error
import matplotlib as mpl
import matplotlib.pyplot as plt
#import seaborn.apionly as sns
import pylab as pl
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import kde
from scipy.stats import pearsonr                                                
from taylorDiagram import TaylorDiagram
from matplotlib import ticker
# import importlib
# importlib.reload(mpl); importlib.reload(plt); importlib.reload(sns)





latex = {}
latex['dthetadt'] =  r'$d \theta / dt $'
latex['dqdt'] =      r'$d q / dt $'
latex['dhdt'] =      r'$d h / dt $'

def abline(slope, intercept,axis):
    """Plot a line from slope and intercept"""
    #axis = plt.gca()
    x_vals = np.array(axis.get_xlim())
    y_vals = intercept + slope * x_vals
    axis.plot(x_vals, y_vals, 'k--')

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

if args.load_globaldata:
    # iniitialize global data
    globaldata = data_global()
    # ...  and load initial data pages
    globaldata.load_datasets(recalc=0)
else:
    globaldata = None

c4gldata = []
    
c4gldata.append(c4gl_interface_soundings( \
                  '/data/gent/vo/000/gvo00090/D2D/data/C4GL/20181030/IOPS_ADV_ITER/',\
                  '/data/gent/vo/000/gvo00090/D2D/data/SOUNDINGS/IOPS/',\
                  globaldata,\
                  refetch_records=False,\
                  obs_filter = True,\
                  tendencies_revised = args.tendencies_revised\
                ))
c4gldata.append(c4gl_interface_soundings( \
                  '/data/gent/vo/000/gvo00090/D2D/data/C4GL/20181017_NOON/IOPS_ADV_ITER/',\
                  '/data/gent/vo/000/gvo00090/D2D/data/SOUNDINGS/IOPS_NOON/',\
                  globaldata,\
                  refetch_records=False,
                  obs_filter = False,
                  tendencies_revised = args.tendencies_revised
                ))


profiles_morning            = []
profiles_morning_raw        = []
profiles_noon_obs           = []
profiles_noon_obs_raw       = []
profiles_noon_mod           = []
profiles_afternoon_obs      = []
profiles_afternoon_obs_raw  = []
profiles_afternoon_mod      = []
profiles_location      = []
profiles_morning_datetime      = []
profiles_noon_datetime      = []
profiles_afternoon_datetime      = []

for i in range(15):

    if i == 0:
        location = 'HUMPPA'
        c4gldata[0].sel_station(90000)
        c4gldata[1].sel_station(90000)
    elif i == 7:
        location = 'BLLAST'
        c4gldata[0].sel_station(90001)
        c4gldata[1].sel_station(90001)
    elif i == 11:
        location = 'GOAMAZON'
        c4gldata[0].sel_station(90002)
        c4gldata[1].sel_station(90002)

    # if i == 12:
    #     c4gldata[1].next_record()
    #     c4gldata[1].next_record()

        
    profiles_location.append(str(location))
    profiles_morning_datetime.append(c4gldata[0].frames['profiles']['record_yaml_ini'].pars.ldatetime)
    profiles_morning.append( c4gldata[0].frames['profiles']['record_yaml_ini'].air_ap)
    profiles_morning_raw.append( c4gldata[0].frames['profiles']['record_yaml_ini'].air_balloon)

    profiles_noon_datetime.append(c4gldata[1].frames['profiles']['record_yaml_obs_afternoon'].pars.ldatetime)

    if i in [12,13]:
        profiles_noon_datetime.append(c4gldata[1].frames['profiles']['record_yaml_obs_afternoon'].pars.ldatetime)
        profiles_noon_obs.append(c4gldata[1].frames['profiles']['record_yaml_obs_afternoon'].air_ap)
        profiles_noon_obs_raw.append(c4gldata[1].frames['profiles']['record_yaml_obs_afternoon'].air_balloon)
        profiles_noon_mod.append(c4gldata[1].frames['profiles']['record_yaml_mod'].air_ap)
        profiles_noon_obs[-1].theta = np.nan
        profiles_noon_obs_raw[-1].theta = np.nan
        profiles_noon_mod[-1].theta = np.nan
        profiles_noon_obs[-1].q = np.nan
        profiles_noon_obs_raw[-1].q = np.nan
        profiles_noon_mod[-1].q = np.nan

    else:
        profiles_noon_obs.append( c4gldata[1].frames['profiles']['record_yaml_obs_afternoon'].air_ap)
        profiles_noon_obs_raw.append( c4gldata[1].frames['profiles']['record_yaml_obs_afternoon'].air_balloon)
        profiles_noon_mod.append( c4gldata[1].frames['profiles']['record_yaml_mod'].air_ap)

    profiles_afternoon_datetime.append(c4gldata[0].frames['profiles']['record_yaml_obs_afternoon'].pars.ldatetime)
    profiles_afternoon_obs.append( c4gldata[0].frames['profiles']['record_yaml_obs_afternoon'].air_ap)
    profiles_afternoon_obs_raw.append(c4gldata[0].frames['profiles']['record_yaml_obs_afternoon'].air_balloon)
    profiles_afternoon_mod.append( c4gldata[0].frames['profiles']['record_yaml_mod'].air_ap)

    if i < 14:
        c4gldata[0].next_record()
        if i not in [9,10]:
            c4gldata[1].next_record()




            #axes[varkey].set_title(latex['d'+varkey+'dt']+' '+units_final,fontsize=12.)                                     


'''



numprofiles = len(profiles_afternoon_mod)


per_row = 4.
VARS = ['theta','q']

for var in VARS:
    fig = plt.figure(figsize=(9,8.5))
    for i in range(numprofiles):
        ax = fig.add_subplot(np.ceil(numprofiles/per_row),per_row,i+1)
    
        if var == 'q':
            fac = 1000.
        else:
            fac = 1.
    
        label = profiles_morning_datetime[i].strftime("%H:%m")+'ST'+' ini'
        plt.plot(profiles_morning[i][var].values*fac,profiles_morning[i].z.values,'b--',label=label)
        #plt.plot(profiles_morning_raw[i][var].values,profiles_morning_raw[i].z.values,'b*')
        
        label = profiles_noon_datetime[i].strftime("%H:%m")+'ST'+' obs'
        #label = profiles_morning_datetime[i].hour()+':'+profiles_morning_datetime[i].minute()+'ST')+' ini')
        plt.plot(profiles_noon_obs[i][var].values*fac,profiles_noon_obs[i].z.values,'g--',label=label)
        #plt.plot(profiles_noon_obs_raw[i][var].values,profiles_noon_obs_raw[i].z.values,'g*')
        label = profiles_noon_datetime[i].strftime("%H:%m")+'ST'+' mod'
        plt.plot(profiles_noon_mod[i][var].values*fac,profiles_noon_mod[i].z.values,'g-',label=label)
        
        
        label = profiles_afternoon_datetime[i].strftime("%H:%m")+'ST'+' obs'
        plt.plot(profiles_afternoon_obs[i][var].values*fac, profiles_afternoon_obs[i].z.values,'r--',label=label)
        #plt.plot(profiles_afternoon_obs_raw[i][var].values,profiles_afternoon_obs_raw[i].z.values,'r*')
        label = profiles_afternoon_datetime[i].strftime("%H:%m")+'ST'+' mod'
        plt.plot(profiles_afternoon_mod[i][var].values*fac, profiles_afternoon_mod[i].z.values,'r-',label=label)
        if var == 'q':
            if i >=(numprofiles - per_row):
                ax.legend(loc=2,fontsize=6)
            else:
                ax.legend(loc=1,fontsize=6)
        else:
            ax.legend(loc=2,fontsize=6)
    
    
        ax.set_ylim(0,3500)
        if var == 'theta':
             ax.set_xlim((283.,315.))
        elif var == 'q':
             ax.set_xlim((0.,21.))
        ax.set_title(profiles_location[i]+' '+ str(profiles_morning_datetime[i].date()))
    
        if (np.mod(i,per_row) !=0):
            ax.set_yticklabels([])
        else:
            ax.set_ylabel('$h$ [$\mathrm{m}$]')
    
        if i < (numprofiles - per_row):
            ax.set_xticklabels([])
        else:
           if var == 'q':
               units_final = r'$q$ [$\mathrm{g\, kg^{-1}}$]'
           elif var == 'theta':
               units_final = r'$\theta$ [$\mathrm{K}$]'
           ax.set_xlabel(units_final)




    #ax.tick_params(
    #axis='x',          # changes apply to the x-axis
    #which='both',      # both major and minor ticks are affected
    #left=(np.mod(i,per_row) ==0),      # ticks along the bottom edge are off
    #top=False,         # ticks along the top edge are off
    #labelbottom= (i==np.ceil(numprofiles/per_row)) # labels along the bottom edge are off
    #)

    fig.tight_layout()
    fig.subplots_adjust(wspace=0.05,right=0.99,hspace=0.25)
    fig.show()



# if bool(args.make_figures):
#     fig = plt.figure(figsize=(10,7))   #width,height
#     i = 1                                                                           
#     axes = {}         
#     axes_taylor = {}         
#     
#     #colors = ['r','g','b','m']
#     colors = ['k']
#     symbols = ['^','x','+']
#     dias = {}
#     
#     for varkey in ['h','theta','q']:                                                    
#         axes[varkey] = fig.add_subplot(2,3,i)                                       
#         #axes_taylor[varkey] = fig.add_subplot(2,3,i+3)                                       
#     
#         #print(obs.std())
#         obs = c4gldata[args.experiments.strip().split()[0]].frames['stats']['records_all_stations_obs_afternoon_stats']['d'+varkey+'dt']
#         STD_OBS = obs.std()
#         dias[varkey] =  TaylorDiagram(1., srange=[0.0,1.7],fig=fig, rect=(230+i+3),label='Reference')
#         dias[varkey]._ax.axis["left"].label.set_text(\
#             "Normalized standard deviation")
#         if i == 1:
#             axes[varkey].annotate('Normalized standard deviation',\
#                         xy= (0.05,0.36),
#                         color='black',
#                         rotation=90.,
#                         xycoords='figure fraction',
#                         weight='normal',
#                         fontsize=10.,
#                         horizontalalignment='center',
#                         verticalalignment='center' ,
#                         #bbox={'edgecolor':'black',
#                         #      'boxstyle':'circle',
#                         #      'fc':koeppen.color,
#                         #      'alpha':1.0}
#                        )
#         # dias[varkey]._ax.axis["left"].axis.set_ticks(np.arange(0.,2.,0.25))
#         # dias[varkey]._ax.axis["left"].axis.set_major_locator(np.arange(0.,2.,0.25))
#         #dias[varkey]._ax.axis["left"].axis.set_ticks(np.arange(0.,2.,0.25))
#         # Q95 = obs.quantile(0.95)
#         # Q95 = obs.quantile(0.90)
#         # Add RMS contours, and label them
#         contours = dias[varkey].add_contours(levels=5, colors='0.5') # 5 levels
#         dias[varkey].ax.clabel(contours, inline=1, fontsize=10, fmt='%.1f')
#         #dia._ax.set_title(season.capitalize())
#     
#         dias[varkey].add_grid()
#     
#     
#         #dia.ax.plot(x99,y99,color='k')
#     
#         
#         for ikey,key in enumerate(args.experiments.strip(' ').split(' ')):
#             mod = c4gldata[key].frames['stats']['records_all_stations_mod_stats']['d'+varkey+'dt']
#             obs = c4gldata[key].frames['stats']['records_all_stations_obs_afternoon_stats']['d'+varkey+'dt']
#             x, y = obs.values,mod.values
#             print(key,len(obs.values))
#     
#             #scores
#             PR = pearsonr(mod,obs)[0]
#             RMSE = rmse(obs,mod)                                               
#             BIAS = np.mean(mod) - np.mean(obs)
#             STD = mod.std()
#             
#             fit = np.polyfit(x,y,deg=1)
# 
#             if varkey == 'q':
#                 axes[varkey].plot(x, fit[0] * x + fit[1],\
#                                   color=colors[ikey],alpha=0.8,lw=2,\
#                                   label=key+", "+\
#                            'RMSE = '+format((RMSE*1000.),'0.2f')+r'$\,  \mathrm{g\,  kg^{-1}\,  h^{-1}}$'+ '\n'+\
#                            'Bias = '+format((BIAS*1000.),'0.2f')+r'$\,  \mathrm{g\,  kg^{-1}\,  h^{-1}}$'+' \n'+\
#                            r'$R$ = '+format(PR,'0.2f') )
# 
# 
#             elif varkey == 'h':
#                 axes[varkey].plot(x, fit[0] * x + fit[1],\
#                                   color=colors[ikey],alpha=0.8,lw=2,\
#                                   label=key+", "+\
#                             'RMSE = '+format(RMSE,'0.1f')+r'$\,  \mathrm{m\, h^{-1}}$'+'\n'+\
#                             'Bias = '+format(BIAS,'0.1f')+r'$\,  \mathrm{m\, h^{-1}}$'+'\n'+\
#                             r'$R$ = '+format(PR,'0.2f'))
#             else: #theta
#                 axes[varkey].plot(x, fit[0] * x + fit[1],\
#                                   color=colors[ikey],alpha=0.8,lw=2,\
#                                   label=key+", "+\
#                             'RMSE = '+format(RMSE,'0.3f')+r'$\, \mathrm{K\, h^{-1}}$'+'\n'+\
#                             'Bias = '+format(BIAS,'0.3f')+r'$\, \mathrm{K\, h^{-1}}$'+'\n'+\
#                             r'$R$ = '+format(PR,'0.2f'))
# 
#             if varkey == 'q':
#                 annotate_text = \
#                                'RMSE = '+format((RMSE*1000.),'0.2f')+r'$\,  \mathrm{g\,  kg^{-1}\,  h^{-1}}$'+ '\n'+\
#                                'Bias = '+format((BIAS*1000.),'0.2f')+r'$\,  \mathrm{g\,  kg^{-1}\,  h^{-1}}$'+' \n'+\
#                                r'$R$ = '+format(PR,'0.2f')
#                 ann = axes[varkey].annotate(annotate_text, xy=(0.95, .05 ), xycoords='axes fraction',fontsize=9,
#        horizontalalignment='right', verticalalignment='bottom' ,
#         bbox={'edgecolor':'black',
#                           'fc':'white',  
#                               'boxstyle':'square',
#                               'alpha':0.8}
#                                        )
#             elif varkey == 'h':
#                 annotate_text = \
#                                 'RMSE = '+format(RMSE,'0.1f')+r'$\,  \mathrm{m\, h^{-1}}$'+'\n'+\
#                                 'Bias = '+format(BIAS,'0.1f')+r'$\,  \mathrm{m\, h^{-1}}$'+'\n'+\
#                                 r'$R$ = '+format(PR,'0.2f')
#                 ann = axes[varkey].annotate(annotate_text, xy=(0.05, .97 ), xycoords='axes fraction',fontsize=9,
#        horizontalalignment='left', verticalalignment='top' ,
#         bbox={'edgecolor':'black',
#                           'fc':'white',  
#                               'boxstyle':'square',
#                               'alpha':0.8}
#                                        )
#             else:
#                 annotate_text = \
#                                 'RMSE = '+format(RMSE,'0.3f')+r'$\, \mathrm{K\, h^{-1}}$'+'\n'+\
#                                 'Bias = '+format(BIAS,'0.3f')+r'$\, \mathrm{K\, h^{-1}}$'+'\n'+\
#                                 r'$R$ = '+format(PR,'0.2f')
# 
#                 ann = axes[varkey].annotate(annotate_text, xy=(0.05, .97 ), xycoords='axes fraction',fontsize=9,
#        horizontalalignment='left', verticalalignment='top' ,
#         bbox={'edgecolor':'black',
#                           'fc':'white',  
#                               'boxstyle':'square',
#                               'alpha':0.8}
#                                        )
# 
# 
# 
# 
#             
#             # print(STD)
#             # print(PR)
#             dias[varkey].add_sample(STD/STD_OBS, PR,
#                            marker='o', ms=5, ls='',
#                            #mfc='k', mec='k', # B&W
#                            mfc=colors[ikey], mec=colors[ikey], # Colors
#                            label=key)
#     
#         # put ticker position, see
#         # https://matplotlib.org/examples/ticks_and_spines/tick-locators.html 
#         # dia.ax.axis['bottom'].
#         # dia.ax.axis['left'].
#         # dia.ax.axis['left'].
#     
#         i += 1
#     
#     i = 0
#     for varkey in ['h','theta','q']:                                                    
#         for ikey,key in enumerate(args.experiments.strip(' ').split(' ')):
#             istation = 0
#             for icurrent_station,current_station in c4gldata[key].frames['worldmap']['stations'].table.iterrows():
#                 indices =  (c4gldata[key].frames['stats']['records_all_stations_index'].get_level_values('STNID') == current_station.name)
#                 station_mod = c4gldata[key].frames['stats']['records_all_stations_mod_stats']['d'+varkey+'dt'].iloc[indices]
#                 station_obs = c4gldata[key].frames['stats']['records_all_stations_obs_afternoon_stats']['d'+varkey+'dt'].iloc[indices]
#     
#                 axes[varkey].scatter(station_obs,station_mod,marker=symbols[istation],color=colors[ikey])
#                          #  label=key+", "+\
#                          #                    'R = '+str(round(PR[0],3))+', '+\
#                          #                    'RMSE = '+str(round(RMSE,5))+', '+\
#                          #                    'BIAS = '+str(round(BIAS,5)),s=1.,color=colors[ikey])
#     
#     
#     
#             # # pl.scatter(obs,mod,label=key+", "+\
#             # #                              'R = '+str(round(PR[0],3))+', '+\
#             # #                              'RMSE = '+str(round(RMSE,5))+', '+\
#             # #                              'BIAS = '+str(round(BIAS,5)),s=1.,color=colors[ikey])
#                 
#                 dias[varkey].add_sample(station_mod.std()/station_obs.std(),
#                                pearsonr(station_mod,station_obs)[0],#annotate=symbols[istation],
#                                marker=symbols[istation], ms=5, ls='',
#                                mfc='k', mec='k', # B&W
#                                #mfc=colors[ikey], mec=colors[ikey], # Colors
#                                label=key)
# 
#                 istation += 1
#     
#             if varkey == 'q':
#                 units_final = r'[$g\, kg^{-1}\, h^{-1}$]'
#             elif varkey == 'theta':
#                 units_final = r'[$K\, h^{-1}$]'
#             elif varkey == 'h':
#                 units_final = r'[$m\, h^{-1}$]'
#     
#             axes[varkey].set_xlabel('Observed')     
#             axes[varkey].set_title(latex['d'+varkey+'dt']+' '+units_final,fontsize=12)                                     
# 
# 
#         # if varkey == 'q':
#         #     print('get_xlim not working well...STRANGE')
#         #     limits =  [np.percentile(nani,1),np.percentile(nani,99)]
#         # else:
#         #     limits =  [np.percentile(nani,1.0),np.percentile(nani,99.0)]
# 
# 
#         if i==0:                                    
#             axes[varkey].set_ylabel('Modelled')                                            
#         i +=1
#           
#         axes[varkey].set_aspect('equal')
#         low  = c4gldata[key].frames['stats']['records_all_stations_mod_stats']['d'+varkey+'dt'].min()
#         high  = c4gldata[key].frames['stats']['records_all_stations_mod_stats']['d'+varkey+'dt'].max()
# 
#         low  = np.min([low,c4gldata[key].frames['stats']['records_all_stations_obs_afternoon_stats']['d'+varkey+'dt'].min()])
#         high  = np.max([high,c4gldata[key].frames['stats']['records_all_stations_obs_afternoon_stats']['d'+varkey+'dt'].max()])
# 
#         low = low - (high - low)*0.1
#         high = high + (high - low)*0.1
#         axes[varkey].set_xlim([low,high])
#         axes[varkey].set_ylim([low,high])
#         abline(1,0,axis=axes[varkey])
#         if varkey == 'q':
#             ticks = ticker.FuncFormatter(lambda x, pos:
#                                          '{0:g}'.format(x*1000.))
#             axes[varkey].xaxis.set_major_formatter(ticks)
#             axes[varkey].yaxis.set_major_formatter(ticks)
#     
#     
#     # # legend for different forcing simulations (colors)
#     # ax = fig.add_axes([0.05,0.00,0.15,0.15]) #[*left*, *bottom*, *width*,    *height*]
#     # leg = []
#     # for ikey,key in enumerate(args.experiments.strip(' ').split(' ')):
#     #     leg1, = ax.plot([],colors[ikey]+'s' ,markersize=10)
#     #     leg.append(leg1)
#     # ax.axis('off')
#     # #leg1 =
#     # ax.legend(leg,list(args.experiments.strip(' ').split(' ')),loc=2,fontsize=10)
#     
#     
#     # legend for different stations (symbols)
#     ax = fig.add_axes([0.08,-0.02,0.15,0.15]) #[*left*, *bottom*, *width*,    *height*]
#     leg = []
#     isymbol = 0
#     for icurrent_station,current_station in c4gldata[key].frames['worldmap']['stations'].table.iterrows():
#         leg1, = ax.plot([],'k'+symbols[isymbol] ,markersize=10)
#         leg.append(leg1)
#         isymbol += 1
#     
#     # symbol for all stations
#     leg1, = ax.plot([],'ko',markersize=10)
#     leg.append(leg1)
#     
#     
#     ax.axis('off')
#     ax.legend(leg,['HUMPPA','BLLAST','GOAMAZON','All'],loc=2,fontsize=10,ncol=4)
#     
#     
#     fig.subplots_adjust(top=0.95,bottom=0.20,left=0.08,right=0.94,hspace=0.28,wspace=0.29)
#     
#     
#     #pl.legend(leglist,('EMI:WOC','EMI:MED','EMI:BEC'),loc=2,fontsize=16,prop={'family':
#     #figfn = '/user/data/gent/gvo000/gvo00090/D2D/archive/report/iops_eval_report.png'
#     
#     if args.figure_filename is not None:
#         fig.savefig(args.figure_filename,dpi=200); print("Image file written to:",args.figure_filename)
#         fig.savefig(args.figure_filename.replace('png','pdf')); print("Image file written to:", args.figure_filename)
#     fig.show()  












