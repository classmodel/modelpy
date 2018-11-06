
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
parser.add_argument('--tendencies_revised',default=False)
parser.add_argument('--obs_filter',default='True')
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


if args.experiments_labels is None:
    keylabels = args.experiments.strip().split(' ')
else:
    keylabels = args.experiments_labels.strip().split(';')



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
                      args.path_forcing,\
                      globaldata,\
                      refetch_records=False,
                      tendencies_revised = args.tendencies_revised,
                      obs_filter = (args.obs_filter == 'True')
                    )

if args.make_figures:
    # the lines below activate TaylorPlots but it is disabled for now
    fig = plt.figure(figsize=(10,7))   #width,height
    i = 1                                                                           
    axes = {}         
    axes_taylor = {}         
    
    colors = ['r','g','b','m','y','c']
    symbols = ['*','x','+','o']
    dias = {}
    
    varkeys = ['h','theta','q']
    for varkey in varkeys:                                                    
        axes[varkey] = fig.add_subplot(2,3,i)                                       
        #axes_taylor[varkey] = fig.add_subplot(2,3,i+3)                                       
    
        #print(obs.std())
        dias[varkey] =  TaylorDiagram(1., srange=[0.0,1.7],fig=fig, rect=(230+i+3),label='Reference')
        if i == 0:
            dias[varkey]._ax.axis["left"].label.set_text(\
                "Standard deviation (model) / Standard deviation (observations)")
            # dias[varkey]._ax.axis["left"].axis.set_ticks(np.arange(0.,2.,0.25))
            # dias[varkey]._ax.axis["left"].axis.set_major_locator(np.arange(0.,2.,0.25))
        #dias[varkey]._ax.axis["left"].axis.set_ticks(np.arange(0.,2.,0.25))
        # Q95 = obs.quantile(0.95)
        # Q95 = obs.quantile(0.90)
        # Add RMS contours, and label them
        contours = dias[varkey].add_contours(levels=5, colors='0.5') # 5 levels
        dias[varkey].ax.clabel(contours, inline=1, fontsize=10, fmt='%.1f')
        #dia._ax.set_title(season.capitalize())
    
        dias[varkey].add_grid()
    
    
        #dia.ax.plot(x99,y99,color='k')
    
        
        for ikey,key in enumerate(args.experiments.strip(' ').split(' ')):
            # cc = c4gldata[key].frames['stats']['records_all_stations_ini']['cc']
            # clearsky = (cc < 0.05)
            # mod = c4gldata[key].frames['stats']['records_all_stations_mod_stats'].loc[clearsky]['d'+varkey+'dt']
            # obs = c4gldata[key].frames['stats']['records_all_stations_obs_afternoon_stats'].loc[clearsky]['d'+varkey+'dt']
            mod = c4gldata[key].frames['stats']['records_all_stations_mod_stats']['d'+varkey+'dt']
            obs = c4gldata[key].frames['stats']['records_all_stations_obs_afternoon_stats']['d'+varkey+'dt']
            x, y = obs.values,mod.values
            print(key,len(obs.values))
    
            STD_OBS = obs.std()
            #scores
            PR = pearsonr(mod,obs)[0]
            RMSE = rmse(obs,mod)                                               
            BIAS = np.mean(mod) - np.mean(obs)
            STD = mod.std()
            
            # fit = np.polyfit(x,y,deg=1)
            # axes[varkey].plot(x, fit[0] * x + fit[1],\
            #                   color=colors[ikey],alpha=0.8,lw=2,\
            #                   label=key+", "+\
            #                               'R = '+str(round(PR,3))+', '+\
            #                               'RMSE = '+str(round(RMSE,5))+units['d'+varkey+'dt']+', '+\
            #                               'BIAS = '+str(round(BIAS,5))+units['d'+varkey+'dt'] )
            # axes[varkey].legend(fontsize=5)
            
            # print(STD)
            # print(PR)
            dias[varkey].add_sample(STD/STD_OBS, PR,
                           marker='o', ms=5, ls='',
                           #mfc='k', mec='k', # B&W
                           mfc=colors[ikey], mec=colors[ikey], # Colors
                           label=key)
    
        # put ticker position, see
        # https://matplotlib.org/examples/ticks_and_spines/tick-locators.html 
        # dia.ax.axis['bottom'].
        # dia.ax.axis['left'].
        # dia.ax.axis['left'].
    
        i += 1
    
    i = 0
    for varkey in ['h','theta','q']:                                                    
        ikey = 0
        key = list(args.experiments.strip().split(' '))[ikey]
        # cc = c4gldata[key].frames['stats']['records_all_stations_ini']['cc']
        # clearsky = (cc < 0.05)
    
        # mod = c4gldata[key].frames['stats']['records_all_stations_mod_stats'].loc[clearsky]['d'+varkey+'dt']
        # obs = c4gldata[key].frames['stats']['records_all_stations_obs_afternoon_stats'].loc[clearsky]['d'+varkey+'dt']
        mod = c4gldata[key].frames['stats']['records_all_stations_mod_stats']['d'+varkey+'dt']
        obs = c4gldata[key].frames['stats']['records_all_stations_obs_afternoon_stats']['d'+varkey+'dt']
    
    
        nbins=40       
        x, y = obs.values,mod.values
        
        xi, yi = np.mgrid[x.min():x.max():nbins*1j, y.min():y.max():nbins*1j]
        zi = np.zeros_like(xi)*np.nan       
        for ibin in range(nbins):
            xmin = x.min() + ibin * (x.max() - x.min())/nbins
            xmax = xmin + (x.max() - x.min())/nbins
            in_bin = ((x >= xmin) & (x < xmax))
            ybin = y[in_bin]
            xbin = x[in_bin]
            if len(ybin) > 20:
                k = kde.gaussian_kde((ybin))
                zi[ibin] = k(np.vstack([yi[ibin].flatten()]))
        zi = zi/np.sum(zi,axis=1)[:,np.newaxis]
        zi_int = zi.cumsum(axis=1) 
                     #  label=key+", "+\
                     #                    'R = '+str(round(PR[0],3))+', '+\
                     #                    'RMSE = '+str(round(RMSE,5))+', '+\
                     #                    'BIAS = '+str(round(BIAS,5)),s=1.,color=colors[ikey])
        axes[varkey].contour(xi, yi, zi_int.reshape(xi.shape),levels=[0.16,0.5,0.84] ,
                colors=['darkred','lightgreen','darkred'],linewidths=[1,2,1])
        axes[varkey].contourf(xi, yi, zi_int.reshape(xi.shape),levels=[0.16,0.84] ,
                colors=['darkred'],alpha=0.5,)
        nanxi = xi[zi != np.nan]
        axes[varkey].set_xlim((nanxi.min(),nanxi.max()))
        axes[varkey].set_ylim((nanxi.min(),nanxi.max()))
        print(varkey,(nanxi.min(),nanxi.max()))
    
    
        latex = {}
        latex['dthetadt'] =  r'$d \theta / dt $'
        latex['dqdt'] =      r'$d q / dt $'
        latex['dhdt'] =      r'$d h / dt $'
    
        axes[varkey].set_xlabel('observations')     
        axes[varkey].set_title(latex['d'+varkey+'dt']+' ['+units['d'+varkey+'dt']+']')                                     
    
        PR = pearsonr(mod,obs)[0]
        RMSE = rmse(obs,mod)                                               
        BIAS = np.mean(mod) - np.mean(obs)
        STD = mod.std()
    
        axes[varkey].scatter(obs,mod, label='(only) '+key+", "+\
                                      'R = '+str(round(PR,3))+', '+\
                                      'RMSE = '+str(round(RMSE,5))+units['d'+varkey+'dt']+', '+\
                                      'BIAS = '+str(round(BIAS,5))+units['d'+varkey+'dt'] ,\
                             s=0.1,alpha=0.14,color='k')
        axes[varkey].legend(fontsize=5)
        



        axes[varkey].set_xlabel('observations')     
        if i==0:                                    
            axes[varkey].set_ylabel('model')                                            
        abline(1,0,axis=axes[varkey])
        i +=1
    
    
    
    # legend for different forcing simulations (colors)
    ax = fig.add_axes([0.05,0.00,0.15,0.15]) #[*left*, *bottom*, *width*,    *height*]
    leg = []
    for ikey,key in enumerate(args.experiments.strip().split(' ')):
        leg1, = ax.plot([],colors[ikey]+'o' ,markersize=10)
        leg.append(leg1)
    ax.axis('off')
    #leg1 =
    ax.legend(leg,list(args.experiments.strip().split(' ')),loc=2,fontsize=10)
    
    
    # # legend for different stations (symbols)
    # ax = fig.add_axes([0.25,0.00,0.15,0.15]) #[*left*, *bottom*, *width*,    *height*]
    # leg = []
    # isymbol = 0
    # for icurrent_station,current_station in c4gldata[key].frames['worldmap']['stations'].table.iterrows():
    #     leg1, = ax.plot([],'k'+symbols[isymbol] ,markersize=10)
    #     leg.append(leg1)
    #     isymbol += 1
    # 
    # # symbol for all stations
    # leg1, = ax.plot([],'ko',markersize=10)
    # leg.append(leg1)
    
    
    # ax.axis('off')
    # ax.legend(leg,['HUMPPA','BLLAST','GOAMAZON','All'],loc=2,fontsize=10)
    
    
    fig.subplots_adjust(top=0.95,bottom=0.20,left=0.08,right=0.94,hspace=0.28,wspace=0.29)
    
    
    #pl.legend(leglist,('EMI:WOC','EMI:MED','EMI:BEC'),loc=2,fontsize=16,prop={'family':
    # figfn = '/user/data/gent/gvo000/gvo00090/D2D/archive/report/global_eval_report_cs.png'
    # fig.savefig(figfn,dpi=200); print("Image file written to:", figfn)
    
    if args.figure_filename is not None:
        fig.savefig(args.figure_filename,dpi=200); print("Image file written to:",args.figure_filename)
    fig.show()  

    if bool(args.show_control_parameters):

        import seaborn as sns

        pkmn_type_colors = [
                                            '#A0A0A0',  # Poison
                                            '#78C850',  # Grass
                                            '#F08030',  # Fire
                                            '#6890F0',  # Water
                                            '#F08030',  # Fire
                                            '#C03028',  # Fighting
                                            '#F85888',  # Psychic
                                            '#A8B820',  # Bug
                                            '#A8A878',  # Normal
                                            '#F8D030',  # Electric
                                            '#E0C068',  # Ground
                                            '#EE99AC',  # Fairy
                                            '#B8A038',  # Rock
                                            '#705898',  # Ghost
                                            '#98D8D8',  # Ice
                                            '#7038F8',  # Dragon
                                           ]



        sns.set_style('whitegrid')
        #sns.set()
        fig = pl.figure(figsize=(11,7))
        i = 1
        axes = {}
        data_all = pd.DataFrame()
        data_input = pd.DataFrame()
        
        
        
        # #for varkey in ['theta','q']:     
        # EF =\
        #     c4gldata[key].frames['stats']['records_all_stations_ini'].BR/(1.+\
        #     c4gldata[key].frames['stats']['records_all_stations_ini'].BR)
        # EF[EF<0] = np.nan
        # EF[EF>1] = np.nan
        
        # c4gldata[key].frames['stats']['records_all_stations_ini']['EF'] = EF
        
        ikey = 0
        key = list(args.experiments.strip().split(' '))[ikey]
        data_all = pd.DataFrame()

        tempdatamodstats = pd.DataFrame(c4gldata[key].frames['stats']['records_all_stations_obs_afternoon_stats'].copy())
        tempdatamodstats["source"] = "Soundings"
        tempdatamodstats["source_index"] = "Soundings"

        ini_ref = pd.DataFrame(c4gldata[key].frames['stats']['records_all_stations_ini'].copy())
        tempdataini_this = pd.DataFrame(ini_ref.copy())

        tempdatamodstats['dates']= tempdataini_this.ldatetime.dt.date
        tempdatamodstats['STNID']= tempdataini_this.STNID
        tempdatamodstats['source']= "Soundings"
        tempdatamodstats['source_index']= "Soundings"
        tempdatamodstats.set_index(['source_index','STNID','dates'],inplace=True)
        #print('hello')

        tempdataini = pd.DataFrame(ini_ref)
        tempdataini["source"] = "Soundings"
        tempdataini["source_index"] = "Soundings"
        tempdataini = tempdataini.set_index(['source_index','STNID','dates'])
        #print('hello2')


        data_all = pd.concat([data_all,tempdatamodstats],axis=0)
        data_input = pd.concat([data_input,tempdataini],axis=0)
        #print(data_input.shape)
        #print(data_all.shape)

            
        for ikey,key in enumerate(list(args.experiments.strip().split(' '))):
            keylabel = keylabels[ikey]

            tempdatamodstats = pd.DataFrame(c4gldata[key].frames['stats']['records_all_stations_mod_stats'].copy())
            tempdataini_this= pd.DataFrame(c4gldata[key].frames['stats']['records_all_stations_ini'].copy())
            tempdatamodstats['dates']= tempdataini_this.ldatetime.dt.date
            tempdatamodstats['STNID']= tempdataini_this.STNID
            tempdatamodstats['source']= keylabel
            tempdatamodstats['source_index']= keylabel
            tempdatamodstats.set_index(['source_index','STNID','dates'],inplace=True)
            #print('hello')


            tempdataini = pd.DataFrame(ini_ref.copy())
            tempdataini["source"] = keylabel
            tempdataini["source_index"] = keylabel
            tempdataini = tempdataini.set_index(['source_index','STNID','dates'])
    

            #print('hello2')
            index_intersect = tempdataini.index.intersection(tempdatamodstats.index)
            #print('hello3')

            tempdataini = tempdataini.loc[index_intersect]
            #print('hello4')
            tempdatamodstats = tempdatamodstats.loc[index_intersect]
            #print('hello5')


            # data[varkey] = tempdatamodstats['d'+varkey+'dt']
            data_all = pd.concat([data_all,tempdatamodstats],axis=0)
            data_input = pd.concat([data_input, tempdataini],axis=0)
            #print(data_input.shape)
            #print(data_all.shape)

        data_input.cc = data_input.cc.clip(0.,+np.inf)

        for varkey in ['h','theta','q']:
            varkey_full = 'd'+varkey+'dt ['+units[varkey]+'/h]'
            data_all = data_all.rename(columns={'d'+varkey+'dt':varkey_full})
            
        data_input['advt_tropo'] = data_input['advt_tropo'] * 3600.
        data_all['advt_tropo'] = data_input['advt_tropo']
            #print(data_input.shape)
            #print(data_all.shape)
        #print('hello6')
        #print(data_all.columns)
        #print('hello7')
        i = 1
        for varkey in ['h','theta','q']:
            input_keys =['wg','advt_tropo']
            for input_key in input_keys:
                varkey_full = 'd'+varkey+'dt ['+units[varkey]+'/h]'

                #print('hello8')
                #print(data_input.shape)
                #print(data_all.shape)
                units['advt_tropo'] = 'K/h'
                input_key_full = input_key + "["+units[input_key]+"]"
                data_all[input_key_full] = pd.cut(x=data_input[input_key].values,bins=8,precision=2)
                data_input[input_key_full] = pd.cut(x=data_input[input_key].values,bins=8,precision=2,)
                #print('hello9')
                #print(data_input.shape)
                #print(data_all.shape)
                
                qvalmax = data_all[varkey_full].quantile(0.999)
                qvalmin = data_all[varkey_full].quantile(0.001)
                select_data = (data_all[varkey_full] >= qvalmin) & (data_all[varkey_full] < qvalmax)
                #print('hello11')
                data_all = data_all[select_data]
                #print('hello12')
                data_input = data_input[select_data.values]
                #print('hello13')
                #print(data_input.shape)
                #print(data_all.shape)
                #print('hello10')
                
                sns.set(style="ticks", palette="pastel")
                ax = fig.add_subplot(3,len(input_keys),i)
                #sns.violinplot(x=input_key_full,y=varkey_full,data=data_all,hue='source',linewidth=2.,palette="muted",split=True,inner='quart') #,label=key+", R = "+str(round(PR[0],3)),data=data)       
                
                #ax.set_title(input_key_full)
                sb = sns.boxplot(x=input_key_full, y=varkey_full, hue="source",
                                 palette=pkmn_type_colors,
                                # palette=["m", "g",'r','b'],
                                 linewidth=1.2, data=data_all,sym='')
                if i ==1:
                     plt.legend(loc='upper right',fontsize=7.)
                else:
                     ax.get_legend().set_visible(False)
                #     plt.legend('off')
                if i >= 5:
                    #ax.set_xticklabels(labels=['['+str(i)+','+str(i+1)+'[' for i in list(range(0,7))]+['[7,8]'])

                    ax.set_xticklabels(labels=ax.get_xticklabels(),rotation=45.,ha='right')
                else:
                    ax.set_xticklabels([])
                    ax.set_xlabel('')

                if np.mod(i,len(input_keys)) != 0:
                    ax.set_yticklabels([])
                    ax.set_ylabel('')

                if varkey == 'q':
                    ticks = ticker.FuncFormatter(lambda x, pos:
                                                 '{0:g}'.format(x*1000.))
                    #ax.xaxis.set_major_formatter(ticks)
                    ax.yaxis.set_major_formatter(ticks)

                    ax.set_ylabel(latex['d'+varkey+'dt']+' ['+r'$10^{-3} \times $'+units['d'+varkey+'dt']+']')        
                else:
                    ax.set_ylabel(latex['d'+varkey+'dt']+' ['+units['d'+varkey+'dt']+']')        


                for j,artist in enumerate(ax.artists):
                    if np.mod(j,len(list(args.experiments.strip().split(' ')))+1) !=0:
                        # Set the linecolor on the artist to the facecolor, and set the facecolor to None
                        #print(j,artist)
                        col = artist.get_facecolor()
                        #print(j,artist)
                        artist.set_edgecolor(col)
                        #print(j,artist)
                        artist.set_facecolor('None')
                
                        # Each box has 6 associated Line2D objects (to make the whiskers, fliers, etc.)
                        # Loop over them here, and use the same colour as above
                        
                        for k in range(j*5,j*5+5):
                            line = ax.lines[k]
                            line.set_color(col)
                            line.set_mfc(col)
                            line.set_mec(col)
                
                # Also fix the legend
                j = 0
                for legpatch in ax.get_legend().get_patches():
                    if j > 0:

                        col = legpatch.get_facecolor()
                        legpatch.set_edgecolor(col)
                        legpatch.set_facecolor('None')
                    j +=1




                #ax.grid()
                #sns.despine(offset=10, trim=True)
                i +=1
        fig.tight_layout()
        fig.subplots_adjust( bottom=0.12,left=0.15,top=0.99,right=0.99,wspace=0.05,hspace=0.05,)
        if args.figure_filename_2 is not None:
            fig.savefig(args.figure_filename_2,dpi=200); print("Image file written to:", args.figure_filename_2)
        fig.show()



