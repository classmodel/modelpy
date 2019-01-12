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

sns.reset_orig()



lookup_symbols= {
 'A':'equatorial',
 'B':'arid',
 'C':'warm temperate',
 'D':'snow',
 'E':'polar',
 'W':'desert',
 'S':'steppe',
 'f':'fully humid',
 's':'summer dry',
 'w':'winter dry',
 'm':'monsoonal',
 'h':'hot arid',
 'k':'cold arid',
 'a':'hot summer',
 'b':'warm summer',
 'c':'cool summer',
 'd':'extremely continental',
 'F':'polar frost',
 'T':'polar tundra',
 } 

key = args.experiments.strip(' ').split(' ')[0]
xrkoeppen = xr.open_dataset('/user/data/gent/gvo000/gvo00090/EXT/data/KOEPPEN/Koeppen-Geiger.nc')
koeppenlookuptable = pd.DataFrame()
koeppenlookuptable['KGCID'] = pd.Series(xrkoeppen['KGCID'])

from matplotlib.patches import Rectangle,FancyBboxPatch

def abline(slope, intercept,axis):
    """Plot a line from slope and intercept"""
    #axis = plt.gca()
    xlim = axis.get_xlim()
    ylim = axis.get_ylim()
    x_vals = np.array([np.min([xlim[0],ylim[0]]),np.max([xlim[1],ylim[1]])])
    y_vals = intercept + slope * x_vals
    axis.plot(x_vals, y_vals, 'k--')

KGCID=    ['Af', 'Am', 'As', 'Aw', 'BSh', 'BSk', 'BWh', 'BWk', 'Cfa', 'Cfb','Cfc', 'Csa', 'Csb', 'Csc', 'Cwa','Cwb', 'Cwc', 'Dfa', 'Dfb', 'Dfc','Dfd', 'Dsa', 'Dsb', 'Dsc', 'Dsd','Dwa', 'Dwb', 'Dwc', 'Dwd', 'EF','ET', 'Ocean'] 
KGCcolors=["#960000", "#FF0000", "#FF6E6E", "#FFCCCC", "#CC8D14", "#CCAA54", "#FFCC00", "#FFFF64", "#007800", "#005000", "#003200", "#96FF00", "#00D700", "#00AA00", "#BEBE00", "#8C8C00", "#5A5A00", "#550055", "#820082", "#C800C8", "#FF6EFF", "#646464", "#8C8C8C", "#BEBEBE", "#E6E6E6", "#6E28B4", "#B464FA", "#C89BFA", "#C8C8FF", "#6496FF", "#64FFFF", "#F5FFFF"]

def brightness(rrggbb):
    """ W3C brightness definition
        input:
            hexadecimal color in the format:
            #RRGGBB
        output: value between 0 and 1
    """
    print(rrggbb)
    rr = int(rrggbb[1:3],16)/int('FF',16)
    gg = int(rrggbb[3:5],16)/int('FF',16)
    bb = int(rrggbb[5:7],16)/int('FF',16)
    #rr = math.floor(rrggbb/10000.)
    #gg = math.floor((rrggbb - rr*10000.)/100.)
    #bb = rrggbb - rr*10000 - gg*100
    return (rr * 299. + gg * 587. + bb * 114.) / 1000.

kgccolors = {}
for iKGC,KGCname in enumerate(KGCID):
    kgccolors[KGCname] = [KGCcolors[iKGC],'white' if (brightness(KGCcolors[iKGC])<0.5) else 'black']

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
kgcnames = {
    'Dfa':'snow \n fully humid \n hot summer',
    'Cfb':'green'       ,
    'BSk':'4'      ,
    'Csb':'5'      ,
    'Cfa':'darkgreen' ,     
    'BWh':'6'      ,
    'Aw' :'7'     ,
    'Dwc':'8'     ,
    'Dfb':'9'     ,
    #'Dfa':'',
}
for KGCID in list(pd.Series(xrkoeppen['KGCID'])):
    if KGCID not in kgcnames.keys():
        kgcnames[KGCID] = KGCID
    if KGCID not in kgccolors.keys():
        kgccolors[KGCID] = ['k','k']


koeppenlookuptable['color'] = ""
koeppenlookuptable['textcolor'] = ""
koeppenlookuptable['name'] = ""
for ikoeppen,koeppen in koeppenlookuptable.iterrows():
    print(ikoeppen)
    print(koeppen.KGCID)
    print(kgccolors[koeppen.KGCID])
    koeppenlookuptable['color'].loc[ikoeppen] = kgccolors[koeppen.KGCID][0]
    koeppenlookuptable['textcolor'].loc[ikoeppen] = kgccolors[koeppen.KGCID][1]
    koeppenlookuptable['name'].loc[ikoeppen] = kgcnames[koeppen.KGCID]



c4gldata[key].frames['stats']['records_all_stations_ini']['KGCname'] =  \
    c4gldata[key].frames['stats']['records_all_stations_ini']['KGC'].map(koeppenlookuptable['KGCID'])

print('sort the climate classes according to the amount ')
koeppenlookuptable['amount'] = ""

#exclude_koeppen = ['Dfc','Cwb']

for ikoeppen,koeppen in koeppenlookuptable.iterrows():

    print(ikoeppen,':',koeppen)
    kgc_select = (c4gldata[key].frames['stats']['records_all_stations_ini']['KGCname'] == koeppen['KGCID'])
    print(np.sum(kgc_select))
    koeppenlookuptable.iloc[ikoeppen]['amount'] = np.sum(kgc_select)

koeppenlookuptable = koeppenlookuptable.sort_values('amount',ascending=False)
include_koeppen = list(koeppenlookuptable.KGCID)


if args.make_figures:
    fig = plt.figure(figsize=(11,8.0))   #width,height
    i = 1                                                                           
    axes = {}         
    axes_taylor = {}         
    
    colors = ['r','g','b','m','y','purple','orange','sienna','navy']
    symbols = ['*','x','+']
    dias = {}


    i = 1
    for varkey in ['h','theta','q']:                                                    
        dias[varkey] =  TaylorDiagram(1., srange=[0.0,1.7],fig=fig, rect=(230+i+3),label='Reference')
        axes[varkey] = fig.add_subplot(2,3,i)                                       
        #axes_taylor[varkey] = fig.add_subplot(2,3,i+3)                                       
    
        #print(obs.std())
        dias[varkey]._ax.axis["left"].label.set_text(\
             "Normalized root mean square error")
        if i == 1:
            axes[varkey].annotate('Normalized standard deviation',\
                        xy= (0.05,0.37),
                        color='black',
                        rotation=90.,
                        xycoords='figure fraction',
                        weight='normal',
                        fontsize=10.,
                        horizontalalignment='center',
                        verticalalignment='center' ,
                        #bbox={'edgecolor':'black',
                        #      'boxstyle':'circle',
                        #      'fc':koeppen.color,
                        #      'alpha':1.0}
                       )


            # ticks = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x-1.))
            # #ax.axis["left"].axis.set_major_formatter(ticks)

            # # dias[varkey]._ax.axis["left"].axis.set_ticks(np.arange(0.,2.,0.25))
            # dias[varkey]._ax.axis["left"].axis.set_major_formatter(ticks)
        #dias[varkey]._ax.axis["left"].axis.set_ticks(np.arange(0.,2.,0.25))
        # Q95 = obs.quantile(0.95)
        # Q95 = obs.quantile(0.90)
        # Add RMS contours, and label them
        contours = dias[varkey].add_contours(levels=5, colors='0.5') # 5 levels
        dias[varkey].ax.clabel(contours, inline=1, fontsize=10, fmt='%.1f')
        #dia._ax.set_title(season.capitalize())
    
        dias[varkey].add_grid()
    
    
        #dia.ax.plot(x99,y99,color='k')
        i += 1
    
    i = 1
    for varkey in ['h','theta','q']:                                                    
        #for ikey,key in enumerate(args.experiments.strip(' ').split(' ')):
        for ikey,key in enumerate(args.experiments.strip(' ').split(' ')[:1]):
            # cc = c4gldata[key].frames['stats']['records_all_stations_ini']['cc']
            # clearsky = (cc < 0.05)
            # mod = c4gldata[key].frames['stats']['records_all_stations_mod_stats'].loc[clearsky]['d'+varkey+'dt']
            # obs = c4gldata[key].frames['stats']['records_all_stations_obs_afternoon_stats'].loc[clearsky]['d'+varkey+'dt']
            mod = c4gldata[key].frames['stats']['records_all_stations_mod_stats']['d'+varkey+'dt']
            obs = c4gldata[key].frames['stats']['records_all_stations_obs_afternoon_stats']['d'+varkey+'dt']

            print ('filtering classes that have sufficient samples: ', include_koeppen)
            filter_classes = (c4gldata[key].frames['stats']['records_all_stations_ini'].KGCname.isin(include_koeppen))
            mod = mod.loc[filter_classes]
            obs = obs.loc[filter_classes]
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
            dias[varkey].add_sample(STD/STD_OBS, PR,\
                           marker='o',ls='', mfc='white',mec='black',
                           zorder=-100,
                           ms=3.5*np.sqrt(np.sum(np.array(koeppenlookuptable.amount.values,dtype=np.float)))/\
                                np.mean(np.sqrt(np.array(koeppenlookuptable.amount.values,dtype=np.float)))
                           # annotate=koeppen.KGCID, color=koeppen.textcolor,weight='bold',fontsize=5.,\
                           # bbox={'edgecolor':'black','boxstyle':'circle','fc':koeppen.color,'alpha':0.7}
                           )
            dias[varkey].add_sample(STD/STD_OBS, PR,\
                           marker='o',ls='', mfc='none',mec='black',
                           zorder=700,
                           ms=3.5*np.sqrt(np.sum(np.array(koeppenlookuptable.amount.values,dtype=np.float)))/\
                                np.mean(np.sqrt(np.array(koeppenlookuptable.amount.values,dtype=np.float)))
                           # annotate=koeppen.KGCID, color=koeppen.textcolor,weight='bold',fontsize=5.,\
                           # bbox={'edgecolor':'black','boxstyle':'circle','fc':koeppen.color,'alpha':0.7}
                           )
            dias[varkey].add_sample(STD/STD_OBS, PR,\
                           marker='o',ls='', mfc='none',mec='black',
                           zorder=700,
                           ms=1.
                           # annotate=koeppen.KGCID, color=koeppen.textcolor,weight='bold',fontsize=5.,\
                           # bbox={'edgecolor':'black','boxstyle':'circle','fc':koeppen.color,'alpha':0.7}
                           )
            # dias[varkey].add_sample(STD/STD_OBS, PR,\
            #                    annotate='All', alpha=1.,color='black',weight='bold',fontsize=5.,\
            #                 zorder=100,\
            #                         bbox={'edgecolor':'black','boxstyle':'circle','alpha':0.0}\
            #                 )
    
        # put ticker position, see
        # https://matplotlib.org/examples/ticks_and_spines/tick-locators.html 
        # dia.ax.axis['bottom'].
        # dia.ax.axis['left'].
        # dia.ax.axis['left'].
    
        i += 1

    i = 1
    for varkey in ['h','theta','q']:                                                    

        for ikey,key in enumerate(args.experiments.strip(' ').split(' ')[:1]):
            icolor = 0
            for ikoeppen,koeppen in koeppenlookuptable.iterrows():
                if koeppen.amount >= 200:
                    print(ikoeppen,':',koeppen)
                    kgc_select = (c4gldata[key].frames['stats']['records_all_stations_ini']['KGCname'] == koeppen['KGCID'])
                    koeppen_mod = c4gldata[key].frames['stats']['records_all_stations_mod_stats']['d'+varkey+'dt'][kgc_select]
                    koeppen_obs = c4gldata[key].frames['stats']['records_all_stations_obs_afternoon_stats']['d'+varkey+'dt'][kgc_select]
    
                    #axes[varkey].scatter(koeppen_obs,koeppen_mod,marker=symbols[ikoeppen],color=colors[ikey])
                             #  label=key+", "+\
                             #                    'R = '+str(round(PR[0],3))+', '+\
                             #                    'RMSE = '+str(round(RMSE,5))+', '+\
                             #                    'BIAS = '+str(round(BIAS,5)),s=1.,color=colors[ikey])
    
    
    
                # # pl.scatter(obs,mod,label=key+", "+\
                # #                              'R = '+str(round(PR[0],3))+', '+\
                # #                              'RMSE = '+str(round(RMSE,5))+', '+\
                # #                              'BIAS = '+str(round(BIAS,5)),s=1.,color=colors[ikey])
                    
                    dias[varkey].add_sample(koeppen_mod.std()/koeppen_obs.std(),
                                   pearsonr(koeppen_mod,koeppen_obs)[0],
                                   marker='o',linewidth=0.5,
                                            mfc=koeppen.color,mec='black',#koeppen.color,
                                            zorder=300+icolor,
                                   ms=3.5*np.sqrt(koeppen.amount)/np.mean(np.sqrt(np.array(koeppenlookuptable.amount.values,dtype=np.float)))
                                   # annotate=koeppen.KGCID, color=koeppen.textcolor,weight='bold',fontsize=5.,\
                                   # bbox={'edgecolor':'black','boxstyle':'circle','fc':koeppen.color,'alpha':0.7}
                                   )
                    dias[varkey].add_sample(koeppen_mod.std()/koeppen_obs.std(),
                                   pearsonr(koeppen_mod,koeppen_obs)[0],
                                   marker='o',linewidth=0.5,
                                            mfc=koeppen.color,mec='black',#koeppen.color,
                                            zorder=301+icolor, ms=1
                                   # annotate=koeppen.KGCID, color=koeppen.textcolor,weight='bold',fontsize=5.,\
                                   # bbox={'edgecolor':'black','boxstyle':'circle','fc':koeppen.color,'alpha':0.7}
                                   )


                    # dias[varkey].add_sample(koeppen_mod.std()/koeppen_obs.std(),
                    #                pearsonr(koeppen_mod,koeppen_obs)[0],
                    #                         marker='o',linewidth=0.5, mfc='none',mec=str(koeppen.color),
                    #                         zorder=600+icolor,
                    #                ms=10.*np.sqrt(koeppen.amount)/np.mean(np.sqrt(np.array(koeppenlookuptable.amount.values,dtype=np.float)))
                    #                # annotate=koeppen.KGCID, color=koeppen.textcolor,weight='bold',fontsize=5.,\
                    #                # bbox={'edgecolor':'black','boxstyle':'circle','fc':koeppen.color,'alpha':0.7}
                    #                )

                    icolor += 1








    
            latex = {}
            latex['dthetadt'] =  r'$d \theta / dt $'
            latex['dqdt'] =      r'$d q / dt $'
            latex['dhdt'] =      r'$d h / dt $'
    
            axes[varkey].set_xlabel('Observed')     


            if varkey == 'q':
                units_final = r'[$g\, kg^{-1}\, h^{-1}$]'
            elif varkey == 'theta':
                units_final = r'[$K\, h^{-1}$]'
            elif varkey == 'h':
                units_final = r'[$m\, h^{-1}$]'

            axes[varkey].set_title(latex['d'+varkey+'dt']+' '+units_final,fontsize=12.)                                     
        if i==1:                                    
            axes[varkey].set_ylabel('Modelled')                                            
        abline(1,0,axis=axes[varkey])
        i +=1



    axph = fig.add_axes([0,0,1,1])
    axph.xaxis.set_visible(False)
    axph.yaxis.set_visible(False)
    axph.set_zorder(1000)
    axph.patch.set_alpha(0.0)
    
    axph.add_patch(Rectangle((0.006,0.01), 0.99, 0.13, alpha=0.71,
                           facecolor='white',edgecolor='grey'))


    idx = 0
    koeppenlookuptable_sel = koeppenlookuptable[koeppenlookuptable.amount >= 200].sort_values('KGCID',ascending=True)

    for ikoepp,koepp in koeppenlookuptable_sel.iterrows():
        xy_box = ((0.059+np.floor(idx/4)*0.335),0.127- 0.09*((np.mod(idx,4)/(len(koeppenlookuptable_sel)/4))))
        xy_text = list(xy_box)
        xy_text[0] += 0.009
        xy_text[1] -= 0.00
        xy_color = list(xy_box)
        xy_color[0] += -0.045
        xy_color[1] -= 0.019
        axph.add_patch(Rectangle(xy_color,0.049,0.023,
                                      edgecolor='black', 
                               #       boxstyle='round', 
                                      #xycoords='figure fraction',
                                      fc=koepp.color,
                                      alpha=1.0))
    
        axph.annotate(koepp.KGCID,
                    xy= xy_box,
                    color='white' if (brightness(koepp.color)<0.5) else 'black', 
                    family='monospace',
                    xycoords='figure fraction',
                    weight='bold',
                    fontsize=10.,
                    horizontalalignment='right',
                    verticalalignment='top' ,
                    # bbox={'edgecolor':'black',
                    #       'boxstyle':'round',
                    #       'fc':koepp.color,
                    #       'alpha':1.0}
                   )
        full_name = ""
        if koepp.KGCID is not "Ocean":
            for char in koepp.KGCID:
                full_name += lookup_symbols[char] + ' - '
            full_name = full_name[:-3]
        axph.annotate(full_name,xy=xy_text,color='k'
                    ,fontsize=8,xycoords='figure fraction',weight='bold',horizontalalignment='left',verticalalignment='top')

        idx +=1

    

    i = 1
    for varkey in ['h','theta','q']:                                                    
        ikey = 0
        key = list(args.experiments.strip().split(' '))[ikey]
        keylabel = keylabels[ikey]
        # cc = c4gldata[key].frames['stats']['records_all_stations_ini']['cc']
        # clearsky = (cc < 0.05)
    
        # mod = c4gldata[key].frames['stats']['records_all_stations_mod_stats'].loc[clearsky]['d'+varkey+'dt']
        # obs = c4gldata[key].frames['stats']['records_all_stations_obs_afternoon_stats'].loc[clearsky]['d'+varkey+'dt']
    
        mod = c4gldata[key].frames['stats']['records_all_stations_mod_stats']['d'+varkey+'dt']
        obs = c4gldata[key].frames['stats']['records_all_stations_obs_afternoon_stats']['d'+varkey+'dt']
        print ('filtering classes that have sufficient samples: ', include_koeppen)
        filter_classess = (c4gldata[key].frames['stats']['records_all_stations_ini'].KGCname.isin(include_koeppen))
        mod = mod.loc[filter_classes]
        obs = obs.loc[filter_classes]
    
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
        axes[varkey].contour(xi, yi, zi_int.reshape(xi.shape),levels=[0.25,0.5,0.75] ,
                colors=['darkred','lightgreen','darkred'],linewidths=[1,2,1])
        axes[varkey].contourf(xi, yi, zi_int.reshape(xi.shape),levels=[0.25,0.75] ,
                colors=['darkred'],alpha=0.5,)

        # if varkey == 'q':
        nani = np.concatenate([xi[zi != np.nan],yi[zi != np.nan]])
        # axes[varkey].set_ylim((np.percentile(nani,20),np.percentile(nani,80)))
        # #axes[varkey].set_ylim((nani.min(),nani.max()))
        # print(varkey,(nani.min(),nani.max()))
    
    
        latex = {}
        latex['dthetadt'] =  r'$d \theta / dt $'
        latex['dqdt'] =      r'$d q / dt $'
        latex['dhdt'] =      r'$d h / dt $'
    
        axes[varkey].set_xlabel('observations')     

        if varkey == 'q':
            units_final = r'[$\mathrm{g\, kg^{-1}\, h^{-1}}$]'
        elif varkey == 'theta':
            units_final = r'[$\mathrm{K\, h^{-1}}$]'
        elif varkey == 'h':
            units_final = r'[$\mathrm{m\, h^{-1}}$]'

        if varkey == 'q':
            axes[varkey].set_title(latex['d'+varkey+'dt']+' '+units_final,fontsize=15)        
        elif varkey == 'theta':
            axes[varkey].set_title(latex['d'+varkey+'dt']+' '+units_final,fontsize=15)
        elif varkey == 'h':
            axes[varkey].set_title(latex['d'+varkey+'dt']+' '+units_final,fontsize=15)
       #  c4gldata[key].frames['stats']['records_all_stations_ini']['KGCname']

        PR = pearsonr(mod,obs)[0]
        RMSE = rmse(obs,mod)                                               
        BIAS = np.mean(mod) - np.mean(obs)
        STD = mod.std()
    
        
        axes[varkey].scatter(obs,mod, label='All',s=0.1,alpha=0.14,color='k')



        #axes[varkey].legend(fontsize=5)

        #trans = ax.get_xaxis_transform() # x in data untis, y in axes fraction
        if varkey == 'q':
            annotate_text = \
                           'RMSE = '+format((RMSE*1000.),'0.2f')+r'$\,  \mathrm{g\,  kg^{-1}\,  h^{-1}}$'+ '\n'+\
                           'Bias = '+format((BIAS*1000.),'0.2f')+r'$\,  \mathrm{g\,  kg^{-1}\,  h^{-1}}$'+' \n'+\
                           r'$R$ = '+format(PR,'0.2f')
        elif varkey == 'h':
            annotate_text = \
                            'RMSE = '+format(RMSE,'0.1f')+r'$\,  \mathrm{m\, h^{-1}}$'+'\n'+\
                            'Bias = '+format(BIAS,'0.1f')+r'$\,  \mathrm{m\, h^{-1}}$'+'\n'+\
                            r'$R$ = '+format(PR,'0.2f')
        else:
            annotate_text = \
                            'RMSE = '+format(RMSE,'0.3f')+r'$\, \mathrm{K\, h^{-1}}$'+'\n'+\
                            'Bias = '+format(BIAS,'0.3f')+r'$\, \mathrm{K\, h^{-1}}$'+'\n'+\
                            r'$R$ = '+format(PR,'0.2f')


        ann = axes[varkey].annotate(annotate_text, xy=(0.05, .98 ), xycoords='axes fraction',fontsize=9,
       horizontalalignment='left', verticalalignment='top' 
        )

        if varkey == 'q':
            print('get_xlim not working well...STRANGE')
            limits =  [np.percentile(nani,1),np.percentile(nani,99)]
        else:
            limits =  [np.percentile(nani,1.0),np.percentile(nani,99.0)]

        axes[varkey].set_xlabel('Observed')     
        if i==1:                                    
            axes[varkey].set_ylabel('Modelled')                                            
        abline(1,0,axis=axes[varkey])


        i +=1

        #axes[varkey].axis('equal')
        axes[varkey].set_aspect('equal')
        axes[varkey].set_xlim(limits)
        axes[varkey].set_ylim(limits)

        # To specify the number of ticks on both or any single axes
        # plt.locator_params(axis='x', nbins=6)
        #plt.locator_params( nbins=10)
        axes[varkey].xaxis.set_major_locator(ticker.MaxNLocator(4))
        axes[varkey].yaxis.set_major_locator(ticker.MaxNLocator(4))
        # axes[varkey].xaxis.set_major_locator(ticker.MultipleLocator(5))
        # axes[varkey].yaxis.set_major_locator(ticker.MultipleLocator(5))

        if varkey == 'q':
            ticks = ticker.FuncFormatter(lambda x, pos:
                                         '{0:g}'.format(x*1000.))
            axes[varkey].xaxis.set_major_formatter(ticks)
            axes[varkey].yaxis.set_major_formatter(ticks)

        #     # axes[varkey].set_xticklabels(labels=ax.get_xticklabels()*1000.)
        #     # axes[varkey].set_yticklabels(labels=ax.get_yticklabels()*1000.)
    
    
    # # legend for different forcing simulations (colors)
    # ax = fig.add_axes([0.05,0.00,0.15,0.15]) #[*left*, *bottom*, *width*,    *height*]
    # leg = []
    # for ikey,key in enumerate(args.experiments.strip().split(' ')):
    #     leg1, = ax.plot([],colors[ikey]+'o' ,markersize=10)
    #     leg.append(leg1)
    # ax.axis('off')
    # #leg1 =
    # ax.legend(leg,list(args.experiments.strip().split(' ')),loc=2,fontsize=10)
    
    
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
    
    
    fig.subplots_adjust(top=0.95,bottom=0.20,left=0.08,right=0.94,hspace=0.35,wspace=0.29)
    
    
    #pl.legend(leglist,('EMI:WOC','EMI:MED','EMI:BEC'),loc=2,fontsize=16,prop={'family':
    # figfn = '/user/data/gent/gvo000/gvo00090/D2D/archive/report/global_eval_report_cs.png'
    # fig.savefig(figfn,dpi=200); print("Image file written to:", figfn)
    
    if args.figure_filename is not None:
        fig.savefig(args.figure_filename,dpi=200); print("Image file written to:",args.figure_filename)
        fig.savefig(args.figure_filename.replace('png','pdf')); print("Image file written to:", args.figure_filename)
    fig.show()  

    koeppenlookuptable = koeppenlookuptable.sort_index()
    if bool(args.show_control_parameters):


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



        #sns.set()
        #fig = pl.figure(figsize=(11,7))
        i = 1
        #axes = {}
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
        

        tempdatamodstats["source"] = "soundings"
        tempdatamodstats["source_index"] = "soundings"

        ini_ref = pd.DataFrame(c4gldata[key].frames['stats']['records_all_stations_ini'].copy())
        tempdataini_this = pd.DataFrame(ini_ref.copy())

        tempdatamodstats['dates']= tempdataini_this.ldatetime.dt.date
        tempdatamodstats['STNID']= tempdataini_this.STNID
        tempdatamodstats['source']= "soundings"
        tempdatamodstats['source_index']= "soundings"
        tempdatamodstats.set_index(['source_index','STNID','dates'],inplace=True)
        #print('hello')

        tempdataini = pd.DataFrame(ini_ref)
        tempdataini["source"] = "soundings"
        tempdataini["source_index"] = "soundings"
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
            #print('hello5

            # data[varkey] = tempdatamodstats['d'+varkey+'dt']
            data_all = pd.concat([data_all,tempdatamodstats],axis=0)
            data_input = pd.concat([data_input, tempdataini],axis=0)
            #print(data_input.shape)
            #print(data_all.shape)

        data_input.cc = data_input.cc.clip(0.,+np.inf)

        for varkey in ['h','theta','q']:
            varkey_full = 'd'+varkey+'dt ['+units[varkey]+'/h]'
            data_all = data_all.rename(columns={'d'+varkey+'dt':varkey_full})
            data_all['KGCname'] = data_input['KGCname']




            #print(data_input.shape)
            #print(data_all.shape)
        # xrkoeppen = xr.open_dataset('/user/data/gent/gvo000/gvo00090/EXT/data/KOEPPEN/Koeppen-Geiger.nc')
        # lookuptable = pd.Series(xrkoeppen['KGCID'])
        # data_all['KGCname'] = data_input['KGC'].map(lookuptable)
        #print('hello6')
        #print(data_all.columns)
        #print('hello7')

        


        varkeys = ['h','theta','q']
        for varkey in varkeys:
            #input_keys =['wg','cc']
            #for input_key in input_keys:
            varkey_full = 'd'+varkey+'dt ['+units[varkey]+'/h]'

            #print('hello8')
            #print(data_input.shape)
            #print(data_all.shape)
            #input_key_full = input_key + "["+units[input_key]+"]"
            #print('hello9')
            #print(data_input.shape)
            #print(data_all.shape)
            print ('Excluding extreme values from the classes plots')
            qvalmax = data_all[varkey_full].quantile(0.999)
            qvalmin = data_all[varkey_full].quantile(0.001)
            select_data = (data_all[varkey_full] >= qvalmin) & (data_all[varkey_full] < qvalmax)
            #print('hello11')
            data_all = data_all[select_data]
            #print('hello12')
            data_input = data_input[select_data.values]

            data_input = data_input[data_all.KGCname.isin(list(koeppenlookuptable.KGCID))]
            data_all = data_all[data_all.KGCname.isin(list(koeppenlookuptable.KGCID))]
            #print('hello13')
            #print(data_input.shape)
            #print(data_all.shape)
            #print('hello10')
            

        #sns.set(style="ticks", palette="deep")

        
        exppairs = {'obs|ref'     :['soundings',keylabels[0]],
                    'fcap|wilt'   :[keylabels[1] ,keylabels[2]],
                    'allveg|noveg':[keylabels[3],keylabels[4]]
                   }
        current_palette = sns.color_palette('deep')
        exppalettes = {'obs|ref'     :['white','grey'],
                    'fcap|wilt'   :[current_palette[0],current_palette[3]],
                       'allveg|noveg':[current_palette[2],current_palette[8]]
                   }

        data_all['expname'] = ""
        print('making alternative names for legends')
        expnames = {'soundings':'obs',\
                    keylabels[0]:'ref',\
                    keylabels[1]:'dry',\
                    keylabels[2]:'wet',\
                    keylabels[3]:'noveg',\
                    keylabels[4]:'fullveg',\
                   }
        for expname_orig,expname in expnames.items():
            data_all['expname'][data_all['source'] == expname_orig] = expname

        data_all['exppair'] = ""
        for exppairname,exppair in exppairs.items():
            data_all['exppair'][  data_all['source'].isin(exppair)  ] = exppairname

        icolor = 0
        
        koeppenlookuptable = koeppenlookuptable[koeppenlookuptable.amount >= 200]
        koeppenlookuptable = koeppenlookuptable.sort_values('amount',ascending=False)
        koeppenlookuptable = koeppenlookuptable[:9]
        fig, axes = plt.subplots(nrows=len(varkeys)*len(koeppenlookuptable), \
                                 ncols=len(exppairs), \
                                 figsize=(8, 13), #width, height
                                 sharex='col',
                                 #gridspec_kw=dict(height_ratios=(1, 3), 
                                 gridspec_kw=dict(hspace=0.20,wspace=0.08,top=0.94,bottom=0.06,left=0.15,right=0.99))

        data_all['d'+varkey+'dt ['+units[varkey]+'/h]'] *= 1000.  

        icol = 0
        irow = 0
        sns.set_style('whitegrid')
        for ikoeppen,koeppen in koeppenlookuptable.iterrows():
                for exppairname,exppair in exppairs.items():
                    for varkey in varkeys:
                        varkey_full = 'd'+varkey+'dt ['+units[varkey]+'/h]'
                        ax = axes[irow,icol]

                     #   axes[i] = fig.add_subplot(len(varkeys)*len(koeppenlookuptable),len(exppairs),icolor)
                #sns.violinplot(x='KGC',y=varkey_full,data=data_all,hue='source',linewidth=2.,palette="muted",split=True,inner='quart') #,label=key+", R = "+str(round(PR[0],3)),data=data)       
                
                #ax.set_title(input_key_full)
                        current_data = data_all[(data_all['exppair'] == exppairname) & (data_all['KGCname']  == koeppen.KGCID)]
                        sns.violinplot(y='exppair', x=varkey_full,
                                            hue="expname",split=True,
                                 palette=exppalettes[exppairname],
                                # palette=["m", "g",'r','b'],
                                 linewidth=1.0,inner='quart',
                                            data=current_data,sym='',legend=False,ax=ax)
                        ax.legend("")
                        ax.legend_.draw_frame(False)
                        ax.set_yticks([])
                        ax.set_ylabel("")

                        # if varkey == 'q':
                        #     ticks = ticker.FuncFormatter(lambda x, pos:
                        #                                  '{0:g}'.format(x*1000.))
                        #     ax.xaxis.set_major_formatter(ticks)

                        if varkey == 'q':
                            title_final = r'$dq/dt$'
                            xlabel_final = r'[$\mathrm{g\, kg^{-1}\, h^{-1}}$]'
                        elif varkey == 'theta':
                            title_final = r'$d\theta/dt$'
                            xlabel_final = r'[$\mathrm{K\, h^{-1}}$]'
                        elif varkey == 'h':
                            title_final = r'$dh/dt$'
                            xlabel_final = r'[$\mathrm{m\, h^{-1}}$]'


                        ax.set_xlabel("")
                        #sns.despine(left=True, right=True, bottom=False, top=False)
                        if irow == (len(varkeys)*len(koeppenlookuptable)-1):
                            #ax.set_frame_on(False)

                            ax.set_xlabel(xlabel_final)
                            ax.tick_params(top='off', bottom='on', left='off',
                                            right='off', labelleft='off',
                                            labeltop='off',
                                            labelbottom='on'
                                          )
                            ax.spines['top'].set_visible(False)
                            ax.spines['bottom'].set_visible(True)
                            ax.spines['left'].set_visible(True)
                            ax.spines['right'].set_visible(True)
                            #sns.despine(left=True, right=True, bottom=True, top=False)
                        elif irow == 0:
                            ax.set_title(title_final,fontsize=17.)
                            ax.tick_params(top='off', bottom='off', left='off',
                                            right='off', labelleft='off',
                                            labelbottom='off')
                            #ax.set_frame_on(False)
                            # ax.spines['left'].set_visible(True)
                            # ax.spines['right'].set_visible(True)
                            ax.spines['top'].set_visible(True)
                            ax.spines['bottom'].set_visible(False)
                            ax.spines['left'].set_visible(True)
                            ax.spines['right'].set_visible(True)
                            #sns.despine(left=True, right=True, bottom=False, top=True)
                            #ax.axis("off")
                        elif np.mod(irow,len(exppairs)) == 0:
                            ax.tick_params(top='off', bottom='off', left='off',
                                            right='off', labelleft='off',
                                            labelbottom='off')
                            #ax.set_frame_on(False)
                            # ax.spines['left'].set_visible(True)
                            # ax.spines['right'].set_visible(True)
                            ax.spines['top'].set_visible(True)
                            ax.spines['bottom'].set_visible(False)
                            ax.spines['left'].set_visible(True)
                            ax.spines['right'].set_visible(True)
                            #sns.despine(left=True, right=True, bottom=False, top=True)
                            #ax.axis("off")
                        elif np.mod(irow,len(exppairs)) == 2:
                            ax.tick_params(top='off', bottom='on', left='off',
                                            right='off', labelleft='off',
                                            labelbottom='off')
                            #ax.set_frame_on(False)
                            # ax.spines['left'].set_visible(True)
                            # ax.spines['right'].set_visible(True)
                            ax.spines['top'].set_visible(False)
                            ax.spines['bottom'].set_visible(True)
                            ax.spines['left'].set_visible(True)
                            ax.spines['right'].set_visible(True)
                            #sns.despine(left=True, right=True, bottom=False, top=True)
                            #ax.axis("off")
                        else:
                            ax.tick_params(top='off', bottom='off', left='off',
                                            right='off', labelleft='off',
                                            labelbottom='off')
                            #ax.set_frame_on(False)
                            #ax.spines['left'].set_visible(True)
                            #ax.spines['right'].set_visible(True)
                            ax.spines['top'].set_visible(False)
                            ax.spines['bottom'].set_visible(False)
                            ax.spines['left'].set_visible(True)
                            ax.spines['right'].set_visible(True)
                            #ax.axis("off")
                        icol +=1
                    irow +=1
                    icol=0

        idx = 0
        for ikoeppen,koeppen in koeppenlookuptable.iterrows():
            ax.annotate(koeppen.KGCID,
                        xy= (0.01,0.09 + (1.-0.12)*(1. - (idx+.5)/len(koeppenlookuptable))),
                        color=koeppen.textcolor, 
                        family='monospace',
                        xycoords='figure fraction',
                        weight='bold',
                        fontsize=8.,
                        horizontalalignment='top',
                        verticalalignment='left' ,
                        bbox={'edgecolor':'black',
                              'boxstyle':'circle',
                              'fc':koeppen.color,
                              'alpha':1.0}
                       )
            data_select = data_input[(data_input['source'] == 'soundings') & (data_input['KGCname']  == koeppen.KGCID)]
            ax.annotate('#:   '+r'$'+str(koeppen.amount)+' '+'$'+'\n'+\
                        'sm: '+r'$'+str(round(data_select.wg.mean(),2))+'$'+r'$\, \mathrm{m^{3}\, m^{-3}}$'+'\n'+\
                        'vf: '+str(round(data_select.cveg.mean(),2))+'\n'+\
                        'cc:  '+str(round(data_select.cc.mean(),2))+'\n'+\
                        'lat: '+r'$'+str(int(data_select.latitude.mean()))+r'\, ^\circ$'+' \n',\
                        xy= (0.01,0.015 + (1.-0.12)*(1. - (idx+.5)/len(koeppenlookuptable))),
                        color='black',
                        family='monospace',
                        xycoords='figure fraction',
                        weight='normal',
                        fontsize=8.,
                        horizontalalignment='top',
                        verticalalignment='left' ,
                        #bbox={'edgecolor':'black',
                        #      'boxstyle':'circle',
                        #      'fc':koeppen.color,
                        #      'alpha':1.0}
                       )
            idx+=1




            # if i ==1:
            #      plt.legend(loc='upper right',fontsize=7.,frameon=True,framealpha=0.7)
            # else:
            #      ax.get_legend().set_visible(False)
            # #     plt.legend('off')
            # if i >= 3:
            #     idx = 0
            #     for ikoeppen,koeppen in koeppenlookuptable.iterrows():

            #         ax.annotate(koeppen.KGCID,
            #                     xy=((idx+.5)/len(koeppenlookuptable),-0.00),
            #                     color=koeppen.textcolor, 
            #                     xycoords='axes fraction',
            #                     weight='bold',
            #                     fontsize=8.,
            #                     horizontalalignment='center',
            #                     verticalalignment='center' ,
            #                     bbox={'edgecolor':'black',
            #                           'boxstyle':'circle',
            #                           'fc':koeppen.color,
            #                           'alpha':1.0}
            #                    )
            #         idx+=1
            #     ax.set_xticklabels([])#labels=ax.get_xticklabels())
            #     ax.set_xlabel('Köppen climate class')
            # else:
            #     ax.set_xticklabels([])
            #     ax.set_xlabel('')

            # # ax.set_yticklabels([])
            # # ax.set_ylabel('')
            # if varkey == 'q':
            #     ticks = ticker.FuncFormatter(lambda x, pos:
            #                                  '{0:g}'.format(x*1000.))
            #     #ax.xaxis.set_major_formatter(ticks)
            #     ax.yaxis.set_major_formatter(ticks)

            #     ax.set_ylabel(latex['d'+varkey+'dt']+' ['+r'$10^{-3} \times $'+units['d'+varkey+'dt']+']')        
            # else:
            #     ax.set_ylabel(latex['d'+varkey+'dt']+' ['+units['d'+varkey+'dt']+']')        




            # for j,artist in enumerate(ax.artists):
            #     if np.mod(j,len(list(args.experiments.strip().split(' ')))+1) !=0:
            #         # Set the linecolor on the artist to the facecolor, and set the facecolor to None
            #         #print(j,artist)
            #         col = artist.get_facecolor()
            #         #print(j,artist)
            #         artist.set_edgecolor(col)
            #         #print(j,artist)
            #         artist.set_facecolor('None')
            # 
            #         # Each box has 6 associated Line2D objects (to make the whiskers, fliers, etc.)
            #         # Loop over them here, and use the same colour as above
            #         
            #         for k in range(j*5,j*5+5):
            #             line = ax.lines[k]
            #             line.set_color(col)
            #             line.set_mfc(col)
            #             line.set_mec(col)
            # 
            # # Also fix the legend
            # j = 0
            # for legpatch in ax.get_legend().get_patches():
            #     if j > 0:

            #         col = legpatch.get_facecolor()
            #         legpatch.set_edgecolor(col)
            #         legpatch.set_facecolor('None')
            #     j +=1





            #ax.grid()
            #sns.despine(offset=10, trim=True)
        fig.tight_layout()
        # fig.subplots_adjust(
        #     bottom=0.12,left=0.15,top=0.99,right=0.99,wspace=0.10,hspace=0.05,)
        if args.figure_filename_2 is not None:
            fig.savefig(args.figure_filename_2,dpi=200); print("Image file written to:", args.figure_filename_2)

            fig.savefig(args.figure_filename_2.replace('png','pdf')); print("Image file written to:", args.figure_filename_2)
        fig.show()



