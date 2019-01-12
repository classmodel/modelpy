#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  7 10:51:03 2017

@author: Hendrik Wouters

Purpose: provides class routines for ground and atmosphere conditions used for
the CLASS miced-layer model

Usage:
    from data_global import data_global
    from class4gl import class4gl_input
    from data_soundings import wyoming

    # create a data_global object and load initial data pages
    globaldata = data_global()
    globaldata.load_datasets()
    # create a class4gl_input object
    c4gli = class4gl_input()
    # Initialize it with profile data. We need to do this first. Actually this
    # will set the coordinate parameters (datetime, latitude, longitude) in
    # class4gl_input.pars.__dict__, which is required to read point data from
    # the data_global object.

    # open a Wyoming stream for a specific station
    wy_strm = wyoming(STNM=91376)
    # load the first profile
    wy_strm.find_first()
    # load the profile data into the class4gl_input object
    c4gli.get_profile_wyoming(wy_strm)
    
    # and finally, read the global input data for this profile
    c4gli.get_global_input(globaldata)


"""

import netCDF4 as nc4
import numpy as np
import datetime as dt
#you can install with
#import pynacolada as pcd
import pandas as pd
import xarray as xr
import os
import glob
import sys
import errno
import warnings
import logging


#formatter = logging.Formatter()
logging.basicConfig(format='%(asctime)s - \
                               %(name)s - \
                               %(levelname)s - \
                               %(message)s')

class book(object):
    """ this is a class for a dataset spread over multiple files. It has a
    similar purpose  open_mfdataset, but only 1 file (called current 'page')
    one is loaded at a time. This saves precious memory.  """
    def __init__(self,fn,concat_dim = None,debug_level=None):
        self.logger = logging.getLogger('book')
        if debug_level is not None:
            self.logger.setLevel(debug_level)

        # filenames are expanded as a list and sorted by filename
        self.pages = glob.glob(fn); self.pages.sort()
        # In case length of the resulting list is zero, this means no file was found that matches fn. In that case we raise an error.
        if len(self.pages) == 0:
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), fn)
        self.ipage = -1; self.page = None
        self.renames = {} # each time when opening a file, a renaming should be done.
        self.set_page(0)

        # we consider that the outer dimension is the one we concatenate
        self.concat_dim = concat_dim
        if self.concat_dim is None:
            self.concat_dim = self.concat_dim=list(self.page.dims.keys())[0]

    # this wraps the xarray sel-commmand
    def sel(*args, **kwargs):
        for dim in kwargs.keys():
            if dim == self.concat_dim:
                self.browse_page(**{dim: kwargs[dim]})
        return page.sel(*args,**kwargs)


    ## this wraps the xarray class -> some issues with that, so I just copy the sel command (which I do not use yet)
    #def __getattr__(self,attr):
    #    orig_attr = self.page.__getattribute__(attr)
    #    if callable(orig_attr):
    #        def hooked(*args, **kwargs):
    #            for dim in kwargs.keys():
    #                if dim == self.concat_dim:
    #                    self.browse_page(**{dim: kwargs[dim]})
    #
    #            result = orig_attr(*args, **kwargs)
    #            # prevent wrapped_class from becoming unwrapped
    #            if result == self.page:
    #                return self
    #            self.post()
    #            return result
    #        return hooked
    #    else:
    #        return orig_attr

    def set_renames(self,renames):
        #first, we convert back to original names, and afterwards, we apply the update of the renames.
        reverse_renames = dict((v,k) for k,v in self.renames.items())
        self.renames = renames
        if self.page is not None:
            self.page = self.page.rename(reverse_renames)
            self.page = self.page.rename(self.renames)

    def set_page(self,ipage,page=None):
        """ this sets the right page according to ipage:
                - We do not switch the page if we are already at the right one
                - we set the correct renamings (level -> lev, latitude -> lat,
                etc.)
                - The dataset is also squeezed.
        """

        if ((ipage != self.ipage) or (page is not None)):

            if self.page is not None:
                self.page.close()

            self.ipage = ipage
            if page is not None:
                self.page = page
            else:
                if self.ipage == -1:
                   self.page = None
                else:
                    #try:

                    self.logger.info("Switching to page "+str(self.ipage)+': '\
                                     +self.pages[self.ipage])
                    self.page = xr.open_dataset(self.pages[self.ipage])


            # do some final corrections to the dataset to make them uniform
            if self.page is not None:
               if 'latitude' in self.page.dims:
#    sel       f.library[fn] = self.library[fn].rename({'latitude':'lat','longitude':'lon'})

                   self.page = self.page.rename({'latitude':'lat','longitude':'lon'})
               if 'level' in self.page.dims:
                   self.page = self.page.rename({'level':'lev'})

               self.page = self.page.rename(self.renames)
               self.page = self.page.squeeze(drop=True)

    def browse_page(self,rewind=2,**args):

        # at the moment, this is only tested with files that are stacked according to the time dimension.
        dims = args.keys()


        if self.ipage == -1:
            self.set_page(0)

        found = False
        iipage = 0
        startipage = self.ipage - rewind
        while (iipage < len(self.pages)) and not found:
            ipage = (iipage+startipage) % len(self.pages)
            for dim in args.keys():
                this_file = True

                # here we store the datetimes in a directly-readable dictionary, so that we don't need to load it every time again
                if 'dims' not in self.__dict__:
                    self.dims = {}
                if dim not in self.dims.keys():
                    self.dims[dim] = [None]*len(self.pages)

                if self.dims[dim][ipage] is None:
                    self.logger.info('Loading coordinates of dimension "'+dim+\
                                     '" of page "' +str(ipage)+'".')
                    self.set_page(ipage)
                    # print(ipage)
                    # print(dim)
                    # print(dim,self.page[dim].values)
                    self.dims[dim][ipage] = self.page[dim].values

                # determine current time range of the current page
                mindim = self.dims[dim][ipage][0] -(self.dims[dim][ipage][1] - self.dims[dim][ipage][0])/2.
                maxdim = self.dims[dim][ipage][-1] +(self.dims[dim][ipage][-1] - self.dims[dim][ipage][-2])/2.

                if not ((args[dim] >= mindim) and (args[dim] < maxdim )):
                    this_file = False

            if this_file:
                found = True
                self.set_page(ipage)
            else:

                #if ((args[dim] >= self.page[dim].min().values) and (args[dim] < self.page[dim].max().values)):
                #    iipage = len(self.pages) # we stop searching

                iipage += 1

        if not found:
            self.logger.info("Page not found. Setting to page -1")
            #iipage = len(self.pages) # we stop searching further
            self.set_page(-1)

        if self.ipage != -1:
            self.logger.debug("I'm now at page "+ str(self.ipage)+': '+self.pages[self.ipage])
        else:
            self.logger.debug("I'm now at page "+ str(self.ipage))


class data_global(object):
    def __init__(self,sources= {
        'KOEPPEN:KGC'   : '/user/data/gent/gvo000/gvo00090/EXT/data/KOEPPEN/Koeppen-Geiger.nc',
        # # old gleam
        # 'GLEAM:wg'      : '/user/data/gent/gvo000/gvo00090/GLEAM/data/v3.1a/????/SMsurf_*_GLEAM_v3.1a.nc:SMsurf',
        # 'GLEAM:w2'      : '/user/data/gent/gvo000/gvo00090/GLEAM/data/v3.1a/????/SMroot_*_GLEAM_v3.1a.nc:SMroot',
        # 'GLEAM:BR'      : '/user/data/gent/gvo000/gvo00090/GLEAM/data/v3.1a/????/BR_*_GLEAM_v3.1a.nc:BR',
        # 'GLEAM:EF'      : '/user/data/gent/gvo000/gvo00090/GLEAM/data/v3.1a/????/EF_*_GLEAM_v3.1a.nc:EF',
        'GLEAM:wg'      : '/user/data/gent/gvo000/gvo00090/GLEAM/data/GLEAM_v3.2/v3.2a_OUTPUT/????/SMsurf_*_GLEAM_v3.2a.nc:SMsurf',
        'GLEAM:w2'      : '/user/data/gent/gvo000/gvo00090/GLEAM/data/GLEAM_v3.2/v3.2a_OUTPUT/????/SMroot_*_GLEAM_v3.2a.nc:SMroot',
        #'GLEAM:BR'      : '/user/data/gent/gvo000/gvo00090/GLEAM/data/GLEAM_v3.2/v3.2a/????/BR_*_GLEAM_v3.2a.nc:BR',
        'GLEAM:EF'      : '/user/data/gent/gvo000/gvo00090/GLEAM/data/GLEAM_v3.2/v3.2a_OUTPUT/????/EF_*_GLEAM_v3.2a.nc:EF',
        "IGBPDIS:alpha" : "/user/data/gent/gvo000/gvo00090/EXT/data/IGBP-DIS/FRACTIONS_GLEAMv31a.nc",
        "GLAS:z0m"      : "/user/data/gent/gvo000/gvo00090/EXT/data/GLAS/global_canopy_height_0.25.nc:Band1",
        "GLAS:z0h"      : "/user/data/gent/gvo000/gvo00090/EXT/data/GLAS/global_canopy_height_0.25.nc:Band1",
        'IGBPDIS:wsat'  : '/user/data/gent/gvo000/gvo00090/EXT/data/IGBP-DIS/wsat.nc',
        "ERAINT:Ts"  : "/user/data/gent/gvo000/gvo00090/EXT/data/ERA-INTERIM/by_var_nc/stl1_3hourly_xarray/stl1*_3hourly.nc:stl1",
        "ERAINT:Tsoil"  : "/user/data/gent/gvo000/gvo00090/EXT/data/ERA-INTERIM/by_var_nc/stl1_3hourly_xarray/stl1*_3hourly.nc:stl1",
        "ERAINT:T2"     : "/user/data/gent/gvo000/gvo00090/EXT/data/ERA-INTERIM/by_var_nc/stl2_3hourly_xarray/stl2*_3hourly.nc:stl2",
        "ERAINT:cc"     : "/user/data/gent/gvo000/gvo00090/EXT/data/ERA-INTERIM/by_var_nc/tcc_3hourly_xarray/tcc*_3hourly.nc:tcc",
        'IGBPDIS:wfc'   : '/user/data/gent/gvo000/gvo00090/EXT/data/IGBP-DIS/wfc.nc',
        'IGBPDIS:wwilt' : '/user/data/gent/gvo000/gvo00090/EXT/data/IGBP-DIS/wwp.nc:wwp',
        'MOD44B:cveg'   : '/user/data/gent/gvo000/gvo00090/EXT/data/MOD44B/fv.nc:fv',
        #'CERES:cc'      : '/user/data/gent/gvo000/gvo00090/EXT/data/CERES/CERES_SYN1deg-1H_Terra-Aqua-MODIS_Ed4A_Subset*.nc:cldarea_total_1h',
        "DSMW:b"        : "/user/data/gent/gvo000/gvo00090/EXT/data/DSMW/FAO_DSMW_DP.nc:DSMW:b",
        #"DSMW.C1sat"    : "/user/data/gent/gvo000/gvo00090/EXT/data/DSMW/FAO_DSMW_DP.nc:DSMW:C1sat",
        #"DSMW.C2ref"    : "/user/data/gent/gvo000/gvo00090/EXT/data/DSMW/FAO_DSMW_DP.nc:DSMW:C2ref",
        #"DSMW.p"        : "/user/data/gent/gvo000/gvo00090/EXT/data/DSMW/FAO_DSMW_DP.nc:DSMW:p",
        #"DSMW.a"        : "/user/data/gent/gvo000/gvo00090/EXT/data/DSMW/FAO_DSMW_DP.nc:DSMW:a",
        #"DSMW.CGsat"    : "/user/data/gent/gvo000/gvo00090/EXT/data/DSMW/FAO_DSMW_DP.nc:DSMW:CGsat",
        "GIMMS:LAIpixel": "/user/data/gent/gvo000/gvo00090/EXT/data/GIMMS/v2/LAI/gimms-3g.v2.lai.1981-2015_monmean_remapcon_0.25.nc:LAI",
        #'CERES.low': '/user/data/gent/gvo000/gvo00090/vsc42247/EXT/data/CERES/CERES_SYN1deg-1H_Terra-Aqua-MODIS_Ed4A_Subset_*.nc%cldarea_low_1h',
        #'CERES.cc%20000301%20100101': '/user/data/gent/gvo000/gvo00090/vsc42247/EXT/data/CERES/CERES_SYN1deg-1H_Terra-Aqua-MODIS_Ed4A_Subset_$YYYYMMDD_CERES_START-$YYYYMMDD_CERES_END.nc.cldarea_total_1h%cldarea_total_1h'
        "ERAINT:advt_x"     : "/user/data/gent/gvo000/gvo00090/EXT/data/ERA-INTERIM/by_var_nc/advt_x_6hourly/advt_x*_6hourly.nc:advt_x",
        "ERAINT:advt_y"     : "/user/data/gent/gvo000/gvo00090/EXT/data/ERA-INTERIM/by_var_nc/advt_y_6hourly/advt_y*_6hourly.nc:advt_y",
        "ERAINT:advq_x"     : "/user/data/gent/gvo000/gvo00090/EXT/data/ERA-INTERIM/by_var_nc/advq_x_6hourly/advq_x*_6hourly.nc",
        "ERAINT:advq_y"     : "/user/data/gent/gvo000/gvo00090/EXT/data/ERA-INTERIM/by_var_nc/advq_y_6hourly/advq_y*_6hourly.nc",
        "ERAINT:advu_x"     : "/user/data/gent/gvo000/gvo00090/EXT/data/ERA-INTERIM/by_var_nc/advu_x_6hourly/advu_x*_6hourly.nc",
        "ERAINT:advu_y"     : "/user/data/gent/gvo000/gvo00090/EXT/data/ERA-INTERIM/by_var_nc/advu_y_6hourly/advu_y*_6hourly.nc",
        "ERAINT:advv_x"     : "/user/data/gent/gvo000/gvo00090/EXT/data/ERA-INTERIM/by_var_nc/advv_x_6hourly/advv_x*_6hourly.nc",
        "ERAINT:advv_y"     : "/user/data/gent/gvo000/gvo00090/EXT/data/ERA-INTERIM/by_var_nc/advv_y_6hourly/advv_y*_6hourly.nc",
        #"ERAINT:divU_x"     : "/user/data/gent/gvo000/gvo00090/EXT/data/ERA-INTERIM/by_var_nc/divU_x_6hourly/divU_x*_6hourly.nc:__xarray_dataarray_variable__",
        #"ERAINT:divU_y"     : "/user/data/gent/gvo000/gvo00090/EXT/data/ERA-INTERIM/by_var_nc/divU_y_6hourly/divU_y*_6hourly.nc:__xarray_dataarray_variable__",
        "ERAINT:sp"     : "/user/data/gent/gvo000/gvo00090/EXT/data/ERA-INTERIM/by_var_nc/sp_6hourly/sp_*_6hourly.nc",
        "ERAINT:wp"  : '/user/data/gent/gvo000/gvo00090/EXT/data/ERA-INTERIM/by_var_nc/w_6hourly_xarray/w*_6hourly.nc:w',
        #"MSWEP:pr"    :"/user/data/gent/gvo000/gvo00090/EXT/data/MSWEP/MSWEP_v1.2_precip_1979-2015/3hr/raw_data/globe/*.nc:precipitation"
        },debug_level=None):
        self.library = {} #unique references to data sources being used. They can be files that are original on the disks or some unambiguous xarray virtual sources. These references are used in other variables. This way, a file or source cannot be loaded twice (a warning is made if one would try it).
        self.sources = sources
        self.datarefs = {}
        self.datasets = {}
        self.datetime = dt.datetime(1981,1,1)

        self.logger = logging.getLogger('data_global')
        if debug_level is not None:
            self.logger.setLevel(debug_level)
        self.debug_level = debug_level

        warnings.warn('omitting pressure field p and advection')

    def in_library(self,fn):
        if fn not in self.library.keys():
            return False
        else:
            print("Warning: "+fn+" is already in the library.")
            return True

    def add_to_library(self,fn):
        if not self.in_library(fn):
            print("opening: "+fn)
            self.library[fn] = \
                book(fn,concat_dim='time',debug_level=self.debug_level)

            #self.library[fn] = xr.open_mfdataset(fn,concat_dim='time')
            #if 'latitude' in self.library[fn].variables:
            #    self.library[fn] = self.library[fn].rename({'latitude':'lat','longitude':'lon'})


    # default procedure for loading datasets into the globaldata library
    def load_dataset_default(self,input_fn,varssource=None,varsdest=None):
        if type(varssource) is str:
            varssource = [varssource]
        if type(varsdest) is str:
            varsdest = [varsdest]

        self.add_to_library(input_fn)

        if varssource is None:
            varssource = []
            for var in self.sources[input_fn].variables:
                avoid = \
                ['lat','lon','latitude','longitude','time','lev','level']
                if ((len(list(var.shape)) >= 2) & (var not in avoid)): #two-dimensional array
                    varssource.append(var)

        if varsdest is None:
            varsdest = varssource

        #input_fn = "/user/data/gent/gvo000/gvo00090/EXT/data/IGBP-DIS/wsat.nc"
        for ivar,vardest in enumerate(varsdest):
            varsource = varssource[ivar]
            print('setting '+vardest+' as '+varsource+' from '+input_fn)

            if vardest in self.datarefs.keys():
                print("Warning! "+vardest+' is already provided by ',self.datarefs[vardest]+'. \n Overwriting....')
            #self.add_to_library(fn,varsource,vardest)
            if vardest != varsource:
                libkey = input_fn+'.'+varsource+'.'+vardest
                if libkey not in self.library.keys():
                    #self.library[libkey] = self.library[input_fn].rename({varsource:vardest})
                    self.library[libkey] = book(input_fn,\
                                                debug_level=self.debug_level)
                    self.library[libkey].set_renames({varsource: vardest})

                self.datarefs[vardest] = libkey # this is to remember that it was originally varsource in input_fn
                self.datasets[vardest] =self.library[self.datarefs[vardest]]
            else:
                self.datarefs[vardest] = input_fn
                self.datasets[vardest] =self.library[self.datarefs[vardest]]

            # if ((vardest is not None) & (vardest not in self.datasets[vardest].variables)):
            #     print('Warning: '+ vardest "not in " + input_fn)



    def load_datasets(self,sources = None,recalc=0):

        if sources is None:
            sources = self.sources
        for key in sources.keys():
            #datakey,vardest,*args = key.split(':')
            datakey,vardest = key.split(':')
            #print(datakey)

            fnvarsource = sources[key].split(':')
            if len(fnvarsource) > 2:
                #fn,varsource,*fnargs = fnvarsource
                fn,varsource,fnargs = fnvarsource
                fnargs = [fnargs]
            elif len(fnvarsource) > 1:
                #fn,varsource,*fnargs = fnvarsource
                fn,varsource = fnvarsource
                fnargs = []
            else:
                fn = sources[key]
                varsource = vardest
            self.load_dataset(fn,varsource,vardest,datakey,recalc=recalc)

    def load_dataset(self,fn,varsource,vardest,datakey,recalc=0):
            # the default way of loading a 2d dataset
            if datakey in ['CERES','GLEAM','ERAINT','GIMMS']:
                self.load_dataset_default(fn,varsource,vardest)
            elif datakey == 'IGBPDIS':
                if vardest == 'alpha':
                    ltypes = ['W','B','H','TC']
                    for ltype in ltypes:
                        self.load_dataset_default(fn,'f'+ltype,'f'+ltype)
                        ##self.datasets['f'+ltype]['f'+ltype]=  self.datasets['f'+ltype]['f'+ltype].squeeze(drop=True)


                    # landfr = {}
                    # for ltype in ['W','B','H','TC']:
                    #     landfr[ltype] = datasets['f'+ltype]['f'+ltype].values



                    keytemp = 'alpha'
                    fnkeytemp = fn+':IGBPDIS:alpha'
                    if (os.path.isfile(fnkeytemp)) and ( recalc < 6):
                        self.library[fnkeytemp]  = book(fnkeytemp,
                                                        debug_level=self.debug_level)
                        self.datasets[keytemp] = self.library[fnkeytemp]
                        self.datarefs[keytemp] = fnkeytemp
                    else:
                        self.library[fn+':IGBPDIS:alpha'] = xr.Dataset()
                        #self.library[fn+':IGBPDIS:alpha'][keytemp] = xr.zeros_like(self.datasets['IGBPDIS']['IGBPDIS'],dtype=np.float)*np.nan
                        self.library[fn+':IGBPDIS:alpha']['lat'] = self.datasets['fW'].page['lat']
                        self.library[fn+':IGBPDIS:alpha']['lon'] = self.datasets['fW'].page['lon']
                        self.library[fn+':IGBPDIS:alpha'][keytemp] = xr.DataArray(np.zeros(shape=(self.datasets['fW'].page['lon'].shape[0],self.datasets['fW'].page['lat'].shape[0]),dtype=np.float),dims=('lon','lat'))
                        self.datasets[keytemp] = self.library[fn+':IGBPDIS:alpha']
                        self.datarefs[keytemp] =fn+':IGBPDIS:alpha'

                        aweights = {'W':0.075,'TC':0.15,'H':0.22,'B':0.30}

                        alpha=self.library[fn+':IGBPDIS:alpha'][keytemp].values
                        for ltype in ltypes:
                            alpha += self.datasets['f'+ltype].page['f'+ltype].values*aweights[ltype]

                        self.library[fn+':IGBPDIS:alpha'][keytemp].values = alpha
                        print('writing file to: '+fnkeytemp)
                        os.system('rm '+fnkeytemp)
                        self.library[fnkeytemp].to_netcdf(fnkeytemp)
                        self.library[fnkeytemp].close()


                        self.library[fnkeytemp]  = \
                            book(fnkeytemp,debug_level=self.debug_level)
                        self.datasets[keytemp] = self.library[fnkeytemp]
                        self.datarefs[keytemp] = fnkeytemp


                else:
                    self.load_dataset_default(fn,varsource,vardest)


            elif datakey == 'GLAS':
                self.load_dataset_default(fn,varsource,vardest)
                if vardest == 'z0m':
                    self.datasets['z0m'].page['z0m'].values = (self.datasets['z0m'].page['z0m'].values/10.).clip(0.01,None)
                elif vardest == 'z0h':
                    self.datasets['z0h'].page['z0h'].values = (self.datasets['z0h'].page['z0h'].values/100.).clip(0.001,None)
            elif datakey == 'DSMW':


                # Procedure of the thermal properties:
                # 1. determine soil texture from DSMW/10.
                # 2. soil type with look-up table (according to DWD/EXTPAR)
                # 3. Thermal properties used in the force-restore method (Clapp and Hornberger, 1987)
                #    with parameter look-up table from Noilhan and Planton (1989).
                #    Note: The look-up table is inspired on DWD/COSMO

                # to do: implement inheretance, so that the the preliminary output of DSMW or any other dataset can be calculated first



                fnout = fn.replace('*','') # for storing computationally heavy soil properties, instead of calculating everytime
                self.load_dataset_default(fn,'DSMW')
                print('calculating texture')
                SPKEYS = ['tex_coarse', 'tex_medium', 'tex_fine', 'code','undefined']
                TEMP  = {}
                TEMP2 = self.datasets['DSMW'].page['DSMW'].values
                TEMP3 = {}
                for SPKEY in SPKEYS:


                    keytemp = SPKEY+'_values'
                    fnoutkeytemp = fnout+':DSMW:'+keytemp
                    if (os.path.isfile(fnoutkeytemp)) and ( recalc < 5 ):
                        self.library[fn+':DSMW:'+SPKEY+'_values'] = \
                                book(fnoutkeytemp,debug_level=self.debug_level)
                        self.datasets[SPKEY+'_values'] = self.library[fn+':DSMW:'+SPKEY+'_values']
                        self.datarefs[SPKEY+'_values'] =fn+':DSMW:'+SPKEY+'_values'


                    else:
                        #DSMW = self.datasets['DSMW']['DSMW']#   self.input_nc.variables['DSMW'][ilat,ilon]
                        self.library[fn+':DSMW:'+SPKEY+'_values'] = xr.Dataset()
                        self.library[fn+':DSMW:'+SPKEY+'_values']['lat'] = self.datasets['DSMW'].page['lat']
                        self.library[fn+':DSMW:'+SPKEY+'_values']['lon'] = self.datasets['DSMW'].page['lon']
                        self.library[fn+':DSMW:'+SPKEY+'_values'][SPKEY+'_values'] = xr.DataArray(np.zeros(shape=(self.datasets['DSMW'].page['lat'].shape[0],self.datasets['DSMW'].page['lon'].shape[0]),dtype=np.int),dims=('lat','lon'))
                        #self.library[fn+':DSMW:'+SPKEY+'_values'][SPKEY+'_values'] = xr.zeros_like(self.datasets['DSMW']['DSMW'],dtype=(np.int if SPKEY == 'code' else np.float))
                        self.datasets[SPKEY+'_values'] = self.library[fn+':DSMW:'+SPKEY+'_values']
                        self.datarefs[SPKEY+'_values'] =fn+':DSMW:'+SPKEY+'_values'

                        # for faster computation, we need to get it to memory out of Dask.
                        TEMP[SPKEY] = self.datasets[SPKEY+'_values'][SPKEY+'_values'].values
                        TEMP3[SPKEY] = self.datasets['DSMW'].page[SPKEY].values

                # yes, I know I only check the last file.
                if not ((os.path.isfile(fnoutkeytemp)) and ( recalc < 5)):
                    for idx in range(len(self.datasets['DSMW'].page['tex_coarse'].values))[:]:
                        print('idx',idx,SPKEY)
                        SEL = (TEMP2 == idx)
                    #     print(idx,len(TEMP3))
                        for SPKEY in SPKEYS:
                            TEMP[SPKEY][SEL] = TEMP3[SPKEY][idx]

                    for SPKEY in SPKEYS:
                        keytemp = SPKEY+'_values'
                        fnoutkeytemp = fnout+':DSMW:'+keytemp
                        self.datasets[SPKEY+'_values'][SPKEY+'_values'].values = TEMP[SPKEY][:]
                        os.system('rm '+fnoutkeytemp)
                        self.datasets[SPKEY+'_values'].to_netcdf(fnoutkeytemp)
                        self.datasets[SPKEY+'_values'].close()


                        self.library[fn+':DSMW:'+SPKEY+'_values'] = \
                                book(fnoutkeytemp,debug_level=self.debug_level)
                        self.datasets[SPKEY+'_values'] = self.library[fn+':DSMW:'+SPKEY+'_values']
                        self.datarefs[SPKEY+'_values'] =fn+':DSMW:'+SPKEY+'_values'


                keytemp = 'texture'
                fnoutkeytemp=fnout+':DSMW:'+keytemp
                if (os.path.isfile(fnoutkeytemp)) and ( recalc < 3 ):
                    self.library[fnoutkeytemp]  = \
                        book(fnoutkeytemp,debug_level=self.debug_level)
                    self.datasets[keytemp] = self.library[fn+':DSMW:texture']
                    self.datarefs[keytemp] =fn+':DSMW:texture'
                else:
                    self.library[fn+':DSMW:texture'] = xr.Dataset()
                    #self.library[fn+':DSMW:texture'][keytemp] = xr.zeros_like(self.datasets['DSMW']['DSMW'],dtype=np.float)*np.nan
                    self.library[fn+':DSMW:texture']['lat'] = self.datasets['DSMW'].page['lat']
                    self.library[fn+':DSMW:texture']['lon'] = self.datasets['DSMW'].page['lon']
                    self.library[fn+':DSMW:texture'][keytemp] = xr.DataArray(np.zeros(shape=(self.datasets['DSMW'].page['lat'].shape[0],self.datasets['DSMW'].page['lon'].shape[0]),dtype=np.float),dims=('lat','lon'))
                    self.datasets[keytemp] = self.library[fn+':DSMW:texture']
                    self.datarefs[keytemp] =fn+':DSMW:texture'



                    self.datasets[keytemp][keytemp].values = (0.5*self.datasets['tex_medium_values'].page['tex_medium_values'].values+1.0*self.datasets['tex_coarse_values'].page['tex_coarse_values'].values)/(self.datasets['tex_coarse_values'].page['tex_coarse_values'].values+self.datasets['tex_medium_values'].page['tex_medium_values'].values+self.datasets['tex_fine_values'].page['tex_fine_values'].values)

                    zundef = np.array(self.datasets['undefined_values'].page['undefined_values'].values,dtype=np.float)
                    zundef[zundef < 0] = np.nan
                    zsum_tex = self.datasets['tex_coarse_values'].page['tex_coarse_values'].values+self.datasets['tex_medium_values'].page['tex_medium_values'].values+ self.datasets['tex_fine_values'].page['tex_fine_values'].values
                    VALID  = (zsum_tex >= zundef) *( ~np.isnan(zundef))

                    self.datasets[keytemp][keytemp].values[~VALID] = 9012.

                    os.system('rm '+fnoutkeytemp)
                    self.datasets[keytemp].to_netcdf(fnoutkeytemp)
                    self.datasets[keytemp].close()


                    self.library[fnoutkeytemp]  = \
                        book(fnoutkeytemp,debug_level=self.debug_level)
                    self.datasets[keytemp] = self.library[fn+':DSMW:texture']
                    self.datarefs[keytemp] =fn+':DSMW:texture'


                print('calculating texture type')



                keytemp = 'itex'
                fnoutkeytemp=fnout+':DSMW:'+keytemp
                if (os.path.isfile(fnoutkeytemp)) and ( recalc < 2 ):
                    self.library[fnoutkeytemp] = \
                            book(fnoutkeytemp,debug_level=self.debug_level)
                    self.datasets[keytemp] = self.library[fn+':DSMW:itex']
                    self.datarefs[keytemp] =fn+':DSMW:itex'
                else:
                    self.library[fnoutkeytemp] = xr.Dataset()
                    self.library[fnoutkeytemp][keytemp] = xr.zeros_like(self.datasets['DSMW'].page['DSMW'],dtype=np.int)
                    self.datasets[keytemp] = self.library[fn+':DSMW:itex']
                    self.datarefs[keytemp] =fn+':DSMW:itex'

                    X = self.datasets['texture'].page['texture'].values*100
                    X[pd.isnull(X)] = -9


                    self.datasets[keytemp][keytemp].values = X

                    os.system('rm '+fnoutkeytemp)
                    self.datasets['itex'].to_netcdf(fnoutkeytemp)
                    self.datasets['itex'].close()


                    self.library[fnoutkeytemp] = \
                            book(fnoutkeytemp,debug_level=self.debug_level)
                    self.datasets[keytemp] = self.library[fn+':DSMW:itex']
                    self.datarefs[keytemp] =fn+':DSMW:itex'


                keytemp = 'isoil'
                fnoutkeytemp=fnout+':DSMW:'+keytemp
                isoil_reprocessed = False
                if (os.path.isfile(fnoutkeytemp)) and ( recalc < 1):
                    self.library[fn+':DSMW:isoil'] = \
                            book(fnoutkeytemp,debug_level=self.debug_level)
                    self.datasets['isoil'] = self.library[fn+':DSMW:isoil']
                    self.datarefs['isoil'] =fn+':DSMW:isoil'
                else:
                    isoil_reprocessed = True
                    print('calculating soil type')
                    self.library[fn+':DSMW:isoil'] = xr.Dataset()
                    self.library[fn+':DSMW:isoil']['isoil'] = xr.zeros_like(self.datasets['DSMW'].page['DSMW'],dtype=np.int)
                    self.datasets['isoil'] = self.library[fn+':DSMW:isoil']
                    self.datarefs['isoil'] =fn+':DSMW:isoil'

                    #adopted from mo_agg_soil.f90 (EXTPAR3.0)
                    self.datasets['isoil']['isoil'] = xr.zeros_like(self.datasets['DSMW'].page['DSMW'],dtype=np.int)
                    ITEX = self.datasets['itex'].page['itex'].values
                    ISOIL = 9 + 0.*self.datasets['isoil']['isoil'].values
                    LOOKUP = [
                              [-10 ,9],# ocean
                              [0 ,7],# fine textured, clay (soil type 7)
                              [20,6],# medium to fine textured, loamy clay (soil type 6)
                              [40,5],# medium textured, loam (soil type 5)
                              [60,4],# coarse to medium textured, sandy loam (soil type 4)
                              [80,3],# coarse textured, sand (soil type 3)
                              [100,9],# coarse textured, sand (soil type 3)
                            ]
                    for iitex,iisoil in LOOKUP:
                        ISOIL[ITEX > iitex] = iisoil
                        print('iitex,iisoil',iitex,iisoil)


                    #adopted from mo_agg_soil.f90 (EXTPAR3.0)
                    LOOKUP = [
                              [9001, 1 ], # ice, glacier (soil type 1)
                              [9002, 2 ], # rock, lithosols (soil type 2)
                              [9003, 3 ], # salt, set soiltype to sand (soil type 3)
                              [9004, 8 ], # histosol, e.g. peat (soil type 8)
                              [9,    9 ], # undefined (ocean)
                              [9005, 3 ], # shifting sands or dunes, set soiltype to sand (soil type 3)
                              [9000, 9 ], # undefined (inland lake)
                              [9009, 5 ], #  default_soiltype ! undefined (nodata), set soiltype to loam (soil type )
                              [9012, 5 ], #  default_soiltype undefined (dominant part undefined), set soiltype to loam (soil type 5)
                            ]
                    # EXTPAR: soil_code = soil_texslo(soil_unit)%dsmw_code # the legend has some special cases for the "soil_code"
                    CODE_VALUES = self.datasets['code_values'].page['code_values'].values

                    CODE_VALUES[ITEX == 901200] = 9012
                    for icode,iisoil in LOOKUP:
                        ISOIL[CODE_VALUES == icode] = iisoil

                    self.datasets['isoil']['isoil'].values = ISOIL
                    os.system('rm '+fnoutkeytemp)
                    self.datasets[keytemp].to_netcdf(fnoutkeytemp)
                    self.datasets[keytemp].close()
                    print('saved inbetween file to: '+fnoutkeytemp)

                    self.library[fn+':DSMW:isoil'] = \
                            book(fnoutkeytemp,debug_level=self.debug_level)
                    self.datasets['isoil'] = self.library[fn+':DSMW:isoil']
                    self.datarefs['isoil'] =fn+':DSMW:isoil'

                #adopted from data_soil.f90 (COSMO5.0)
                SP_LOOKUP = {
                  # soil type:         ice        rock       sand        sandy      loam         clay        clay        peat        sea        sea
                  # (by index)                                           loam                    loam                                water      ice
                  'cporv'  : [ np.nan, 1.E-10   , 1.E-10   , 0.364     , 0.445     , 0.455     , 0.475     , 0.507     , 0.863     , 1.E-10   , 1.E-10   ],
                  'cfcap'  : [ np.nan, 1.E-10   , 1.E-10   , 0.196     , 0.260     , 0.340     , 0.370     , 0.463     , 0.763     , 1.E-10   , 1.E-10   ],
                  'cpwp'   : [ np.nan, 0.0      , 0.0      , 0.042     , 0.100     , 0.110     , 0.185     , 0.257     , 0.265     , 0.0      ,  0.0     ],
                  'cadp'   : [ np.nan, 0.0      , 0.0      , 0.012     , 0.030     , 0.035     , 0.060     , 0.065     , 0.098     , 0.0      ,  0.0     ],
                  'crhoc'  : [ np.nan, 1.92E6   , 2.10E6   , 1.28E6    , 1.35E6    , 1.42E6    , 1.50E6    , 1.63E6    , 0.58E6    , 4.18E6   , 1.92E6   ],
                  'cik2'   : [ np.nan, 0.0      , 0.0      , 0.0035    , 0.0023    , 0.0010    , 0.0006    , 0.0001    , 0.0002    , 0.0      ,  0.0     ],
                  'ckw0'   : [ np.nan, 0.0      , 0.0      , 479.E-7   , 943.E-8   , 531.E-8   , 764.E-9   , 17.E-9    , 58.E-9    , 0.0      ,  0.0     ],
                  'ckw1'   : [ np.nan, 0.0      , 0.0      , -19.27    , -20.86    , -19.66    , -18.52    , -16.32    , -16.48    , 0.0      ,  0.0     ],
                  'cdw0'   : [ np.nan, 0.0      , 0.0      , 184.E-7   , 346.E-8   , 357.E-8   , 118.E-8   , 442.E-9   , 106.E-9   , 0.0      ,  0.0     ],
                  'cdw1'   : [ np.nan, 0.0      , 0.0      , -8.45     , -9.47     , -7.44     , -7.76     , -6.74     , -5.97     , 0.0      ,  0.0     ],
                  'crock'  : [ np.nan, 0.0      , 0.0      , 1.0       , 1.0       , 1.0       , 1.0       , 1.0       , 1.0       , 0.0      ,  0.0     ],
                  'cala0'  : [ np.nan, 2.26     , 2.41     , 0.30      , 0.28      , 0.25      , 0.21      , 0.18      , 0.06      , 1.0      ,  2.26    ],
                  'cala1'  : [ np.nan, 2.26     , 2.41     , 2.40      , 2.40      , 1.58      , 1.55      , 1.50      , 0.50      , 1.0      ,  2.26    ],
                  'csalb'  : [ np.nan, 0.70     , 0.30     , 0.30      , 0.25      , 0.25      , 0.25      , 0.25      , 0.20      , 0.07     ,  0.70    ],
                  'csalbw' : [ np.nan, 0.00     , 0.00     , 0.44      , 0.27      , 0.24      , 0.23      , 0.22      , 0.10      , 0.00     ,  0.00    ],
                  'ck0di'  : [ np.nan, 1.E-4    , 1.E-4    , 2.E-4     , 2.E-5     , 6.E-6     , 2.E-6     , 1.E-6     , 1.5E-6    , 0.00     ,  0.00    ],
                  'cbedi'  : [ np.nan, 1.00     , 1.00     , 3.5       , 4.8       , 6.1       , 8.6       , 10.0      , 9.0       , 0.00     ,  0.00    ],
                  'csandf' : [ np.nan, 0.0      , 0.0      , 90.       , 65.       , 40.       , 35.       , 15.       , 90.       , 0.00     ,  0.00    ],
                  'cclayf' : [ np.nan, 0.0      , 0.0      , 5.0       , 10.       , 20.       , 35.       , 70.       , 5.0       , 0.00     ,  0.00    ],
                  # Important note: For peat, the unknown values below are set equal to that of loam
                  #supplement Noihhan andf Planton 1989 soil texture parameters for the force-restore method.
                  'b'      : [ np.nan, np.nan   , np.nan   , 4.05      , 4.90      , 5.39      , 8.52      , 11.40     , 5.39    , np.nan   ,  np.nan  ],
                  #error in table 2 of NP89: values need to be multiplied by e-6
                  'CGsat'  : [ np.nan, np.nan   , np.nan   , 3.222e-6     , 3.560e-6     , 4.111e-6     , 3.995e-6     , 3.600e-6     , np.nan    , np.nan   ,  np.nan  ],
                  'p'  :     [ np.nan, np.nan   , np.nan   , 4.        , 4.        , 6.        , 10.       , 12.       , 6.    , np.nan   ,  np.nan  ],

                  'a'  :     [ np.nan, np.nan   , np.nan   , 0.387     , 0.219     , 0.148     , 0.084     , 0.083     , 0.148    , np.nan   ,  np.nan  ],
                  'C1sat'  : [ np.nan, np.nan   , np.nan   , 0.082     , 0.132     , 0.191     , 0.227     , 0.342     , 0.191    , np.nan   ,  np.nan  ],
                  'C2ref'  : [ np.nan, np.nan   , np.nan   , 3.9       , 1.8       , 0.8       , 0.6       , 0.3       , 0.8    , np.nan   ,  np.nan  ],
                }


                # isoil_reprocessed = False
                # if (os.path.isfile(fnoutkeytemp)) and ( recalc < 1):

                #     self.library[fn+':DSMW:isoil'] = \
                #             book(fnoutkeytemp,debug_level=self.debug_level)
                #     self.datasets['isoil'] = self.library[fn+':DSMW:isoil']
                #     self.datarefs['isoil'] =fn+':DSMW:isoil'
                # else:
                #     isoil_reprocessed = True
                #     print('calculating soil type')
                #     self.library[fn+':DSMW:isoil'] = xr.Dataset()
                #     self.library[fn+':DSMW:isoil']['isoil'] = xr.zeros_like(self.datasets['DSMW'].page['DSMW'],dtype=np.int)
                #     self.datasets['isoil'] = self.library[fn+':DSMW:isoil']
                #     self.datarefs['isoil'] =fn+':DSMW:isoil'




                # this should become cleaner in future but let's hard code it for now.
                DSMWVARS = ["b", "C1sat","C2ref","p","a" ]
                print('calculating soil parameter')
                DATATEMPSPKEY = {}
                if (recalc < 1) and (isoil_reprocessed == False): 
                    for SPKEY in DSMWVARS:#SP_LOOKUP.keys():
                        keytemp = SPKEY
                        fnoutkeytemp=fnout+':DSMW:'+keytemp
                        self.library[fn+':DSMW:'+SPKEY] =\
                                book(fnoutkeytemp,debug_level=self.debug_level)
                        self.datasets[SPKEY] = self.library[fnoutkeytemp]
                        self.datarefs[SPKEY] =fnoutkeytemp
                else:
                    for SPKEY in DSMWVARS:#SP_LOOKUP.keys():

                        self.library[fn+':DSMW:'+SPKEY] = xr.Dataset()
                        self.library[fn+':DSMW:'+SPKEY][SPKEY] = xr.zeros_like(self.datasets['DSMW'].page['DSMW'],dtype=np.float)
                        self.datasets[SPKEY] = self.library[fn+':DSMW:'+SPKEY]
                        self.datarefs[SPKEY] =fn+':DSMW:'+SPKEY
                        DATATEMPSPKEY[SPKEY] = self.datasets[SPKEY][SPKEY].values
                    ISOIL = self.datasets['isoil'].page['isoil'].values
                    print(np.where(ISOIL>0.))
                    for i in range(11):
                        SELECT = (ISOIL == i)
                        for SPKEY in DSMWVARS:#SP_LOOKUP.keys():
                            DATATEMPSPKEY[SPKEY][SELECT] = SP_LOOKUP[SPKEY][i]

                    for SPKEY in DSMWVARS:#SP_LOOKUP.keys():
                        self.datasets[SPKEY][SPKEY].values = DATATEMPSPKEY[SPKEY]

                        os.system('rm '+fn+':DSMW:'+SPKEY)
                        self.datasets[SPKEY].to_netcdf(fn+':DSMW:'+SPKEY)
                        self.datasets[SPKEY].close()
                        print('saved inbetween file to: '+fn+':DSMW:'+SPKEY)

                        self.library[fn+':DSMW:'+SPKEY] = \
                                book(fn+':DSMW:'+SPKEY,debug_level=self.debug_level)
                        self.datasets[SPKEY] = self.library[fn+':DSMW:'+SPKEY]
                        self.datarefs[SPKEY] =fn+':DSMW:'+SPKEY


            else:
                self.load_dataset_default(fn,varsource,vardest)






#
#                 # only print the last parameter value in the plot
#
#                 #inputs.append(cp.deepcopy(class_settings))
#                 #var = 'cala'
#                 #class_settings.__dict__[var] = np.float(SP['cala0'])
#                 #valnew = class_settings.__dict__[var]
#                 #labels.append(var+': '+format(valold,"0.2g")+'->'+format(valnew,"0.2g"))
#
#                 #inputs.append(cp.deepcopy(class_settings))
#                 #var = 'crhoc'
#                 #class_settings.__dict__[var] = np.float(SP['crhoc'])
#                 #valnew = class_settings.__dict__[var]
#                 #labels.append(var+': '+format(valold,"0.2g")+'->'+format(valnew,"0.2g"))
#
#     key = "CERES"
#     if ((kwargs == {}) or ((key in kwargs.keys()) and (kwargs[key]))):
#
#         CERES_start_date = dt.datetime(2000,3,1)
#         DT_CERES_START = (CERES_start_date + dt.timedelta(days=(int((class_settings.datetime - CERES_start_date ).days/61) * 61)))
#         DT_CERES_END   = DT_CERES_START +dt.timedelta(days=60)
#
#         input_fn = "/user/data/gent/gvo000/gvo00090/EXT/data/CERES/CERES_SYN1deg-1H_Terra-Aqua-MODIS_Ed4A_Subset_"+DT_CERES_START.strftime("%Y%m%d")+"-"+DT_CERES_END.strftime("%Y%m%d")+".nc"
#         print("Reading afternoon cloud cover for "+str(class_settings.datetime)+" from "+input_fn)
#
#         var = 'cc'
#
#         input_nc = nc4.Dataset(input_fn,'r')
#
#         idatetime = np.where(np.array(pcd.ncgetdatetime(input_nc))  >= class_settings.datetime)[0][0]
#         idatetime_end = np.where(np.array(pcd.ncgetdatetime(input_nc))  < (class_settings.datetime+dt.timedelta(hours=int(class_settings.runtime/3600.))))[0][-1]
#
#         ilat = np.where(input_nc.variables['lat'][:] >= class_settings.lat)[0][-1]
#         ilon = np.where(input_nc.variables['lon'][:] >= class_settings.lon)[0][0]
#         print(class_settings.lat,class_settings.lon)
#
#         class_settings.__dict__[var] = np.nanmean(input_nc.variables['cldarea_total_1h'][idatetime:idatetime_end,ilat,ilon])/100.
#
#         input_nc.close()
#


#     key = "GIMMS"
#     if ((kwargs == {}) or ((key in kwargs.keys()) and (kwargs[key]))):
#
#
#         input_fn = "/user/data/gent/gvo000/gvo00090/EXT/data/GIMMS/v2/LAI/gimms-3g.v2.lai.1981-2015_monmean.nc"
#         print("Reading Leag Area Index from "+input_fn)
#         var = 'LAI'
#
#         #plt.plot
#
#         input_nc = nc4.Dataset(input_fn,'r')
#
#         #idatetime = np.where(np.array(pcd.ncgetdatetime(input_nc))  >= class_settings.datetime)[0][0]
#         idatetime = np.where(np.array(pcd.ncgetdatetime(input_nc))  >= class_settings.datetime)[0][0]
#
#         ilatitude = np.where(input_nc.variables['lat'][:] >= class_settings.lat)[0][-1]
#         ilongitude = np.where(input_nc.variables['lon'][:] >= class_settings.lon)[0][0]
#
#         # divide by cveg, since it only reflects the LAI for the vegetation fraction and not for the entire (satellite) grid cell
#
#         print('Warning! Dividing by cveg, which is: '+str(class_settings.cveg))
#         tarray = np.array(input_nc.variables['LAI'][:,ilatitude,ilongitude])/class_settings.cveg
#
#         if np.isnan(tarray[idatetime]):
#             print("interpolating GIMMS cveg nan value")
#
#             mask = np.isnan(tarray)
#             if np.where(mask)[0].shape[0] < 0.25*mask.shape[0]:
#                 tarray[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), tarray[~mask])
#             else:
#                 print("Warning. Could not interpolate GIMMS cveg nan value")
#
#         class_settings.__dict__[var] = tarray[idatetime]
#
#         input_nc.close()
#
#     key = "IGBPDIS_ALPHA"
#     if ((kwargs == {}) or ((key in kwargs.keys()) and (kwargs[key]))):
#
#         var = 'alpha'
#
#         input_fn = "/user/data/gent/gvo000/gvo00090/EXT/data/IGBP-DIS/FRACTIONS_GLEAMv31a.nc"
#         print("Reading albedo from "+input_fn)
#
#         input_nc = nc4.Dataset(input_fn,'r')
#         ilat = np.where(input_nc.variables['lat'][:] >= class_settings.lat)[0][-1]
#         ilon = np.where(input_nc.variables['lon'][:] >= class_settings.lon)[0][0]
#
#
#         landfr = {}
#         for ltype in ['W','B','H','TC']:
#             landfr[ltype] = input_nc.variables['f'+ltype][0,ilon,ilat]
#
#         aweights = {'W':0.075,'TC':0.15,'H':0.22,'B':0.30}
#
#         alpha=0.
#         for ltype in landfr.keys():
#             alpha += landfr[ltype]*aweights[ltype]
#
#
#         class_settings.__dict__[var] = alpha
#         input_nc.close()
#
#
#     key = "ERAINT_ST"
#     if ((kwargs == {}) or ((key in kwargs.keys()) and (kwargs[key]))):
#
#         input_fn = '/user/data/gent/gvo000/gvo00090/EXT/data/ERA-INTERIM/by_var_nc/stl1_3hourly/stl1_'+str(class_settings.datetime.year)+"_3hourly.nc"
#         print("Reading soil temperature from "+input_fn)
#
#         var = 'Tsoil'
#         input_nc = nc4.Dataset(input_fn,'r')
#
#         idatetime = np.where(np.array(pcd.ncgetdatetime(input_nc))  >= class_settings.datetime)[0][0]
#
#         ilatitude = np.where(input_nc.variables['latitude'][:] >= class_settings.lat)[0][-1]
#         ilongitude = np.where(input_nc.variables['longitude'][:] >= class_settings.lon)[0][0]
#
#
#         class_settings.__dict__[var] = input_nc.variables['stl1'][idatetime,ilatitude,ilongitude]
#
#         input_fn = '/user/data/gent/gvo000/gvo00090/EXT/data/ERA-INTERIM/by_var_nc/stl2_3hourly/stl2_'+str(class_settings.datetime.year)+"_3hourly.nc"
#         var = 'T2'
#
#         input_nc = nc4.Dataset(input_fn,'r')
#
#         idatetime = np.where(np.array(pcd.ncgetdatetime(input_nc))  >= class_settings.datetime)[0][0]
#
#         ilatitude = np.where(input_nc.variables['latitude'][:] >= class_settings.lat)[0][-1]
#         ilongitude = np.where(input_nc.variables['longitude'][:] >= class_settings.lon)[0][0]
#
#
#         class_settings.__dict__[var] = input_nc.variables['stl2'][idatetime,ilatitude,ilongitude]
#
#
#         input_nc.close()
#
#
#
#     #inputs.append(cp.deepcopy(class_settings))
#     #var = 'T2'
#     #valold = class_settings.__dict__[var]
#     #
#     #class_settings.__dict__[var] = 305.
#     #class_settings.__dict__['Tsoil'] = 302.
#     #valnew = class_settings.__dict__[var]
#     #labels.append(var+': '+format(valold,"0.2g")+'->'+format(valnew,"0.2g"))
#
#
#
#     #inputs.append(cp.deepcopy(class_settings))
#     #
#     #var = 'Lambda'
#     #valold = class_settings.__dict__[var]
#
#     ## I presume that the skin layer conductivity scales with both LAI and vegetation fraction, which seems ~ valid according to table 10.6 in CLASS-book.
#     ## I need to ask Chiel.
#     ## I extrapolate from Lambda value of grass with Lambda = 5.9 W m-2 K-1, LAI = 2 and cveg = 0.85
#     #
#     #valnew = 5.9 / 2. / 0.85 * class_settings.__dict__['LAI'] * class_settings.__dict__['cveg']
#     #class_settings.__dict__[var] = valnew
#     #labels.append(var+': '+format(valold,"0.2g")+'->'+format(valnew,"0.2g"))
#
#
#
#     key = "GLAS"
#     if ((kwargs == {}) or ((key in kwargs.keys()) and (kwargs[key]))):
#
#         input_fn = "/user/data/gent/gvo000/gvo00090/EXT/data/GLAS/global_canopy_height_0.25.nc"
#         print("Reading canopy height for determining roughness length from "+input_fn)
#         var = 'z0m'
#
#
#         #plt.plot
#
#         input_nc = nc4.Dataset(input_fn,'r')
#
#         ilat = np.where(input_nc.variables['lat'][:] >= class_settings.lat)[0][0]
#         ilon = np.where(input_nc.variables['lon'][:] >= class_settings.lon)[0][0]
#
#         testval = np.float64(input_nc.variables['Band1'][ilat,ilon])/10.
#
#         lowerlimit = 0.01
#         if testval < lowerlimit:
#             print('forest canopy height very very small. We take a value of '+str(lowerlimit))
#             class_settings.__dict__[var] = lowerlimit
#         else:
#             class_settings.__dict__[var] = testval
#
#         class_settings.__dict__['z0h'] =  class_settings.__dict__['z0m']/10.
#
#
#         input_nc.close()





