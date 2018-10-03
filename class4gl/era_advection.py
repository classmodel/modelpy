

    self.get_idx_in_dataset(globaldata)
    xrin = globaldata['t'].page['t']

    # 3. prepare artificial xarray input datasets which will allow to make gradient calculations along the W-E directions with xarray on the fly with just one rule
    # 3.1 we make 'left' and 'right' datasets which will be substracted for calculating gradients
    # dataset of ilon = 1,2,3...
    xrleft = xrin.sel(lon=xrin.lon[:-25])
    # dataset of ilon = 0,1,2...
    xrright = xrin.sel(lon=xrin.lon[25:])
    
    # 3.2 The output will be on a staggered grid with the lon coordinate to the half-level calculated just hereafer in 3.3. 
    #     Still, we will need to full-level longitude values in the xarray dataset for calculating the grid spacing for the gradients.
    xrright['flon'] = (('lon',), xrright.lon.values)
    xrleft['flon'] = (('lon',), xrleft.lon.values)
    
    
    # 3.3 In order to make xarray doing the calculation advection correctly, the 'left' and 'right' values that we need on each grid cell requires equal underlying longitude coordinate values. 
    xrleft.lon.values = (xrleft.lon.values+xrright.lon.values)/2.
    xrright.lon.values = xrleft.lon.values


    # 4. We do similar preparations for S-N direction. Please note that the advection results for S-N and W-E direction are on different grids, that are also different from the original grid.
    xrbottom = xrin.sel(lat=xrin.lat[:-25])
    xrtop = xrin.sel(lat=xrin.lat[25:])
    xrtop['flat'] = (('lat',), xrtop.lat.values)
    xrbottom['flat'] = (('lat',), xrbottom.lat.values)
    xrbottom.lat.values = (xrbottom.lat.values+xrtop.lat.values)/2.
    xrtop.lat.values = xrbottom.lat.values
    
    
    dia_earth = 40000000.

    # for input variables (COSMO naming)
    VARS_COSMO = ['QV','U','V','T']
    # for output variables (ECMWF naming)
    vars_ECMWF = ['q','u','v','t']


    # some netcdf polishing: add units and description to netcdf output. 
    units = dict(
          advq_x='kg kg-1 s-1', advq_y='kg kg-1 s-1',
          advt_x='K s-1',       advt_y='K s-1',
          advu_x='m s-2',       advu_y='m s-2',
          advv_x='m s-2',       advv_y='m s-2',
          divU_x='s-1',         divU_y='s-1',
                )
    long_names = dict(
          advq_x='zonal advection of specific humidity',        
          advt_x='zonal advection of heat',                     
          advu_x='zonal advection of zonal wind component',     
          advv_x='zonal advection of meridional wind component',
          divU_x='horizontal wind divergence in the zonal direction',
          advq_y='meridional advection of specific humidity',
          advt_y='meridional advection of heat',                            
          advu_y='meridional advection of zonal wind component',            
          advv_y='meridional advection of meridional wind component',
          divU_y='horizontal wind divergence in the meridional direction',
                )
    #print((xrtop.flat - xrbottom.flat)/360.*dia_earth)
    # 5. loop over each variable

    # make the selections
    xrleft_sel = xrleft.isel(time=itimes,lat=ilats,lon=ilons)
    xrright_sel = xrright.isel(time=itimes,lat=ilats,lon=ilons)
    xrtop_sel = xrtop.isel(time=itimes,lat=ilats,lon=ilons)
    xrbottom_sel = xrbottom.isel(time=itimes,lat=ilats,lon=ilons)

    for ivar,var in enumerate(vars_ECMWF):
        VAR = VARS_COSMO[ivar]
    

        dims = globaldata.datasets[key].page[key].dims
        namesmean = list(dims)
        namesmean.remove('lev')
        idxmean = [dims.index(namemean) for namemean in namesmean]
        # over which dimensions we take a mean:
        dims = globaldata.datasets[key].page[key].dims
        namesmean = list(dims)
        namesmean.remove('lev')
        idxmean = [dims.index(namemean) for namemean in namesmean]
        #6. actual calculation for the W-E direction
        #######################################################
        print('calculation of advection')
        #######################################################

        if var == 't':
            self.update(source='era-interim_calc',pars={'adv'+var+'_x':\


        ( - (xrright_sel.p**(Rdcp) * xrright_sel.u*xrright[VAR] -
             xrleft_sel.p**(Rdcp) * xrleft_sel.u*rleft[VAR]) /\
                ((xrright.flon - xrleft.flon) /360.*dia_earth *np.cos(xrright.lat/180.*np.pi)) /\
                ((xrright.p**(Rdcp)+xrleft.p**(Rdcp))/2.)



        self.update(source='era-interim_calc',pars={'adv'+var+'_x':\
                               (- (xrright_sel.u*xrright_sel[var] - 
                                   xrleft_sel.u * xrleft_sel[var]) 
                                  /\
                                  ((xrright_sel.flon - xrleft_sel.flon) /360.*dia_earth 
                                   *np.cos(xrright_sel.lat/180.*np.pi))).mean(axis=tuple(idxmean)).values *1.\
        self.update(source='era-interim_calc',pars={'adv'+var+'_x':\
                                +\
                               (- (xrtop_sel.u*xrtop_sel[var] - 
                                   xrbottom_sel.u * xrbottom_sel[var]) 
                                  /\
                                  ((xrtop_sel.flon - xrbottom_sel.flon)\
                                   /360.*dia_earth).mean(axis=tuple(idxmean)).values *1.
                                                   })

































          # over which dimensions we take a mean:
          dims = globaldata.datasets[key].page[key].dims
          namesmean = list(dims)
          namesmean.remove('lev')
          idxmean = [dims.index(namemean) for namemean in namesmean]
          
          value = \
          globaldata.datasets[key].page[key].isel(time=itimes,
                                                  lat=ilats,lon=ilons).mean(axis=tuple(idxmean)).values * 1.

          # Ideally, source should be equal to the datakey of globaldata.library 
          # or globaldata.datasets (eg., DSMW, IGBP-DIS, ERA-INTERIM etc.) 
          #  but therefore the globaldata class requires a revision to make this work
          self.update(source='globaldata',air_ac=pd.DataFrame({key:list(value)})) 

