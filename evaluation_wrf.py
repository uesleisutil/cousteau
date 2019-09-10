#!/home/uesleisutil/anaconda3/bin/python python
# -*- coding: utf-8 -*-

"""
File name:      evaluation_wrf.py
Author:         Ueslei Adriano Sutil
Email:          uesleisutil1@gmail.com
Created:        03 April 2019
Last modified:  17 April 2019
Version:        2.0.1
Python:         3.7.1

Evaluate WRF output using:
    - Bias (Contour);
    - Root Mean Square Error (RMSE; Contour);
    - Normalized Root Mean Square Error (NRMSE; Contour);
    - Mean Absolute Error (MAE; Contour);

Compare: 
    - ERA5 (Hersbach et al., 2018; https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=overview):
        Temperature at 2 meters height (°C);
        Wind Speed at 10 m height (m/s);
        Sea Level Pressure (hPa);

    - CFSR (Saha et al., 2010; https://rda.ucar.edu/datasets/ds093.0/):
        Temperature at 2 meters height (°C);
        Wind Speed at 10 m height (m/s);
        Sea Level Pressure (hPa);
        
    - MSWEP (Beck et al., 2019; http://www.gloh2o.org/):
        Daily Precipitation (mm).

Post-process de ROMS simulation to match with the databases:
    cdo seltimestep,169/336 wrf.nc wrf_ts.nc
    cdo daymean wrf_ts.nc wrf_ts_daymean_mswep.nc
    cdo timselmean,6 wrf_ts.nc wrf_ts_6hour_cfsr.nc
"""
  
# Library import.
import numpy as np
import matplotlib.pyplot as plt
from   mpl_toolkits.basemap import Basemap
import netCDF4
from   roms_libs import *
import pyresample
import cmocean
import subprocess, os
from wrf import destagger,getvar
from sty import bg, rs
from progress.bar import IncrementalBar
matplotlib.use("Agg")

# Customizations.
bbox                = [-52.5, -45.5, -30.5, -25.5]
lonbounds           = [-52.5,-45.5] 
latbounds           = [-30.5,-25.5]
wrf_file          = '/media/ueslei/Ueslei/INPE/PCI/Projetos/SC_2008/Outputs/normal/wrf.nc'
wrf_file_cfsr     = '/media/ueslei/Ueslei/INPE/PCI/Projetos/SC_2008/Outputs/normal/wrf_6h_ncks.nc'
wrf_file_mswep    = '/media/ueslei/Ueslei/INPE/PCI/Projetos/SC_2008/Outputs/normal/wrf_daymean_ncks.nc'
era_file          = '/media/ueslei/Ueslei/INPE/PCI/Projetos/SC_2008/Dados/Evaluation/ERA5/era5.nc'
cfsr_file         = '/media/ueslei/Ueslei/INPE/PCI/Projetos/SC_2008/Dados/Evaluation/CFSR/cfsr.nc'
mswep_file        = '/media/ueslei/Ueslei/INPE/PCI/Projetos/SC_2008/Dados/Evaluation/MSWEP/mswep.nc'

clevs_bias_t2     = np.arange(-7,2.05,0.01)
ticks_bias_t2     = np.arange(min(clevs_bias_t2),max(clevs_bias_t2),1)
clevs_rmse_t2     = np.arange(0,9.05,0.01)
ticks_rmse_t2     = np.arange(min(clevs_rmse_t2),max(clevs_rmse_t2),1)
clevs_nrmse_t2    = np.arange(0,0.405,0.001)
ticks_nrmse_t2    = np.arange(min(clevs_nrmse_t2),max(clevs_nrmse_t2),0.05)
clevs_mae_t2      = np.arange(0,8.02,0.01)
ticks_mae_t2      = np.arange(min(clevs_mae_t2),max(clevs_mae_t2),1)

clevs_bias_wind   = np.arange(-5,4.1,0.01)
ticks_bias_wind   = np.arange(min(clevs_bias_wind),max(clevs_bias_wind),1)
clevs_rmse_wind   = np.arange(0,6.1,0.01)
ticks_rmse_wind   = np.arange(min(clevs_rmse_wind),max(clevs_rmse_wind),1)
clevs_nrmse_wind  = np.arange(0,3.03,0.001)
ticks_nrmse_wind  = np.arange(min(clevs_nrmse_wind),max(clevs_nrmse_wind),0.3)
clevs_mae_wind    = np.arange(0,5.02,0.01)
ticks_mae_wind    = np.arange(min(clevs_mae_wind),max(clevs_mae_wind),1)

clevs_bias_slp    = np.arange(-10,4.1,0.01)
ticks_bias_slp    = np.arange(min(clevs_bias_slp),max(clevs_bias_slp),1)
clevs_rmse_slp    = np.arange(0,14.1,0.01)
ticks_rmse_slp    = np.arange(min(clevs_rmse_slp),max(clevs_rmse_slp),1)
clevs_nrmse_slp   = np.arange(0,1.03,0.001)
ticks_nrmse_slp   = np.arange(min(clevs_nrmse_slp),max(clevs_nrmse_slp),0.3)
clevs_mae_slp     = np.arange(0,14.02,0.01)
ticks_mae_slp     = np.arange(min(clevs_mae_slp),max(clevs_mae_slp),1)

clevs_bias_prec   = np.arange(-15,225.1,0.1)
ticks_bias_prec   = np.arange(min(clevs_bias_prec),max(clevs_bias_prec),20)
clevs_rmse_prec   = np.arange(0,230.1,0.1)
ticks_rmse_prec   = np.arange(min(clevs_rmse_prec),max(clevs_rmse_prec),20)
clevs_nrmse_prec  = np.arange(0,230.03,0.1)
ticks_nrmse_prec  = np.arange(min(clevs_nrmse_prec),max(clevs_nrmse_prec),5)
clevs_mae_prec    = np.arange(0,5.02,0.01)
ticks_mae_prec    = np.arange(min(clevs_mae_prec),max(clevs_mae_prec),1)


print(bg.da_cyan+'Which data? (1) ERA5, (2) CFSR or (3) MSWEP?'+bg.rs)      
dataset  = input()
if dataset == '1':
    print(bg.da_cyan+'Evaluate which WRF-ARW variable? (1) Temperature at 2 meters, (2) Wind Speed at 10m, (3) Sea Level Pressure?'+bg.rs)
    contourf_var  = input()
    if contourf_var == '1':
        nc_era      = netCDF4.Dataset(era_file)
        lon_era     = nc_era.variables['longitude'][:]-360
        lat_era     = nc_era.variables['latitude'][:]
        latli       = np.argmin(np.abs(lat_era-latbounds[1]))
        latui       = np.argmin(np.abs(lat_era-latbounds[0])) 
        lonli       = np.argmin(np.abs(lon_era-lonbounds[0]))
        lonui       = np.argmin(np.abs(lon_era-lonbounds[1]))  
        lon_era     = nc_era.variables['longitude'][lonli:lonui]-360
        lat_era     = nc_era.variables['latitude'][latli:latui]
        lon_era,lat_era = np.meshgrid(lon_era,lat_era)
        era_lat_len = len(lat_era[:,0])
        era_lon_len = len(lon_era[0, :])
        era_loop    = len(nc_era.variables['time'][:])
        observed    = np.zeros([era_loop,era_lat_len,era_lon_len])
        bar         = IncrementalBar('Processing Air Temperature at 2 m from ERA5', max=era_loop)
        for i in range(0,era_loop):
            temp_era        = nc_era.variables['t2m'][i,latli:latui,lonli:lonui]-273.15
            observed[i,:,:] = temp_era
            bar.next()
        bar.finish()
                       
        nc_wrf      = netCDF4.Dataset(wrf_file)
        lon_wrf     = nc_wrf.variables['XLONG'][0,0,:]
        lat_wrf     = nc_wrf.variables['XLAT'][0,:,0]               
        latli       = np.argmin(np.abs(lat_wrf-latbounds[1]))
        latui       = np.argmin(np.abs(lat_wrf-latbounds[0])) 
        lonli       = np.argmin(np.abs(lon_wrf-lonbounds[0]))
        lonui       = np.argmin(np.abs(lon_wrf-lonbounds[1]))
        lon_wrf     = lon_wrf[lonli:lonui]
        lat_wrf     = lat_wrf[latui:latli]
        lon_wrf,lat_wrf = np.meshgrid(lon_wrf,lat_wrf)      
        wrf_loop    = len(nc_wrf.variables['Times'][:])
        wrf_lat_len = len(lat_wrf[:,0])
        wrf_lon_len = len(lon_wrf[0,:])
        orig_def    = pyresample.geometry.SwathDefinition(lons=lon_wrf, lats=lat_wrf)
        targ_def    = pyresample.geometry.SwathDefinition(lons=lon_era, lats=lat_era)
        expected    = np.zeros([era_loop,era_lat_len,era_lon_len])
        bar         = IncrementalBar('Processing Air Temperature at 2 m from WRF', max=wrf_loop)
        for i in range(0,wrf_loop):  
            temp_wrf        = nc_wrf.variables['T2'][i,latui:latli,lonli:lonui]-273.15
            expected[i,:,:] = pyresample.kd_tree.resample_gauss(orig_def, temp_wrf, targ_def, radius_of_influence=50000, neighbours=10,sigmas=25000, fill_value=None)
            bar.next()
        bar.finish()
   
    if contourf_var=='2':
        nc_era      = netCDF4.Dataset(era_file)
        lon_era     = nc_era.variables['longitude'][:]-360
        lat_era     = nc_era.variables['latitude'][:]
        latli       = np.argmin(np.abs(lat_era-latbounds[1]))
        latui       = np.argmin(np.abs(lat_era-latbounds[0])) 
        lonli       = np.argmin(np.abs(lon_era-lonbounds[0]))
        lonui       = np.argmin(np.abs(lon_era-lonbounds[1]))  
        lon_era     = nc_era.variables['longitude'][lonli:lonui]-360
        lat_era     = nc_era.variables['latitude'][latli:latui]
        lon_era,lat_era = np.meshgrid(lon_era,lat_era)
        era_lat_len = len(lat_era[:,0])
        era_lon_len = len(lon_era[0, :])
        era_loop    = len(nc_era.variables['time'][:])
        u10_era2    = np.zeros([era_loop,era_lat_len,era_lon_len])
        v10_era2    = np.zeros([era_loop,era_lat_len,era_lon_len]) 
        bar         = IncrementalBar('Processing Wind Magnitude at 10 m from ERA5', max=era_loop)      
        for i in range(0,era_loop):
            u10_era         = nc_era.variables['u10'][i,latli:latui,lonli:lonui]
            v10_era         = nc_era.variables['v10'][i,latli:latui,lonli:lonui]
            u10_era2[i,:,:] = u10_era
            v10_era2[i,:,:] = v10_era  
            bar.next()
        bar.finish()
        nc_wrf      = netCDF4.Dataset(wrf_file)
        lon_wrf     = nc_wrf.variables['XLONG'][0,0,:]
        lat_wrf     = nc_wrf.variables['XLAT'][0,:,0]               
        latli       = np.argmin(np.abs(lat_wrf-latbounds[1]))
        latui       = np.argmin(np.abs(lat_wrf-latbounds[0])) 
        lonli       = np.argmin(np.abs(lon_wrf-lonbounds[0]))
        lonui       = np.argmin(np.abs(lon_wrf-lonbounds[1]))
        lon_wrf     = lon_wrf[lonli:lonui]
        lat_wrf     = lat_wrf[latui:latli]
        lon_wrf,lat_wrf = np.meshgrid(lon_wrf,lat_wrf)      
        wrf_loop    = len(nc_wrf.variables['Times'][168:409])
        wrf_lat_len = len(lat_wrf[:,0])
        wrf_lon_len = len(lon_wrf[0,:])

        orig_def    = pyresample.geometry.SwathDefinition(lons=lon_wrf, lats=lat_wrf)
        targ_def    = pyresample.geometry.SwathDefinition(lons=lon_era, lats=lat_era)

        u10_wrf2   = np.zeros([era_loop,era_lat_len,era_lon_len])
        v10_wrf2   = np.zeros([era_loop,era_lat_len,era_lon_len])       
        bar        = IncrementalBar('Processing Wind Magnitude at 10 m from WRF', max=wrf_loop)    
        for i in range(0,wrf_loop):  
            u10_wrf         = nc_wrf.variables['U10'][i,latui:latli,lonli:lonui]
            u10_wrf2[i,:,:] = pyresample.kd_tree.resample_gauss(orig_def, u10_wrf, targ_def, radius_of_influence=50000, neighbours=10,sigmas=25000, fill_value=None)
            v10_wrf         = nc_wrf.variables['V10'][i,latui:latli,lonli:lonui]
            v10_wrf2[i,:,:] = pyresample.kd_tree.resample_gauss(orig_def, v10_wrf, targ_def, radius_of_influence=50000, neighbours=10,sigmas=25000, fill_value=None)
            bar.next()
        bar.finish()     

        expected = np.sqrt(u10_wrf2**2 + v10_wrf2**2)
        observed = np.sqrt(u10_era2**2 + v10_era2**2)

    if contourf_var=='3':
        nc_era      = netCDF4.Dataset(era_file)
        lon_era     = nc_era.variables['longitude'][:]-360
        lat_era     = nc_era.variables['latitude'][:]
        latli       = np.argmin(np.abs(lat_era-latbounds[1]))
        latui       = np.argmin(np.abs(lat_era-latbounds[0])) 
        lonli       = np.argmin(np.abs(lon_era-lonbounds[0]))
        lonui       = np.argmin(np.abs(lon_era-lonbounds[1]))  
        lon_era     = nc_era.variables['longitude'][lonli:lonui]-360
        lat_era     = nc_era.variables['latitude'][latli:latui]
        lon_era,lat_era = np.meshgrid(lon_era,lat_era)
        era_lat_len = len(lat_era[:,0])
        era_lon_len = len(lon_era[0, :])
        era_loop    = len(nc_era.variables['time'][:])
        observed    = np.zeros([era_loop,era_lat_len,era_lon_len])
        bar         = IncrementalBar('Processing Sea Level Pressure from ERA5', max=era_loop)    
        for i in range(0,era_loop):
            slp_era         = nc_era.variables['msl'][i,latli:latui,lonli:lonui]/100
            observed[i,:,:] = slp_era
            bar.next()
        bar.finish()  

        nc_wrf      = netCDF4.Dataset(wrf_file)
        lon_wrf     = nc_wrf.variables['XLONG'][0,0,:]
        lat_wrf     = nc_wrf.variables['XLAT'][0,:,0]               
        latli       = np.argmin(np.abs(lat_wrf-latbounds[1]))
        latui       = np.argmin(np.abs(lat_wrf-latbounds[0])) 
        lonli       = np.argmin(np.abs(lon_wrf-lonbounds[0]))
        lonui       = np.argmin(np.abs(lon_wrf-lonbounds[1]))
        lon_wrf     = lon_wrf[lonli:lonui]
        lat_wrf     = lat_wrf[latui:latli]
        lon_wrf,lat_wrf = np.meshgrid(lon_wrf,lat_wrf)      
        wrf_loop    = len(nc_wrf.variables['Times'][168:409])
        wrf_lat_len = len(lat_wrf[:,0])
        wrf_lon_len = len(lon_wrf[0,:])
        orig_def    = pyresample.geometry.SwathDefinition(lons=lon_wrf, lats=lat_wrf)
        targ_def    = pyresample.geometry.SwathDefinition(lons=lon_era, lats=lat_era)
        expected    = np.zeros([era_loop,era_lat_len,era_lon_len])
        bar         = IncrementalBar('Processing Sea Level Pressure from WRF', max=wrf_loop)    
        for i in range(0,wrf_loop):  
            slp_wrf         = getvar(nc_wrf, "slp", i, units="hPa",meta=False)[latui:latli,lonli:lonui]
            expected[i,:,:] = pyresample.kd_tree.resample_gauss(orig_def, slp_wrf, targ_def, radius_of_influence=50000, neighbours=10,sigmas=25000, fill_value=None)
            bar.next()
        bar.finish()  
        
if dataset == '2':
    print(bg.da_cyan+'Evaluate which WRF-ARW variable? (1) Temperature at 2 meters, (2) Wind Speed at 10m or (3) Sea Level Pressure?'+bg.rs)
    contourf_var  = input()
    if contourf_var=='1':
        nc_cfsr      = netCDF4.Dataset(cfsr_file)
        lon_cfsr     = nc_cfsr.variables['lon'][:]
        lat_cfsr     = nc_cfsr.variables['lat'][:]
        latli        = np.argmin(np.abs(lat_cfsr-latbounds[1]))
        latui        = np.argmin(np.abs(lat_cfsr-latbounds[0])) 
        lonli        = np.argmin(np.abs(lon_cfsr-lonbounds[0]))
        lonui        = np.argmin(np.abs(lon_cfsr-lonbounds[1]))
        lon_cfsr     = nc_cfsr.variables['lon'][lonli:lonui]
        lat_cfsr     = nc_cfsr.variables['lat'][latli:latui]
        lon_cfsr,lat_cfsr = np.meshgrid(lon_cfsr,lat_cfsr)
        cfsr_loop    = len(nc_cfsr.variables['time'][:])
        cfsr_lat_len = len(lat_cfsr[:,0])
        cfsr_lon_len = len(lon_cfsr[0, :])
        observed     = np.zeros([cfsr_loop,cfsr_lat_len,cfsr_lon_len])
        bar          = IncrementalBar('Processing Air Temperature at 2 m from CFSR', max=cfsr_loop)    
        for i in range(0,cfsr_loop):
            temp_cfsr       = nc_cfsr.variables['TMP_L103'][i,latli:latui,lonli:lonui]-273.15              
            observed[i,:,:] = temp_cfsr
            bar.next()
        bar.finish()  

        nc_wrf      = netCDF4.Dataset(wrf_file_cfsr)
        lon_wrf     = nc_wrf.variables['XLONG'][0,:]
        lat_wrf     = nc_wrf.variables['XLAT'][:,0]               
        latli       = np.argmin(np.abs(lat_wrf-latbounds[1]))
        latui       = np.argmin(np.abs(lat_wrf-latbounds[0])) 
        lonli       = np.argmin(np.abs(lon_wrf-lonbounds[0]))
        lonui       = np.argmin(np.abs(lon_wrf-lonbounds[1]))
        lon_wrf     = lon_wrf[lonli:lonui]
        lat_wrf     = lat_wrf[latui:latli]
        lon_wrf,lat_wrf = np.meshgrid(lon_wrf,lat_wrf)      
        wrf_loop    = nc_wrf.variables['Times'][168:409]
        wrf_loop    = wrf_loop[::6,0]
        wrf_loop    = len(wrf_loop)
        wrf_lat_len = len(lat_wrf[:,0])
        wrf_lon_len = len(lon_wrf[0,:])
        orig_def    = pyresample.geometry.SwathDefinition(lons=lon_wrf, lats=lat_wrf)
        targ_def    = pyresample.geometry.SwathDefinition(lons=lon_cfsr, lats=lat_cfsr)
        expected    = np.zeros([cfsr_loop,cfsr_lat_len,cfsr_lon_len])
        bar         = IncrementalBar('Processing Air Temperature at 2 m from WRF', max=wrf_loop)    
        for i in range(0,wrf_loop):  
            temp_wrf        = nc_wrf.variables['T2'][i,latui:latli,lonli:lonui]-273.15
            expected[i,:,:] = pyresample.kd_tree.resample_gauss(orig_def, temp_wrf, targ_def, radius_of_influence=50000, neighbours=10,sigmas=25000, fill_value=None)
            bar.next()
        bar.finish()

    if contourf_var=='2':
        nc_cfsr      = netCDF4.Dataset(cfsr_file)
        lon_cfsr     = nc_cfsr.variables['lon'][:]
        lat_cfsr     = nc_cfsr.variables['lat'][:]
        latli        = np.argmin(np.abs(lat_cfsr-latbounds[1]))
        latui        = np.argmin(np.abs(lat_cfsr-latbounds[0])) 
        lonli        = np.argmin(np.abs(lon_cfsr-lonbounds[0]))
        lonui        = np.argmin(np.abs(lon_cfsr-lonbounds[1]))
        lon_cfsr     = nc_cfsr.variables['lon'][lonli:lonui]
        lat_cfsr     = nc_cfsr.variables['lat'][latli:latui]
        lon_cfsr,lat_cfsr = np.meshgrid(lon_cfsr,lat_cfsr)
        cfsr_loop    = len(nc_cfsr.variables['time'][:])
        cfsr_lat_len = len(lat_cfsr[:,0])
        cfsr_lon_len = len(lon_cfsr[0, :])
        u10_cfsr2   = np.zeros([cfsr_loop,cfsr_lat_len,cfsr_lon_len])  
        v10_cfsr2   = np.zeros([cfsr_loop,cfsr_lat_len,cfsr_lon_len])  
        bar         = IncrementalBar('Processing Wing Magnitude at 10 m from CFSR', max=cfsr_loop)  
        for i in range(0,cfsr_loop):
            u10_cfsr         = nc_cfsr.variables['U_GRD_L103'][i,latli:latui,lonli:lonui]
            v10_cfsr         = nc_cfsr.variables['V_GRD_L103'][i,latli:latui,lonli:lonui]
            u10_cfsr2[i,:,:] = u10_cfsr
            v10_cfsr2[i,:,:] = v10_cfsr
            bar.next()
        bar.finish()

        nc_wrf      = netCDF4.Dataset(wrf_file_cfsr)
        lon_wrf     = nc_wrf.variables['XLONG'][0,0,:]
        lat_wrf     = nc_wrf.variables['XLAT'][0,:,0]               
        latli       = np.argmin(np.abs(lat_wrf-latbounds[1]))
        latui       = np.argmin(np.abs(lat_wrf-latbounds[0])) 
        lonli       = np.argmin(np.abs(lon_wrf-lonbounds[0]))
        lonui       = np.argmin(np.abs(lon_wrf-lonbounds[1]))
        lon_wrf     = lon_wrf[lonli:lonui]
        lat_wrf     = lat_wrf[latui:latli]
        lon_wrf,lat_wrf = np.meshgrid(lon_wrf,lat_wrf)      
        wrf_loop    = len(nc_wrf.variables['Times'][168:409])
        wrf_lat_len = len(lat_wrf[:,0])
        wrf_lon_len = len(lon_wrf[0,:])
        orig_def    = pyresample.geometry.SwathDefinition(lons=lon_wrf, lats=lat_wrf)
        targ_def    = pyresample.geometry.SwathDefinition(lons=lon_cfsr, lats=lat_cfsr)
        u10_wrf2   = np.zeros([cfsr_loop,cfsr_lat_len,cfsr_lon_len])
        v10_wrf2   = np.zeros([cfsr_loop,cfsr_lat_len,cfsr_lon_len])  
        bar         = IncrementalBar('Processing Wing Magnitude at 10 m from WRF', max=wrf_loop)       
        for i in range(0,wrf_loop):  
            u10_wrf         = nc_wrf.variables['U10'][i,latui:latli,lonli:lonui]
            u10_wrf2[i,:,:] = pyresample.kd_tree.resample_gauss(orig_def, u10_wrf, targ_def, radius_of_influence=50000, neighbours=10,sigmas=25000, fill_value=None)
            v10_wrf         = nc_wrf.variables['V10'][i,latui:latli,lonli:lonui]
            v10_wrf2[i,:,:] = pyresample.kd_tree.resample_gauss(orig_def, v10_wrf, targ_def, radius_of_influence=50000, neighbours=10,sigmas=25000, fill_value=None)
            bar.next()
        bar.finish()       
        expected = np.sqrt(u10_wrf2**2 + v10_wrf2**2)
        observed = np.sqrt(u10_cfsr2**2 + v10_cfsr2**2)

    if contourf_var=='3':
        nc_cfsr      = netCDF4.Dataset(cfsr_file)
        lon_cfsr     = nc_cfsr.variables['lon'][:]
        lat_cfsr     = nc_cfsr.variables['lat'][:]
        latli        = np.argmin(np.abs(lat_cfsr-latbounds[1]))
        latui        = np.argmin(np.abs(lat_cfsr-latbounds[0])) 
        lonli        = np.argmin(np.abs(lon_cfsr-lonbounds[0]))
        lonui        = np.argmin(np.abs(lon_cfsr-lonbounds[1]))
        lon_cfsr     = nc_cfsr.variables['lon'][lonli:lonui]
        lat_cfsr     = nc_cfsr.variables['lat'][latli:latui]
        lon_cfsr,lat_cfsr = np.meshgrid(lon_cfsr,lat_cfsr)
        cfsr_loop    = len(nc_cfsr.variables['time'][:])
        cfsr_lat_len = len(lat_cfsr[:,0])
        cfsr_lon_len = len(lon_cfsr[0, :])
        observed     = np.zeros([cfsr_loop,cfsr_lat_len,cfsr_lon_len])
        bar          = IncrementalBar('Processing Sea Level Pressure from CFSR', max=cfsr_loop)  
        for i in range(0,cfsr_loop):
            slp_cfsr        = nc_cfsr.variables['PRES_L101'][i,latli:latui,lonli:lonui]/100              
            observed[i,:,:] = slp_cfsr
            bar.next()
        bar.finish()   

        nc_wrf      = netCDF4.Dataset(wrf_file_cfsr)
        lon_wrf     = nc_wrf.variables['XLONG'][0,0,:]
        lat_wrf     = nc_wrf.variables['XLAT'][0,:,0]               
        latli       = np.argmin(np.abs(lat_wrf-latbounds[1]))
        latui       = np.argmin(np.abs(lat_wrf-latbounds[0])) 
        lonli       = np.argmin(np.abs(lon_wrf-lonbounds[0]))
        lonui       = np.argmin(np.abs(lon_wrf-lonbounds[1]))
        lon_wrf     = lon_wrf[lonli:lonui]
        lat_wrf     = lat_wrf[latui:latli]
        lon_wrf,lat_wrf = np.meshgrid(lon_wrf,lat_wrf)      
        wrf_loop    = len(nc_wrf.variables['Times'][168:409])
        wrf_lat_len = len(lat_wrf[:,0])
        wrf_lon_len = len(lon_wrf[0,:])
        orig_def    = pyresample.geometry.SwathDefinition(lons=lon_wrf, lats=lat_wrf)
        targ_def    = pyresample.geometry.SwathDefinition(lons=lon_cfsr, lats=lat_cfsr)
        expected    = np.zeros([cfsr_loop,cfsr_lat_len,cfsr_lon_len])
        bar         = IncrementalBar('Processing Sea Level Pressure from WRF', max=wrf_loop)  
        for i in range(0,wrf_loop):  
            slp_wrf         = getvar(nc_wrf, "slp", i, units="hPa",meta=False)[latui:latli,lonli:lonui]
            expected[i,:,:] = pyresample.kd_tree.resample_gauss(orig_def, slp_wrf, targ_def, radius_of_influence=50000, neighbours=10,sigmas=25000, fill_value=None)
            bar.next()
        bar.finish() 

if dataset=='3':
        contourf_var  = 1
        nc_mswep      = netCDF4.Dataset(mswep_file)
        lon_mswep     = nc_mswep.variables['lon'][:]
        lat_mswep     = nc_mswep.variables['lat'][:]
        latli_mswep   = np.argmin(np.abs(lat_mswep-latbounds[1]))
        latui_mswep   = np.argmin(np.abs(lat_mswep-latbounds[0])) 
        lonli_mswep   = np.argmin(np.abs(lon_mswep-lonbounds[0]))
        lonui_mswep   = np.argmin(np.abs(lon_mswep-lonbounds[1]))
        lon_mswep     = nc_mswep.variables['lon'][lonli_mswep:lonui_mswep]
        lat_mswep     = nc_mswep.variables['lat'][latli_mswep:latui_mswep]
        lon_mswep,lat_mswep = np.meshgrid(lon_mswep,lat_mswep)
        mswep_loop    = len(nc_mswep.variables['time'][:])
        mswep_lat_len = len(lat_mswep[:,0])
        mswep_lon_len = len(lon_mswep[0, :])
        observed      = np.zeros([mswep_loop,mswep_lat_len,mswep_lon_len])
        bar           = IncrementalBar('Processing Precipitation from MSWEP', max=mswep_loop)  
        for i in range(0,mswep_loop):
            rain_mswep      = nc_mswep.variables['precipitation'][i,latli_mswep:latui_mswep,lonli_mswep:lonui_mswep]           
            observed[i,:,:] = rain_mswep
            bar.next()
        bar.finish() 

        nc_wrf      = netCDF4.Dataset(wrf_file_mswep)
        lon_wrf     = nc_wrf.variables['XLONG'][0,:]
        lat_wrf     = nc_wrf.variables['XLAT'][:,0]               
        latli       = np.argmin(np.abs(lat_wrf-latbounds[1]))
        latui       = np.argmin(np.abs(lat_wrf-latbounds[0])) 
        lonli       = np.argmin(np.abs(lon_wrf-lonbounds[0]))
        lonui       = np.argmin(np.abs(lon_wrf-lonbounds[1]))
        lon_wrf     = lon_wrf[lonli:lonui]
        lat_wrf     = lat_wrf[latui:latli]
        lon_wrf,lat_wrf = np.meshgrid(lon_wrf,lat_wrf)      
        wrf_loop    = len(nc_wrf.variables['Times'][:])
        wrf_lat_len = len(lat_wrf[:,0])
        wrf_lon_len = len(lon_wrf[0,:])
        orig_def    = pyresample.geometry.SwathDefinition(lons=lon_wrf, lats=lat_wrf)
        targ_def    = pyresample.geometry.SwathDefinition(lons=lon_mswep, lats=lat_mswep)
        rainc_wrf2  = np.zeros([mswep_loop,mswep_lat_len,mswep_lon_len])
        rainnc_wrf2 = np.zeros([mswep_loop,mswep_lat_len,mswep_lon_len])   
        rain_wrf2   = np.zeros([mswep_loop,mswep_lat_len,mswep_lon_len])   
        bar         = IncrementalBar('Processing Precipitation from WRF', max=mswep_loop)  
        for i in range(0,wrf_loop):  
            rainc_wrf          =  nc_wrf.variables['RAINC'][i,latui:latli,lonli:lonui]
            rainc_wrf2[i,:,:]  = pyresample.kd_tree.resample_gauss(orig_def, rainc_wrf, targ_def, radius_of_influence=50000, neighbours=10,sigmas=25000, fill_value=None)
            rainnc_wrf         =  nc_wrf.variables['RAINNC'][i,latui:latli,lonli:lonui]
            rainnc_wrf2[i,:,:] = pyresample.kd_tree.resample_gauss(orig_def, rainnc_wrf, targ_def, radius_of_influence=50000, neighbours=10,sigmas=25000, fill_value=None)
            bar.next()
        bar.finish() 
        expected  = rainc_wrf2+rainnc_wrf2

print(bg.da_cyan+'Which statistical metric? (1) Bias, (2) RMSE, (3) NRMSE or (4) MAE.'+bg.rs)
metric  = input()
if metric=='1':
    val = np.mean(expected-observed,axis=0)
if metric=='2':
    differences = expected-observed
    differences_squared = differences ** 2 
    mean_of_differences_squared = np.average(differences_squared,axis=0)
    val = np.sqrt(mean_of_differences_squared)
if metric=='3':
    differences = expected-observed
    differences_squared = differences ** 2 
    mean_of_differences_squared = np.average(differences_squared,axis=0)
    rmse = np.sqrt(mean_of_differences_squared)
    val  = rmse/np.mean(observed,axis=0)
if metric=='4':
    error = expected-observed
    val   = np.mean(np.abs(error),axis=0)


# Create and plot map.
m    = Basemap(projection='merc',llcrnrlat=bbox[2],urcrnrlat=bbox[3],llcrnrlon=bbox[0],urcrnrlon=bbox[1], lat_ts=30,resolution='i')
fig  = plt.figure(1,figsize=(10,8))
plt.xlabel('Longitude'u' [\N{DEGREE SIGN}]',labelpad=18,size=10)
plt.ylabel('Latitude'u' [\N{DEGREE SIGN}]',labelpad=33,size=10)
ax   = fig.add_subplot(111)
m.drawparallels(np.arange(-90.,120.,1), linewidth=0.00, color='black', labels=[1,0,0,1],labelstyle="N/S",fontsize=10)
m.drawmeridians(np.arange(-180.,180.,1), linewidth=0.00,color='black', labels=[1,0,0,1],labelstyle="N/S",fontsize=10)
m.drawcountries(color = '#000000',linewidth=0.5)
m.drawcoastlines(color = '#000000',linewidth=0.5)

if metric=='1' and contourf_var=='1':
    clevs = clevs_bias_t2
    ticks = ticks_bias_t2
if metric=='2' and contourf_var=='1':
    clevs = clevs_rmse_t2
    ticks = ticks_rmse_t2
if metric=='3' and contourf_var=='1':
    clevs = clevs_nrmse_t2
    ticks = ticks_nrmse_t2
if metric=='4' and contourf_var=='1':
    clevs = clevs_mae_t2
    ticks = ticks_mae_t2
if metric=='1' and contourf_var=='2':
    clevs = clevs_bias_wind
    ticks = ticks_bias_wind
if metric=='2' and contourf_var=='2':
    clevs = clevs_rmse_wind
    ticks = ticks_rmse_wind
if metric=='3' and contourf_var=='2':
    clevs = clevs_nrmse_wind
    ticks = ticks_nrmse_wind
if metric=='4' and contourf_var=='2':
    clevs = clevs_mae_wind
    ticks = ticks_mae_wind
if metric=='1' and contourf_var=='3':
    clevs = clevs_bias_slp
    ticks = ticks_bias_slp
if metric=='2' and contourf_var=='3':
    clevs = clevs_rmse_slp
    ticks = ticks_rmse_slp
if metric=='3' and contourf_var=='3':
    clevs = clevs_nrmse_slp
    ticks = ticks_nrmse_slp
if metric=='4' and contourf_var=='3':
    clevs = clevs_mae_slp
    ticks = ticks_mae_slp
if metric=='1' and dataset=='3':
    clevs = clevs_bias_prec
    ticks = ticks_bias_prec
if metric=='2' and dataset=='3':
    clevs = clevs_rmse_prec
    ticks = ticks_rmse_prec
if metric=='3' and dataset=='3':
    clevs = clevs_nrmse_prec
    ticks = ticks_nrmse_prec
if metric=='4' and dataset=='3':
    clevs = clevs_mae_prec
    ticks = ticks_mae_prec

if dataset == '1':
    if metric=='2' or metric=='3' or metric == '4':
        cmap  = matplotlib.pyplot.jet()
        h1    = m.contourf(lon_era, lat_era, val, clevs,latlon=True,cmap=cmap)
    if metric=='1':
        cmap  = cmocean.cm.balance
        h1    = m.contourf(lon_era, lat_era, val, clevs,latlon=True,cmap=cmap,norm=MidpointNormalize(midpoint=0))          
if dataset =='2':
    if metric=='2' or metric=='3' or metric == '4':
        cmap  = matplotlib.pyplot.jet()
        h1    = m.contourf(lon_cfsr, lat_cfsr, val, clevs,latlon=True,cmap=cmap)   
    if metric=='1':
        cmap  = cmocean.cm.balance
        h1    = m.contourf(lon_cfsr, lat_cfsr, val, clevs,latlon=True,cmap=cmap,norm=MidpointNormalize(midpoint=0)) 
if dataset == '3':
    if metric=='2' or metric=='3' or metric == '4':
        cmap  = matplotlib.pyplot.jet()
        h1    = m.contourf(lon_mswep, lat_mswep, val, clevs,latlon=True,cmap=cmap)   
    if metric=='1':
        cmap  = cmocean.cm.balance
        h1    = m.contourf(lon_mswep, lat_mswep, val, clevs,latlon=True,cmap=cmap,norm=MidpointNormalize(midpoint=0)) 
cax   = fig.add_axes([0.37, 0.025, 0.27, 0.025])     
cb    = fig.colorbar(h1, cax=cax, orientation="horizontal",panchor=(0.5,0.5),shrink=0.3,ticks=ticks,pad=-10.5)

if metric=='1':
    if contourf_var=='1':
        cb.set_label(r'Air Temperature at 2 meters Bias [$^\circ\!$C]', fontsize=9, color='0.2',labelpad=0)
    if contourf_var=='2':
        cb.set_label(r'Wind Speed at 10 meters Bias [m.s⁻¹]', fontsize=9, color='0.2',labelpad=0)
    if contourf_var=='3':
         cb.set_label(r'Sea Level Pressure Bias [hPa]', fontsize=9, color='0.2',labelpad=0) 
    if dataset=='3':
        cb.set_label(r'Daily Precipitation Bias [mm]', fontsize=9, color='0.2',labelpad=0)                
if metric=='2':
    if contourf_var=='1':
        cb.set_label(r'Air Temperature at 2 meters Root Mean Square Error [$^\circ\!$C]', fontsize=5, color='0.2',labelpad=-0.2)
    if contourf_var=='2':
        cb.set_label(r'Wind Speed at 10 meters Root Mean Square Error [m.s⁻¹]', fontsize=5, color='0.2',labelpad=-0.2)
    if contourf_var=='3':
        cb.set_label(r'Sea Level Pressure Root Mean Square Error [hPa]', fontsize=5, color='0.2',labelpad=-0.2)  
    if dataset=='3':
        cb.set_label(r'Daily Precipitation Root Mean Square Error [mm]', fontsize=5, color='0.2',labelpad=-0.2)               
if metric=='3':
    if contourf_var=='1':
        cb.set_label(r'Air Temperature at 2 meters Normalized Root Mean Square Error [$^\circ\!$C]', fontsize=5, color='0.2',labelpad=-0.2)
    if contourf_var=='2':
        cb.set_label(r'Wind Speed at 10 meters Normalized Root Mean Square Error [m.s⁻¹]', fontsize=5, color='0.2',labelpad=-0.2)
    if contourf_var=='3':
        cb.set_label(r'Sea Level Pressure Normalized Root Mean Square Error [hPa]', fontsize=5, color='0.2',labelpad=-0.2)  
    if dataset=='3':
        cb.set_label(r'Daily Precipitation Normalized Root Mean Square Error [mm]', fontsize=5, color='0.2',labelpad=-0.2)               
if metric=='4':
    if contourf_var=='1':
        cb.set_label(r'Air Temperature at 2 meters Mean Absolute Error [$^\circ\!$C]', fontsize=5, color='0.2',labelpad=-0.2)
    if contourf_var=='2':
        cb.set_label(r'Wind Speed at 10 meters Mean Absolute Error [m.s⁻¹]', fontsize=5, color='0.2',labelpad=-0.2)
    if contourf_var=='3':
        cb.set_label(r'Sea Level Pressure Mean Absolute Error [hPa]', fontsize=5, color='0.2',labelpad=-0.2)  
    if dataset=='3':
        cb.set_label(r'Daily Precipitation Mean Absolute Error [mm]', fontsize=5, color='0.2',labelpad=-0.2)               

cb.ax.tick_params(labelsize=9, length=2, color='0.2', labelcolor='0.2',direction='in') 
cb.set_ticks(ticks)
try:
    os.makedirs("wrf_evaluation")
except FileExistsError:
    pass 
if metric == '1' and dataset == '1' and contourf_var == '1':
        plt.savefig('./wrf_evaluation/t2_bias_wrf_era5.png', transparent=False, bbox_inches = 'tight', pad_inches=0, dpi=250)              
if metric == '1' and dataset == '1' and contourf_var == '2':
        plt.savefig('./wrf_evaluation/wind_bias_wrf_era5.png', transparent=False, bbox_inches = 'tight', pad_inches=0, dpi=250)              
if metric == '1' and dataset == '1' and contourf_var == '3':
        plt.savefig('./wrf_evaluation/slp_bias_wrf_era5.png', transparent=False, bbox_inches = 'tight', pad_inches=0, dpi=250)              
if metric == '2' and dataset == '1' and contourf_var == '1':
        plt.savefig('./wrf_evaluation/t2_rmse_wrf_era5.png', transparent=False, bbox_inches = 'tight', pad_inches=0, dpi=250)              
if metric == '2' and dataset == '1' and contourf_var == '2':
        plt.savefig('./wrf_evaluation/wind_rmse_wrf_era5.png', transparent=False, bbox_inches = 'tight', pad_inches=0, dpi=250)              
if metric == '2' and dataset == '1' and contourf_var == '3':
        plt.savefig('./wrf_evaluation/slp_rmse_wrf_era5.png', transparent=False, bbox_inches = 'tight', pad_inches=0, dpi=250)              
if metric == '3' and dataset == '1' and contourf_var == '1':
        plt.savefig('./wrf_evaluation/t2_nrmse_wrf_era5.png', transparent=False, bbox_inches = 'tight', pad_inches=0, dpi=250)              
if metric == '3' and dataset == '1' and contourf_var == '2':
        plt.savefig('./wrf_evaluation/wind_nrmse_wrf_era5.png', transparent=False, bbox_inches = 'tight', pad_inches=0, dpi=250)              
if metric == '3' and dataset == '1' and contourf_var == '3':
        plt.savefig('./wrf_evaluation/slp_nrmse_wrf_era5.png', transparent=False, bbox_inches = 'tight', pad_inches=0, dpi=250)              
if metric == '4' and dataset == '1' and contourf_var == '1':
        plt.savefig('./wrf_evaluation/t2_mae_wrf_era5.png', transparent=False, bbox_inches = 'tight', pad_inches=0, dpi=250)              
if metric == '4' and dataset == '1' and contourf_var == '2':
        plt.savefig('./wrf_evaluation/wind_mae_wrf_era5.png', transparent=False, bbox_inches = 'tight', pad_inches=0, dpi=250)              
if metric == '4' and dataset == '1' and contourf_var == '3':
        plt.savefig('./wrf_evaluation/slp_mae_wrf_era5.png', transparent=False, bbox_inches = 'tight', pad_inches=0, dpi=250)    
if metric == '1' and dataset == '2' and contourf_var == '1':
        plt.savefig('./wrf_evaluation/t2_bias_wrf_cfsr.png', transparent=False, bbox_inches = 'tight', pad_inches=0, dpi=250)              
if metric == '1' and dataset == '2' and contourf_var == '2':
        plt.savefig('./wrf_evaluation/wind_bias_wrf_cfsr.png', transparent=False, bbox_inches = 'tight', pad_inches=0, dpi=250)              
if metric == '1' and dataset == '2' and contourf_var == '3':
        plt.savefig('./wrf_evaluation/slp_bias_wrf_cfsr.png', transparent=False, bbox_inches = 'tight', pad_inches=0, dpi=250)              
if metric == '2' and dataset == '2' and contourf_var == '1':
        plt.savefig('./wrf_evaluation/t2_rmse_wrf_cfsr.png', transparent=False, bbox_inches = 'tight', pad_inches=0, dpi=250)              
if metric == '2' and dataset == '2' and contourf_var == '2':
        plt.savefig('./wrf_evaluation/wind_rmse_wrf_cfsr.png', transparent=False, bbox_inches = 'tight', pad_inches=0, dpi=250)              
if metric == '2' and dataset == '2' and contourf_var == '3':
        plt.savefig('./wrf_evaluation/slp_rmse_wrf_cfsr.png', transparent=False, bbox_inches = 'tight', pad_inches=0, dpi=250)              
if metric == '3' and dataset == '2' and contourf_var == '1':
        plt.savefig('./wrf_evaluation/t2_nrmse_wrf_cfsr.png', transparent=False, bbox_inches = 'tight', pad_inches=0, dpi=250)              
if metric == '3' and dataset == '2' and contourf_var == '2':
        plt.savefig('./wrf_evaluation/wind_nrmse_wrf_cfsr.png', transparent=False, bbox_inches = 'tight', pad_inches=0, dpi=250)              
if metric == '3' and dataset == '2' and contourf_var == '3':
        plt.savefig('./wrf_evaluation/slp_nrmse_wrf_cfsr.png', transparent=False, bbox_inches = 'tight', pad_inches=0, dpi=250)              
if metric == '4' and dataset == '2' and contourf_var == '1':
        plt.savefig('./wrf_evaluation/t2_mae_wrf_cfsr.png', transparent=False, bbox_inches = 'tight', pad_inches=0, dpi=250)              
if metric == '4' and dataset == '2' and contourf_var == '2':
        plt.savefig('./wrf_evaluation/wind_mae_wrf_cfsr.png', transparent=False, bbox_inches = 'tight', pad_inches=0, dpi=250)              
if metric == '4' and dataset == '2' and contourf_var == '3':
        plt.savefig('./wrf_evaluation/slp_mae_wrf_cfsr.png', transparent=False, bbox_inches = 'tight', pad_inches=0, dpi=250)      
if metric == '1' and dataset == '3':
        plt.savefig('./wrf_evaluation/prec_bias_wrf_mswep.png', transparent=False, bbox_inches = 'tight', pad_inches=0, dpi=250)  
if metric == '2' and dataset == '3':
        plt.savefig('./wrf_evaluation/prec_rmse_wrf_mswep.png', transparent=False, bbox_inches = 'tight', pad_inches=0, dpi=250)  
if metric == '3' and dataset == '3':
        plt.savefig('./wrf_evaluation/prec_nrmse_wrf_mswep.png', transparent=False, bbox_inches = 'tight', pad_inches=0, dpi=250)  
if metric == '4' and dataset == '3':
        plt.savefig('./wrf_evaluation/prec_mae_wrf_mswep.png', transparent=False, bbox_inches = 'tight', pad_inches=0, dpi=250) 


