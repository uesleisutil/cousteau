#!/home/uesleisutil/anaconda3/bin/python python
# -*- coding: utf-8 -*-
"""
File name:      roms_horizontal_evaluation.py
Author:         Ueslei Adriano Sutil
Email:          uesleisutil1@gmail.com
Created:        22 March 2019
Last modified:  03 April 2019
Version:        2.0
License:        GPL
Python:         3.7.1

Evaluate ROMS output using Root-Mean-Square-Error (RMSE; Contour), Mean-Absolute-Percentage-Error (MAPE; Contour) and Bias(Contour):
Compare: 
         - MUR (Chin et al., 2017; https://podaac.jpl.nasa.gov/Multi-scale_Ultra-high_Resolution_MUR-SST):
           Sea Surface Temperature (°C)

         - OSCAR (Bonjean & Lagerloef, 2002; https://podaac.jpl.nasa.gov/dataset/OSCAR_L4_OC_third-deg):
           Current Speed at surface (m.s⁻¹) 

         - GLORYS12V1 (Fernandez & Lellouch, 2018;http://marine.copernicus.eu/services-portfolio/access-to-products/?option=com_csw&view=details&product_id=GLOBAL_REANALYSIS_PHY_001_030):
           Sea Surface Temperature (°C) 
           Current Speed at surface (m.s⁻¹)
           
Post-process de ROMS simulation to match with the databases.
What I do:
    - MUR and GLORYS (Daily data):
        ncks -v temp,u,v -d s_rho,49,49 roms_avg.nc roms_evaluation.nc
        cdo daymean roms_evaluation.nc roms_evaluation_mean.nc
        cdo splitday roms_evaluation_mean.nc roms_evaluation_mean
        cdo cat roms_evaluation_mean01 roms_evaluation_mean02 ... roms_evaluation_final.nc
    - OSCAR (Each 5 days data):
        ncks -v u,v -d s_rho,49,49 roms_avg.nc roms_evaluation.nc
        cdo daymean roms_evaluation.nc roms_evaluation_mean.nc
        cdo splitday roms_evaluation_mean.nc roms_evaluation_mean
        cdo cat roms_evaluation_mean21 roms_evaluation_mean26 roms_evaluation_final.nc
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
from sty import bg, rs
from progress.bar import IncrementalBar
matplotlib.use("Agg")

# Customizations.
print(bg.da_cyan+'Which project? (1) SC_2008, (2) ATLEQ or (3) Antartic.'+bg.rs)
project = input()
if project =='1':
    bbox                = [-53, -40, -32, -23]
    lonbounds           = [-53,-40] 
    latbounds           = [-32,-23]

    clevs_rmse_sst      = np.arange(0,3.51,0.01)
    ticks_rmse_sst      = np.array([0,0.5,1,1.5,2,2.5,3,3.5,4]) 
    clevs_mape_sst      = np.arange(0,21,0.1)
    ticks_mape_sst      = np.array([0,5,10,15,20]) 
    clevs_bias_sst      = np.arange(-4,4.01,0.01)
    ticks_bias_sst      = np.array([-4,-3,-2,-1,0,1,2,3,4]) 

    clevs_rmse_cur      = np.arange(0,0.751,0.01)
    ticks_rmse_cur      = np.array([0,0.25,0.50,0.75]) 
    clevs_mape_cur      = np.arange(0,51,0.1)
    ticks_mape_cur      = np.array([0,5,10,15,20,25,30,35,40,45,50]) 
    clevs_bias_cur      = np.arange(-0.51,0.51,0.01)
    ticks_bias_cur      = np.array([0.5,-0.25,0,0.25,0.5])    

    roms_mur_glorys_dir = '/media/ueslei/Ueslei/INPE/PCI/Projetos/SC_2008/Outputs/normal/roms_ts_daymean.nc'
    roms_oscar_dir      = '/media/ueslei/Ueslei/INPE/PCI/Projetos/SC_2008/Outputs/normal/roms_oscar_eval.nc'
    mur_dir             = '/media/ueslei/Ueslei/INPE/PCI/Projetos/SC_2008/Dados/Evaluation/MUR/mur.nc' 
    glorys_dir          = '/media/ueslei/Ueslei/INPE/PCI/Projetos/SC_2008/Dados/Evaluation/Glorys/glorys.nc'  
    oscar_dir           = '/media/ueslei/Ueslei/INPE/PCI/Projetos/SC_2008/Dados/Evaluation/OSCAR/oscar.nc'  

if project=='2':
    bbox                = [-60,15,-15.1,15.1]
    lonbounds           = [-60,15] 
    latbounds           = [-15.1,15.1]
    clevs_rmse_sst      = np.arange(0,4.01,0.01)
    ticks_rmse_sst      = np.array([0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5]) 
    clevs_mae_sst       = np.arange(-4,4.05,0.01)
    ticks_mae_sst       = np.arange(min(clevs_mae_sst),max(clevs_mae_sst),1)

    clevs_rmse_cur      = np.arange(0,1.11,0.01)
    ticks_rmse_cur      = np.array([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1]) 
    clevs_mae_cur       = np.arange(-1,1.01,0.01)
    ticks_mae_cur       = np.array([-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1]) 

    roms_mur_glorys_dir = '/media/ueslei/Ueslei/INPE/PCI/ATLEQ/Outputs/roms.nc'
    glorys_dir          = '/media/ueslei/Ueslei/INPE/PCI/ATLEQ/Dados/Evaluation/Glorys/glorys.nc'  
    
print(bg.da_cyan+'Evaluate which ROMS variable? (1) Sea Surface Temperature or (2) Current Speed at Surface.'+bg.rs)
contourf_var  = input()
if contourf_var=='1':
    # Open ROMS file then calculate its average through time.
    roms_file    = roms_mur_glorys_dir
    nc_roms      = netCDF4.Dataset(roms_file)
    lon_rho      = nc_roms.variables['lon_rho'][:]
    lat_rho      = nc_roms.variables['lat_rho'][:]
    i0,i1,j0,j1  = bbox2ij(lon_rho,lat_rho,bbox)
    lon_roms     = lon_rho[j0:j1, i0:i1]
    lat_roms     = lat_rho[j0:j1, i0:i1]
    
    print(bg.da_cyan+'Evaluate Sea Surface Temperature using: (1) MUR or (2) Glorys.'+bg.rs)
    choose_data  = input()
    if choose_data == '1':
        roms_loop    = len(nc_roms.variables['ocean_time'][:])
        roms_lat_len = len(lat_rho[j0:j1,0])
        roms_lon_len = len(lon_rho[0, i0:i1])
        temp_roms2   = np.zeros([roms_loop,roms_lat_len,roms_lon_len])
        bar         = IncrementalBar('Processing Sea Surface Temperature from ROMS', max=roms_loop)
        for i in range(0,roms_loop):
            temp_roms         = nc_roms.variables['temp'][i,0,j0:j1, i0:i1]
            temp_roms2[i,:,:] = temp_roms
            temp_roms2[temp_roms2==0] = np.nan
            bar.next()
        bar.finish()

    if choose_data == '2':
        pass
    if choose_data=='1':
        # Open MUR file then calculate its average through time.
        mur_file    = mur_dir 
        nc_mur      = netCDF4.Dataset(mur_file)
        lon_mur     = nc_mur.variables['lon'][:]
        lat_mur     = nc_mur.variables['lat'][:]
        latli       = np.argmin( np.abs( lat_mur - latbounds[1] ) )
        latui       = np.argmin( np.abs( lat_mur - latbounds[0] ) ) 
        lonli       = np.argmin( np.abs( lon_mur - lonbounds[0] ) )
        lonui       = np.argmin( np.abs( lon_mur - lonbounds[1] ) )  
        lon_mur     = nc_mur.variables['lon'][lonli:lonui]
        lat_mur     = nc_mur.variables['lat'][latui:latli]
        lon_mur,lat_mur = np.meshgrid(lon_mur,lat_mur)
        mur_loop    = len(nc_mur.variables['time'][:])    
        mur_lat_len = len(lat_rho[j0:j1,0])
        mur_lon_len = len(lon_rho[0, i0:i1])
        temp_mur2   = np.zeros([mur_loop,mur_lat_len,mur_lon_len])
        targ_def    = pyresample.geometry.SwathDefinition(lons=lon_roms, lats=lat_roms)
        orig_def    = pyresample.geometry.SwathDefinition(lons=lon_mur, lats=lat_mur)
        bar         = IncrementalBar('Processing Sea Surface Temperature from MUR', max=mur_loop)        
        for i in range(0,mur_loop):
            temp_mur         = nc_mur.variables['analysed_sst'][i,latui:latli,lonli:lonui]
            temp_mur         = temp_mur-273.15
            temp_mur         = pyresample.kd_tree.resample_gauss(orig_def, temp_mur, targ_def, radius_of_influence=500000, neighbours=20,sigmas=250000, fill_value=None)
            temp_mur2[i,:,:] = temp_mur
            temp_mur2[temp_mur2==0] = np.nan
            bar.next()
        bar.finish()

    if choose_data=='2':
        glorys_file    = glorys_dir 
        nc_glorys      = netCDF4.Dataset(glorys_file)
        lon_glorys     = nc_glorys.variables['longitude'][:]
        lat_glorys     = nc_glorys.variables['latitude'][:]
        latli          = np.argmin( np.abs( lat_glorys - latbounds[1] ) )
        latui          = np.argmin( np.abs( lat_glorys - latbounds[0] ) ) 
        lonli          = np.argmin( np.abs( lon_glorys - lonbounds[0] ) )
        lonui          = np.argmin( np.abs( lon_glorys - lonbounds[1] ) )  
        lon_glorys     = nc_glorys.variables['longitude'][lonli:lonui]
        lat_glorys     = nc_glorys.variables['latitude'][latui:latli]
        lon_glorys,lat_glorys = np.meshgrid(lon_glorys,lat_glorys)
        glorys_loop    = len(nc_glorys.variables['time'][:])
        glorys_lat_len = len(lat_glorys[:,0])
        glorys_lon_len = len(lon_glorys[0, :])
        temp_glorys2   = np.zeros([glorys_loop,glorys_lat_len,glorys_lon_len])

        roms_loop      = len(nc_roms.variables['ocean_time'][:])          
        roms_lat_len   = len(lat_glorys[:,0])
        roms_lon_len   = len(lon_glorys[0, :])
        temp_roms2     = np.zeros([roms_loop,roms_lat_len,roms_lon_len])
    
        orig_def       = pyresample.geometry.SwathDefinition(lons=lon_roms, lats=lat_roms)
        targ_def       = pyresample.geometry.SwathDefinition(lons=lon_glorys, lats=lat_glorys)
        bar            = IncrementalBar('Processing Sea Surface Temperature from MUR', max=glorys_loop)  
        for i in range(0,glorys_loop):
            temp_glorys         = nc_glorys.variables['thetao'][i,0,latui:latli,lonli:lonui]
            temp_glorys2[i,:,:] = temp_glorys
            temp_roms           = nc_roms.variables['temp'][i,0,j0:j1, i0:i1]
            temp_roms           = pyresample.kd_tree.resample_gauss(orig_def, temp_roms, targ_def, radius_of_influence=500000, neighbours=20,sigmas=250000, fill_value=None)
            temp_roms2[i,:,:]   = temp_roms
            bar.next()
        bar.finish()

    print(bg.da_cyan+'Which statistical metric? (1) Root Mean Square Error, (2) Mean Absolute Percentage Error or (3) Bias.'+bg.rs)
    metric  = input()
    if metric=='1':
        if choose_data == '1':
            differences         = temp_roms2-temp_mur2
            differences_squared = differences ** 2 
            mean_of_differences_squared = np.average(differences_squared,axis=0)
            val                 = np.sqrt(mean_of_differences_squared)
        if choose_data == '2':
            differences         = temp_roms2-temp_glorys2
            differences_squared = differences ** 2 
            mean_of_differences_squared = np.average(differences_squared,axis=0)
            val                 = np.sqrt(mean_of_differences_squared)
    if metric=='2':
        if choose_data =='1':
            val = np.abs((temp_mur2-temp_roms2)/temp_mur2).mean(axis=0)*100
        if choose_data =='2':
            val = np.abs((temp_glorys2 - temp_roms2) / temp_glorys2).mean(axis=0)*100 
    if metric=='3':
        if choose_data =='1':
             mean_roms  = np.average(temp_roms2,axis=0) 
             mean_mur   = np.average(temp_mur2,axis=0)
             val        = mean_roms-mean_mur   
        if choose_data =='2':
             mean_roms   = np.average(temp_roms2,axis=0) 
             mean_glorys = np.average(temp_glorys2,axis=0)
             val         = mean_roms-mean_glorys

    # Create and plot map.
    m    = Basemap(projection='merc',llcrnrlat=bbox[2],urcrnrlat=bbox[3],llcrnrlon=bbox[0],urcrnrlon=bbox[1], lat_ts=10,resolution='i')
    fig  = plt.figure(1,figsize=(10,8))
    plt.xlabel('Longitude'u' [\N{DEGREE SIGN}]',labelpad=18,size=10)
    plt.ylabel('Latitude'u' [\N{DEGREE SIGN}]',labelpad=33,size=10)
    ax   = fig.add_subplot(111)
    m.drawparallels(np.arange(-90.,120.,1.), linewidth=0.00, color='black', labels=[1,0,0,1],labelstyle="N/S",fontsize=10)
    m.drawmeridians(np.arange(-180.,180.,1.), linewidth=0.00,color='black', labels=[1,0,0,1],labelstyle="N/S",fontsize=10)
    m.drawcountries(color = '#ffffff',linewidth=0.5)
    m.drawcoastlines(color = '#ffffff',linewidth=0.5)
    m.fillcontinents(color = '#000000')
    if metric=='1':
        clevs = clevs_rmse_sst
        ticks = ticks_rmse_sst
        cmap  = cmocean.cm.thermal
    if metric=='2':
        clevs = clevs_mape_sst
        ticks = ticks_mape_sst 
        cmap  = cmocean.cm.thermal
    if metric=='3':
        clevs = clevs_bias_sst
        ticks = ticks_bias_sst 
        cmap  = cmocean.cm.balance

    if choose_data == '1':
        if metric=='1' or metric=='2':
            h1    = m.contourf(lon_roms, lat_roms, val, clevs,latlon=True,cmap=cmap,extend="both")
        if metric=='3':
            h1    = m.contourf(lon_roms, lat_roms, val, clevs,latlon=True,cmap=cmap,norm=MidpointNormalize(midpoint=0),extend="both") 
    if choose_data =='2':
        if metric =='1':
            h1    = m.contourf(lon_glorys, lat_glorys, val, clevs,latlon=True,cmap=cmap,extend="both") 
        if metric =='2' or metric =='3':
            h1    = m.contourf(lon_glorys, lat_glorys, val, clevs,latlon=True,cmap=cmap,norm=MidpointNormalize(midpoint=0),extend="both")      

    cax   = fig.add_axes([0.37, 0.025, 0.27, 0.025])     
    cb    = fig.colorbar(h1, cax=cax, orientation="horizontal",panchor=(0.5,0.5),shrink=0.3,ticks=ticks)
    if metric=='1':
        cb.set_label(r'Sea Surface Temperature Root Mean Square Error [$^\circ\!$C]', fontsize=9, color='0.2',labelpad=0)
    if metric=='2':
        cb.set_label(r'Sea Surface Temperature Mean Absolute Percentage Error [%]', fontsize=9, color='0.2',labelpad=-0.5)
    if metric=='3':
        cb.set_label(r'Sea Surface Temperature Bias [$^\circ\!$C]', fontsize=9, color='0.2',labelpad=-0.5)        
    cb.ax.tick_params(labelsize=9, length=2, color='0.2', labelcolor='0.2',direction='in') 
    cb.set_ticks(ticks)
    try:
        os.makedirs("roms_evaluation")
    except FileExistsError:
        pass 
    if metric=='1' and choose_data == '1':
        plt.savefig('/media/ueslei/Ueslei/Scripts/Python/roms_evaluation/sst_rmse_roms_mur.png', transparent=False, bbox_inches = 'tight', pad_inches=0, dpi=250)          
    if metric=='2' and choose_data =='1':
        plt.savefig('/media/ueslei/Ueslei/Scripts/Python/roms_evaluation/sst_mape_roms_mur.png', transparent=False, bbox_inches = 'tight', pad_inches=0, dpi=250)
    if metric=='3' and choose_data =='1':
        plt.savefig('/media/ueslei/Ueslei/Scripts/Python/roms_evaluation/sst_bias_roms_mur.png', transparent=False, bbox_inches = 'tight', pad_inches=0, dpi=250)          
    if metric=='1' and choose_data == '2':
        plt.savefig('/media/ueslei/Ueslei/Scripts/Python/roms_evaluation/sst_rmse_roms_glorys.png', transparent=False, bbox_inches = 'tight', pad_inches=0, dpi=250)          
    if metric=='2' and choose_data =='2':
        plt.savefig('/media/ueslei/Ueslei/Scripts/Python/roms_evaluation/sst_mape_roms_glorys.png', transparent=False, bbox_inches = 'tight', pad_inches=0, dpi=250)
    if metric=='3' and choose_data =='2':
        plt.savefig('/media/ueslei/Ueslei/Scripts/Python/roms_evaluation/sst_bias_roms_glorys.png', transparent=False, bbox_inches = 'tight', pad_inches=0, dpi=250)     

if contourf_var=='2': 
    print(bg.da_cyan+'Evaluate Sea Surface Currents from: (1) OSCAR or (2) GLORYS.'+bg.rs)
    choose_data  = input()
    if choose_data == '1':
        # Open OSCAR file then calculate its average through time.
        oscar_file    = oscar_dir   
        nc_oscar      = netCDF4.Dataset(oscar_file)
        lon_oscar     = nc_oscar.variables['longitude'][:]
        lon_oscar     = lon_oscar-360
        lat_oscar     = nc_oscar.variables['latitude'][:]
        latli         = np.argmin(np.abs(lat_oscar-latbounds[1]))
        latui         = np.argmin(np.abs(lat_oscar-latbounds[0])) 
        lonli         = np.argmin(np.abs(lon_oscar-lonbounds[0]))
        lonui         = np.argmin(np.abs(lon_oscar-lonbounds[1]))  
        lon_oscar     = nc_oscar.variables['longitude'][lonli:lonui]-360
        lat_oscar     = nc_oscar.variables['latitude'][latli:latui]
        lon_oscar,lat_oscar = np.meshgrid(lon_oscar,lat_oscar)
        oscar_lon_len = len(lon_oscar[0,:])
        oscar_lat_len = len(lat_oscar[:,0])
        oscar_loop    = len(nc_oscar.variables['time'])
        u_oscar2      = np.zeros([oscar_loop,oscar_lat_len,oscar_lon_len])
        v_oscar2      = np.zeros([oscar_loop,oscar_lat_len,oscar_lon_len])   
        bar           = IncrementalBar('Processing Ocean Currents from Oscar', max=oscar_loop)        
        for i in range(0,oscar_loop):
            u_oscar = nc_oscar.variables['u'][i,0,latli:latui,lonli:lonui]
            u_oscar2[i,:,:] = u_oscar
            u_oscar2[u_oscar2==0] = np.nan
            v_oscar = nc_oscar.variables['v'][i,0,latli:latui,lonli:lonui]
            v_oscar2[i,:,:] = v_oscar
            v_oscar2[v_oscar2==0] = np.nan  
            bar.next()
        bar.finish()
    if choose_data == '2':
        glorys_file    = glorys_dir   
        nc_glorys      = netCDF4.Dataset(glorys_file)
        lon_glorys     = nc_glorys.variables['longitude'][:]
        lat_glorys     = nc_glorys.variables['latitude'][:]
        latli          = np.argmin( np.abs( lat_glorys - latbounds[1] ) )
        latui          = np.argmin( np.abs( lat_glorys - latbounds[0] ) ) 
        lonli          = np.argmin( np.abs( lon_glorys - lonbounds[0] ) )
        lonui          = np.argmin( np.abs( lon_glorys - lonbounds[1] ) )  
        lon_glorys     = nc_glorys.variables['longitude'][lonli:lonui]
        lat_glorys     = nc_glorys.variables['latitude'][latui:latli]
        lon_glorys,lat_glorys = np.meshgrid(lon_glorys,lat_glorys)
        glorys_loop    = len(nc_glorys.variables['time'][:])
        glorys_lat_len = len(lat_glorys[:,0])
        glorys_lon_len = len(lon_glorys[0, :])
        u_glorys2      = np.zeros([glorys_loop,glorys_lat_len,glorys_lon_len])
        v_glorys2      = np.zeros([glorys_loop,glorys_lat_len,glorys_lon_len])
        bar            = IncrementalBar('Processing Ocean Currents from Glorys', max=glorys_loop)    
        for i in range(0,glorys_loop):
            u_glorys         = nc_glorys.variables['uo'][i,0,latui:latli,lonli:lonui]
            u_glorys2[i,:,:] = u_glorys
            v_glorys         = nc_glorys.variables['vo'][i,0,latui:latli,lonli:lonui]
            v_glorys2[i,:,:] = v_glorys           
            bar.next()
        bar.finish()

    # Open ROMS file then calculate its average through time.
    if choose_data == '1':
        roms_file = roms_oscar_dir
    if choose_data == '2':
        roms_file = roms_mur_glorys_dir     
    nc_roms      = netCDF4.Dataset(roms_file)
    lon_rho      = nc_roms.variables['lon_rho'][:]
    lat_rho      = nc_roms.variables['lat_rho'][:]
    i0,i1,j0,j1  = bbox2ij(lon_rho,lat_rho,bbox)
    lon_roms     = lon_rho[j0:j1, i0:i1]
    lat_roms     = lat_rho[j0:j1, i0:i1]
    u_roms       = nc_roms.variables['u'][:,0,j0:j1, i0:i1] 
    u_roms       = np.average(u_roms,axis=0)     
    v_roms       = nc_roms.variables['v'][:,0,j0:j1, i0:i1]     
    v_roms       = np.average(v_roms,axis=0) 
    lat_roms     = lat_rho[j0:j1, i0:i1]
    roms_loop    = len(nc_roms.variables['ocean_time'][:])
    if choose_data == '1':
        u_roms2      = np.zeros([roms_loop,oscar_lat_len,oscar_lon_len])
        v_roms2      = np.zeros([roms_loop,oscar_lat_len,oscar_lon_len])       
        orig_def     = pyresample.geometry.SwathDefinition(lons=lon_roms, lats=lat_roms)
        targ_def     = pyresample.geometry.SwathDefinition(lons=lon_oscar, lats=lat_oscar)
    if choose_data == '2': 
        u_roms2      = np.zeros([roms_loop,glorys_lat_len,glorys_lon_len])
        v_roms2      = np.zeros([roms_loop,glorys_lat_len,glorys_lon_len])       
        orig_def     = pyresample.geometry.SwathDefinition(lons=lon_roms, lats=lat_roms)
        targ_def     = pyresample.geometry.SwathDefinition(lons=lon_glorys, lats=lat_glorys)
        bar          = IncrementalBar('Processing Ocean Currents from ROMS', max=roms_loop)    
    for i in range(0,roms_loop):
        u_roms         = nc_roms.variables['u'][i,0,j0:j1,i0:i1]
        u_roms         = pyresample.kd_tree.resample_gauss(orig_def, u_roms, targ_def, radius_of_influence=500000, neighbours=20,sigmas=250000, fill_value=None)
        u_roms2[i,:,:] = u_roms
        v_roms         = nc_roms.variables['v'][i,0,j0:j1,i0:i1]
        v_roms         = pyresample.kd_tree.resample_gauss(orig_def, v_roms, targ_def, radius_of_influence=500000, neighbours=20,sigmas=250000, fill_value=None)
        v_roms2[i,:,:] = v_roms

    # Calculate speed from U and V vectors.
    if choose_data == '1':
        cur_roms  = np.sqrt(u_roms2**2 + v_roms2**2)
        cur_oscar = np.sqrt(u_oscar2**2 + v_oscar2**2)
        print(bg.da_cyan+'Which statistical metric? (1) Root Mean Square Error, (2) Mean Absolute Error or (3) Bias.'+bg.rs)
        metric  = input()
        if metric=='1':
            differences         = cur_roms-cur_oscar
            differences_squared = differences ** 2 
            mean_of_differences_squared = np.average(differences_squared,axis=0)
            val                 = np.sqrt(mean_of_differences_squared)
        if metric=='2':
            val = np.abs((cur_oscar - cur_roms) / cur_oscar).mean(axis=0)*100 
        if metric=='3':
            mean_oscar = np.average(cur_oscar,axis=0) 
            mean_roms  = np.average(cur_roms,axis=0)
            val        = mean_roms-mean_oscar

    if choose_data == '2':
        cur_roms   = np.sqrt(u_roms2**2 + v_roms2**2)
        cur_glorys = np.sqrt(u_glorys2**2 + v_glorys2**2)
        print(bg.da_cyan+'Which statistical metric? (1) Root Mean Square Error, (2) Mean Absolute Percentage Error or (3) Bias.'+bg.rs)
        metric  = input()
        if metric=='1':
            differences         = cur_roms-cur_glorys
            differences_squared = differences ** 2 
            mean_of_differences_squared = np.average(differences_squared,axis=0)
            val                 = np.sqrt(mean_of_differences_squared)
        if metric=='2':
            val = np.abs((cur_glorys - cur_roms) / cur_glorys).mean(axis=0)*100 
        if metric=='3':
            mean_roms   = np.average(cur_roms,axis=0) 
            mean_glorys = np.average(cur_glorys,axis=0)
            val         = mean_roms-mean_glorys
                         
    # Create and plot map.
    m    = Basemap(projection='merc',llcrnrlat=bbox[2],urcrnrlat=bbox[3],llcrnrlon=bbox[0],urcrnrlon=bbox[1], lat_ts=30,resolution='i')
    fig  = plt.figure(1,figsize=(10,8))
    plt.xlabel('Longitude'u' [\N{DEGREE SIGN}]',labelpad=18,size=10)
    plt.ylabel('Latitude'u' [\N{DEGREE SIGN}]',labelpad=33,size=10)
    ax   = fig.add_subplot(111)
    m.drawparallels(np.arange(-90.,120.,1.), linewidth=0.00, color='black', labels=[1,0,0,1],labelstyle="N/S",fontsize=10)
    m.drawmeridians(np.arange(-180.,180.,1.), linewidth=0.00,color='black', labels=[1,0,0,1],labelstyle="N/S",fontsize=10)
    m.drawcountries(color = '#ffffff',linewidth=0.5)
    m.drawcoastlines(color = '#ffffff',linewidth=0.5)
    m.fillcontinents(color = '#000000')
    if metric=='1':
        clevs = clevs_rmse_cur
        ticks = ticks_rmse_cur
        cmap  =  matplotlib.pyplot.jet()   #cmocean.cm.thermal
    if metric=='2':
        clevs = clevs_mape_cur
        ticks = ticks_mape_cur
        cmap  = cmocean.cm.balance
    if metric=='3':
        clevs = clevs_bias_cur
        ticks = ticks_bias_cur
        cmap  = cmocean.cm.balance

    if choose_data == '1':
        if metric == '1' or metric=='2':
            h1    = m.contourf(lon_oscar, lat_oscar, val, clevs,latlon=True,cmap=cmap,extend="both")
        if metric == '3':
            h1    = m.contourf(lon_oscar, lat_oscar, val, clevs,latlon=True,cmap=cmap,norm=MidpointNormalize(midpoint=0),extend="both")          
    if choose_data =='2':
        if metric == '1' or metric=='2':
            h1    = m.contourf(lon_glorys, lat_glorys, val, clevs,latlon=True,cmap=cmap,extend="both")   
        if metric == '3':
            h1    = m.contourf(lon_glorys, lat_glorys, val, clevs,latlon=True,cmap=cmap,norm=MidpointNormalize(midpoint=0),extend="both")              
    cax   = fig.add_axes([0.37, 0.025, 0.27, 0.025])     
    cb    = fig.colorbar(h1, cax=cax, orientation="horizontal",panchor=(0.5,0.5),shrink=0.3,ticks=ticks)
    if metric=='1':
        cb.set_label(r'Surface Currents Root Mean Square Error [m.s⁻¹]', fontsize=10, color='0.2',labelpad=-0.5)
    if metric=='2':
        cb.set_label(r'Current at Surface Mean Absolute Percentage Error [%]', fontsize=10, color='0.2',labelpad=-0.5)
    if metric=='3':
        cb.set_label(r'Current at Surface Bias [m.s⁻¹]', fontsize=10, color='0.2',labelpad=-0.5)      
    cb.ax.tick_params(labelsize=10, length=2, color='0.2', labelcolor='0.2',direction='in') 
    cb.set_ticks(ticks)
    try:
        os.makedirs("roms_evaluation")
    except FileExistsError:
        pass 
    
    if metric=='1' and choose_data == '1':
        plt.savefig('/media/ueslei/Ueslei/Scripts/Python/roms_evaluation/sscs_rmse_roms_oscar.png', transparent=False, bbox_inches = 'tight', pad_inches=0, dpi=250)              
    if metric=='2' and choose_data == '1':
        plt.savefig('/media/ueslei/Ueslei/Scripts/Python/roms_evaluation/sscs_mape_roms_oscar.png', transparent=False, bbox_inches = 'tight', pad_inches=0, dpi=250)  
    if metric=='3' and choose_data == '1':
        plt.savefig('/media/ueslei/Ueslei/Scripts/Python/roms_evaluation/sscs_bias_roms_oscar.png', transparent=False, bbox_inches = 'tight', pad_inches=0, dpi=250)             
    if metric=='1' and choose_data == '2':
        plt.savefig('/media/ueslei/Ueslei/Scripts/Python/roms_evaluation/sscs_rmse_roms_glorys.png', transparent=False, bbox_inches = 'tight', pad_inches=0, dpi=250)              
    if metric=='2' and choose_data == '2':
        plt.savefig('/media/ueslei/Ueslei/Scripts/Python/roms_evaluation/sscs_mape_roms_glorys.png', transparent=False, bbox_inches = 'tight', pad_inches=0, dpi=250)    
    if metric=='3' and choose_data == '2':
        plt.savefig('/media/ueslei/Ueslei/Scripts/Python/roms_evaluation/sscs_bias_roms_glorys.png', transparent=False, bbox_inches = 'tight', pad_inches=0, dpi=250)    



