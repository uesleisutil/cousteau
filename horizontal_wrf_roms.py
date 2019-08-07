#!/home/uesleisutil/anaconda3/bin/python python
# -*- coding: utf-8 -*-

"""
File name:      roms_wrf_horizontal.py
Author:         Ueslei Adriano Sutil
Email:          uesleisutil1@gmail.com
Created:        27 February 2019
Last modified:  07 August 2019
Version:        3.3
Python:         3.7.1

Creates horizontal plots from ROMS (his) and WRF-ARW outputs, displaying:
    - Temperature (Contourf; °C);
    - Salinity (Contourf; PSU);
    - Latent Heat Flux (Contourf; W/m-2);
    - 200 and 1000 meters bathymetry (Contour; m);
    - Wind vectors at 10 m (Vector; m/s);
    - Ocean current at surface (Vector; m/s);
    - Sea Level Pressure (Contourf; hPa);
    - Accumulated Total Precipitation (Contourf; mm);
    - Sea Ice Thickness (Contourf; %).

bbox = [lon_min,lon_max,lat_min,lat_max]
"""

import numpy                  as np
import matplotlib.pyplot      as plt
from   mpl_toolkits.basemap   import Basemap
import netCDF4
from   roms_libs              import *
import cmocean
from   OceanLab.utils         import download_bathy
import os
from   wrf                    import to_np, getvar, latlon_coords, smooth2d, extract_times
from   scipy.ndimage.filters  import gaussian_filter
from   sty                    import bg, rs
from   datetime               import datetime
import pandas                 as pd
from   progress.bar import IncrementalBar
# 2. Customizations.

print(bg.da_cyan+'Which project? (1) SC_2008, (2) ATLEQ or (3) Antartic.'+bg.rs)
project = input()
if project=='1':

    print(bg.da_cyan+'Which experiment? (1) -100%, (2) -80%, (3) -60%, (4) -40%, (5), -20%, (6) Normal, (7) +20%, (8) +40%, (9) +60 (10) +80 or (11) +100%.'+bg.rs)
    exp = input()
    if exp=='1':
        roms_file = '/media/ueslei/Ueslei/INPE/PCI/Projetos/SC_2008/Outputs/cold_sst_100/roms.nc'
        wrf_file  = '/media/ueslei/Ueslei/INPE/PCI/Projetos/SC_2008/Outputs/cold_sst_100/wrf.nc'
    if exp=='2':
        roms_file = '/media/ueslei/Ueslei/INPE/PCI/Projetos/SC_2008/Outputs/cold_sst_80/roms.nc'
        wrf_file  = '/media/ueslei/Ueslei/INPE/PCI/Projetos/SC_2008/Outputs/cold_sst_80/wrf.nc'
    if exp=='3':
        roms_file = '/media/ueslei/Ueslei/INPE/PCI/Projetos/SC_2008/Outputs/cold_sst_60/roms.nc'
        wrf_file  = '/media/ueslei/Ueslei/INPE/PCI/Projetos/SC_2008/Outputs/cold_sst_60/wrf.nc' 
    if exp=='4':
        roms_file = '/media/ueslei/Ueslei/INPE/PCI/Projetos/SC_2008/Outputs/cold_sst_40/roms.nc'
        wrf_file  = '/media/ueslei/Ueslei/INPE/PCI/Projetos/SC_2008/Outputs/cold_sst_40/wrf.nc' 
    if exp=='5':
        roms_file = '/media/ueslei/Ueslei/INPE/PCI/Projetos/SC_2008/Outputs/cold_sst_20/roms.nc'
        wrf_file  = '/media/ueslei/Ueslei/INPE/PCI/Projetos/SC_2008/Outputs/cold_sst_20/wrf.nc' 
    if exp=='6':
        roms_file = '/media/ueslei/Ueslei/INPE/PCI/Projetos/SC_2008/Outputs/normal/roms.nc'
        wrf_file  = '/media/ueslei/Ueslei/INPE/PCI/Projetos/SC_2008/Outputs/normal/wrf.nc' 
    if exp=='7':
        roms_file = '/media/ueslei/Ueslei/INPE/PCI/ProjetosSC_2008/Outputs/warm_20/roms.nc'
        wrf_file  = '/media/ueslei/Ueslei/INPE/PCI/Projetos/SC_2008/Outputs/warm_20/wrf.nc' 
    if exp=='8':
        roms_file = '/media/ueslei/Ueslei/INPE/PCI/Projetos/SC_2008/Outputs/warm_40/roms.nc'
        wrf_file  = '/media/ueslei/Ueslei/INPE/PCI/Projetos/SC_2008/Outputs/warm_40/wrf.nc'
    if exp=='9':
        roms_file = '/media/ueslei/Ueslei/INPE/PCI/Projetos/SC_2008/Outputs/warm_60/roms.nc'
        wrf_file  = '/media/ueslei/Ueslei/INPE/PCI/Projetos/SC_2008/Outputs/warm_60/wrf.nc' 
    if exp=='10':
        roms_file = '/media/ueslei/Ueslei/INPE/PCI/Projetos/SC_2008/Outputs/warm_80/roms.nc'
        wrf_file  = '/media/ueslei/Ueslei/INPE/PCI/Projetos/SC_2008/Outputs/warm_80/wrf.nc' 
    if exp=='11':
        roms_file = '/media/ueslei/Ueslei/INPE/PCI/Projetos/SC_2008/Outputs/warm_100/roms.nc'
        wrf_file  = '/media/ueslei/Ueslei/INPE/PCI/Projetos/SC_2008/Outputs/warm_100/wrf.nc' 
    bbox          = [-56, -42, -34., -22.2]
    initloop      = 144
    zlev          = -1 # Last sigma layer corresponds to surface.
    plot_var      = True 
    plot_bathy    = False
    plot_currents = False
    plot_wind     = False
    plot_slp      = False
    create_video  = False
    nc_roms       = netCDF4.Dataset(roms_file)
    nc_wrf        = netCDF4.Dataset(wrf_file)
    atemp         = nc_roms.variables['temp']
    ntimes        = len(atemp)
if project=='2':
    roms_file     = '/media/ueslei/Ueslei/INPE/PCI/Projetos/ATLEQ/Outputs/roms.nc'
    wrf_file      = '/media/ueslei/Ueslei/INPE/PCI/Projetos/ATLEQ/Outputs/wrf.nc'
    bbox          = [-60,15,-15.1,15.1]
    initloop      = 0
    zlev          = -1 # Last sigma layer corresponds to surface.
    plot_var      = True 
    plot_bathy    = False
    plot_currents = False
    plot_wind     = False
    plot_slp      = False
    create_video  = True
    nc_roms       = netCDF4.Dataset(roms_file)
    nc_wrf        = netCDF4.Dataset(wrf_file)
    atemp         = nc_roms.variables['temp']
    ntimes        = len(atemp)
if project=='3':
    roms_file     = '/media/ueslei/Ueslei/INPE/PCI/Projetos/Antartic/Outputs/roms.nc'
    wrf_file      = '/media/ueslei/Ueslei/INPE/PCI/Projetos/Antartic/Outputs/wrf.nc'
    bbox          = [-179., 179, -79., -45.]
    initloop      = 0
    zlev          = -1 # Last sigma layer corresponds to surface.
    plot_var      = True 
    plot_bathy    = False
    plot_currents = False
    plot_wind     = False
    plot_slp      = False
    create_video  = True
    nc_roms       = netCDF4.Dataset(roms_file)
    nc_wrf        = netCDF4.Dataset(wrf_file)
    atemp         = nc_roms.variables['temp']
    ntimes        = len(atemp)

if project =='1' or project=='2':
    print(bg.da_cyan+'Which contour? (1) Sea Surface Temperature, (2) Salinity, (3) Precipitation or (4) Heat Fluxes.'+bg.rs)
    contourf_var  = input()
if project =='3':
    print(bg.da_cyan+'Which contour? (1) Sea Surface Temperature, (2) Salinity, (3) Precipitation, (4) Heat Fluxes or (5) Fraction of cell covered by ice.'+bg.rs)
    contourf_var  = input() 

# 3. Start looping through time
bar = IncrementalBar('', max=ntimes)
for i in range(0,ntimes):
#for i in range(0,1):
    # 3.1. Get variables and their resources.
    #tvar        = nc_roms.variables['ocean_time']
    #timestr     = netCDF4.num2date(tvar[i], tvar.units).strftime('%b %d, %Y %H:%M')
    timestr1 = extract_times(nc_wrf,timeidx=None,meta=False,do_xtime=False)
    timestr11 = timestr1[i]
    timestr  = pd.to_datetime(timestr11, format="%b %d %Y %H:%M")

    if plot_wind==True:
        uvmet10 = getvar(nc_wrf, "uvmet10", i, units="m s-1")
        lats_wrf, lons_wrf = latlon_coords(uvmet10)
    if plot_slp==True:
        slp         = getvar(nc_wrf, "slp", i, units="hPa")
        slp         = smooth2d(slp,30)
        lats_wrf, lons_wrf = latlon_coords(slp)
        clevs_slp   = np.arange(900,1030,3)
    if plot_currents==True:
        lon_rho     = nc_roms.variables['lon_rho'][:]
        lat_rho     = nc_roms.variables['lat_rho'][:]
        i0,i1,j0,j1 = bbox2ij(lon_rho,lat_rho,bbox)
        lon_var     = lon_rho[j0:j1, i0:i1]
        lat_var     = lat_rho[j0:j1, i0:i1]
        lon         = lon_rho[(j0+1):(j1-1), (i0+1):(i1-1)]
        lat         = lat_rho[(j0+1):(j1-1), (i0+1):(i1-1)]
        u           = nc_roms.variables['u'][i, zlev, j0:j1, i0:(i1-1)]
        v           = nc_roms.variables['v'][i, zlev, j0:(j1-1), i0:i1]
        mask        = 1 - nc_roms.variables['mask_rho'][(j0+1):(j1-1), (i0+1):(i1-1)]
        ang         = nc_roms.variables['angle'][(j0+1):(j1-1), (i0+1):(i1-1)] 
        u           = shrink(u, mask.shape) #Average U to rho points.
        v           = shrink(v, mask.shape) # Average V points to rho points.
        u, v        = rot2d(u, v, ang) # Rotate grid oriented U and V to East and West.
    if plot_var==True:
        if project=='1' and contourf_var=='1':
            lon_rho     = nc_roms.variables['lon_rho'][:]
            lat_rho     = nc_roms.variables['lat_rho'][:]
            i0,i1,j0,j1 = bbox2ij(lon_rho,lat_rho,bbox)
            lon_var     = lon_rho[j0:j1, i0:i1]
            lat_var     = lat_rho[j0:j1, i0:i1]
            var         = nc_roms.variables['temp'][i, zlev,  j0:j1, i0:i1]
            clevs       = np.arange(18,26.1,0.01)
            ticks       = np.arange(min(clevs),max(clevs),2)  
            cmap        = cmocean.cm.balance    # matplotlib.pyplot.jet() 
        if project=='2' and contourf_var=='2':
            lon_rho     = nc_roms.variables['lon_rho'][:]
            lat_rho     = nc_roms.variables['lat_rho'][:]
            i0,i1,j0,j1 = bbox2ij(lon_rho,lat_rho,bbox)
            lon_var     = lon_rho[j0:j1, i0:i1]
            lat_var     = lat_rho[j0:j1, i0:i1]
            var        = nc_roms.variables['temp'][i, zlev,  j0:j1, i0:i1]
            clevs       = np.arange(19,31.1,0.01)
            ticks       = np.arange(min(clevs),max(clevs),1)  
            cmap        = matplotlib.pyplot.jet()  
        if project=='3' and contourf_var=='1':
            lon_rho     = nc_roms.variables['lon_rho'][:]
            lat_rho     = nc_roms.variables['lat_rho'][:]
            i0,i1,j0,j1 = bbox2ij(lon_rho,lat_rho,bbox)
            lon_var     = lon_rho[j0:j1, i0:i1]
            lat_var     = lat_rho[j0:j1, i0:i1]
            var         = nc_roms.variables['temp'][i, zlev,  j0:j1, i0:i1]
            clevs       = np.arange(-4,8.1,0.01)
            ticks       = np.arange(min(clevs),max(clevs),2)  
            cmap        = matplotlib.pyplot.jet()     
        if contourf_var=='2':
            i0,i1,j0,j1 = bbox2ij(lon_rho,lat_rho,bbox)
            lon_rho = nc_roms.variables['lon_rho'][:]
            lat_rho = nc_roms.variables['lat_rho'][:]
            var     = nc_roms.variables['salt'][i, zlev,  j0:j1, i0:i1]
            lon_var = lon_rho[j0:j1, i0:i1]
            lat_var = lat_rho[j0:j1, i0:i1]
            clevs   = np.arange(31,38.01,0.05)
            cmap    = matplotlib.pyplot.viridis()  
            ticks   = np.array([31,32,33,34,35,36,37,38])
        if contourf_var=='3':
            rainnc   = nc_wrf.variables['RAINNC'][i,:,:]
            rainc    = nc_wrf.variables['RAINC'][i,:,:]
            lat_rain = nc_wrf.variables['XLAT'][i,:,0]
            lon_rain = nc_wrf.variables['XLONG'][i,0,:]
            var      = rainnc+rainc    
            clevs    = np.arange(0,501,5)
            ticks    = np.array([0,100,200,300,400,500])  
            cmap     = matplotlib.pyplot.jet()  #cmocean.cm.rain    
        if contourf_var=='4': 
            lon_rho     = nc_roms.variables['lon_rho'][:]
            lat_rho     = nc_roms.variables['lat_rho'][:]
            i0,i1,j0,j1 = bbox2ij(lon_rho,lat_rho,bbox)
            lon_var     = lon_rho[j0:j1, i0:i1]
            lat_var     = lat_rho[j0:j1, i0:i1]
            var1        = nc_roms.variables['latent'][i, j0:j1, i0:i1]
            var2        = nc_roms.variables['sensible'][i, j0:j1, i0:i1]
            var         = abs(var1+var2) # Change to positive flux (atmos to ocean)
            clevs       = np.arange(0,600.1,5)
            ticks       = np.array([0,100,200,300,400,500,600,700,800,900])   
            cmap        = matplotlib.pyplot.jet()   
        if contourf_var=='5':
            lon_rho     = nc_roms.variables['lon_rho'][:]
            lat_rho     = nc_roms.variables['lat_rho'][:]
            i0,i1,j0,j1 = bbox2ij(lon_rho,lat_rho,bbox)
            lon_var     = lon_rho[j0:j1, i0:i1]
            lat_var     = lat_rho[j0:j1, i0:i1]
            var         = nc_roms.variables['aice'][i, j0:j1, i0:i1]*100         
            clevs       = np.arange(0,101,0.1)
            ticks       = np.arange(min(clevs),max(clevs),20) 
            cmap        = cmocean.cm.ice_r
    # 3.2. Create a figure.
    if project=='1' or project=='2':
        m = Basemap(projection='merc',llcrnrlat=bbox[2],urcrnrlat=bbox[3],llcrnrlon=bbox[0],urcrnrlon=bbox[1], lat_ts=30,resolution='i')
    if project=='3':
        m = Basemap(width=4700000,height=3300000,
            rsphere=(6378137.00,6356752.3142),\
            resolution='l',area_thresh=500.,projection='lcc',\
            lat_1=-45.,lat_2=-65,lat_0=-65,lon_0=-35.)     
    fig  = plt.figure(1,figsize=(10,8))
    if project=='1' or project=='2':
        plt.xlabel('Longitude'u' [\N{DEGREE SIGN}]',labelpad=11,size=6)
        plt.ylabel('Latitude'u' [\N{DEGREE SIGN}]',labelpad=19,size=6)
    ax   = fig.add_subplot(111)
    if project=='1' or project=='2':
        plt.title(timestr, fontsize=6)
    if project=='3':
        plt.title(timestr, fontsize=9, pad=25)
        
    # 3.3. Add coastline, continents and lat/lon.
    if project=='1':
        m.drawparallels(np.arange(-90.,120.,1.), linewidth=0.00, color='black', labels=[1,0,0,1],labelstyle="N/S",fontsize=6)
        m.drawmeridians(np.arange(-180.,180.,2.), linewidth=0.00,color='black', labels=[1,0,0,1],labelstyle="N/S",fontsize=6)
    if project=='2':
        m.drawparallels(np.arange(-90.,120.,5), linewidth=0.00, color='black', labels=[1,0,0,1],labelstyle="N/S",fontsize=6)
        m.drawmeridians(np.arange(-180.,180.,5), linewidth=0.00,color='black', labels=[1,0,0,1],labelstyle="N/S",fontsize=6)
    if project=='3':
        m.drawmeridians(np.arange(0,360,10),labels=[1,1,1,0],linewidth=0.1,fontsize=9,labelstyle="N/S")
        m.drawparallels(np.arange(-90,90,5),linewidth=0.1,fontsize=9,labelstyle="N/S")
    if contourf_var=='3':
        m.drawcountries(color = '#000000')
        m.drawcoastlines(color = '#000000')
    else:
        m.fillcontinents(color = '#000000')
        m.drawcountries(color = '#ffffff',linewidth=0.5)
        m.drawcoastlines(color = '#ffffff',linewidth=0.5)

    # 3.4. Plot Current vector.
    if plot_currents==True:
        x_rho, y_rho = m(lon,lat)
        if project=='1':
            nsub  = 22
            scale = 0.045
        if project=='2':
            nsub  = 20
            scale = 0.03
        if project=='3':
            nsub  = 20
            scale = 0.07          
        C = ax.quiver(x_rho[::nsub,::nsub],y_rho[::nsub,::nsub],u[::nsub,::nsub],v[::nsub,::nsub],alpha=0.5,scale=1/scale, zorder=1e35, width=0.0025,color='black',pivot='middle')
        if project=='1':
            qk = ax.quiverkey(C, .16, -0.125, 0.5, ' Sea Surface Current\n 0.5 m.s⁻¹ ', coordinates='axes',color='black',labelsep=0.05, labelcolor='black',alpha=0.4,fontproperties={'size': '5'})
        if project=='2':
            qk = ax.quiverkey(C, .22, -0.25, 0.5, ' Sea Surface Current\n 0.5 m.s⁻¹ ', coordinates='axes',color='black',labelsep=0.05, labelcolor='black',alpha=0.4,fontproperties={'size': '5'})
        if project=='3':
            qk = ax.quiverkey(C, .2, -0.11, 0.5, ' Sea Surface Current\n 0.5 m.s⁻¹ ', coordinates='axes',color='black',labelsep=0.05, labelcolor='black',alpha=0.4,fontproperties={'size': '6'})
    # 3.5. Plot wind speed at 10 meters.
    if plot_wind==True:
        x, y    = m(to_np(lons_wrf), to_np(lats_wrf))
        if project=='1':
            spacing = 22
            scale   = 0.03
        if project=='2':
            spacing = 40
            scale   = 0.03
        if project=='3':
            spacing = 25
            scale   = 0.02
        W = ax.quiver(x[::spacing,::spacing], y[::spacing,::spacing], to_np(uvmet10[0,::spacing, ::spacing]),to_np(uvmet10[1,::spacing, ::spacing]),pivot='middle',scale=8/scale, zorder=1e35, width=0.005,color='gray',headlength=3, headaxislength=2.8 )
        if project=='1':
            wk = ax.quiverkey(W, 0.81, -0.125, 10, ' Wind Vector at 10 m\n 10 m.s⁻¹ ', coordinates='axes',color='#444444',labelsep=0.05, labelcolor='black',fontproperties={'size': '5'})
        if project=='2':
            wk = ax.quiverkey(W, .75, -.25, 10, ' Wind Vector at 10 m\n 10 m.s⁻¹ ', coordinates='axes',color='#444444',labelsep=0.05, labelcolor='black',fontproperties={'size': '5'})
        if project=='3':
            wk = ax.quiverkey(W, .75, -.25, 10, ' Wind Vector at 10 m\n 10 m.s⁻¹ ', coordinates='axes',color='#444444',labelsep=0.05, labelcolor='black',fontproperties={'size': '5'})
        
    # 3.6. Plot the bathymetry then apply a gaussian filter to smooth it.
    if plot_bathy==True:
        bathy_levels     = [-1000,-200]
        bLON,bLAT,BAT    = download_bathy(lnd=bbox[0],lnu=-39,ltd=bbox[2],ltu=bbox[3])
        Ct               = m.contour(gaussian_filter(bLON,3),gaussian_filter(bLAT,3),gaussian_filter(BAT,3),bathy_levels,colors='black',latlon=True,linewidths=0.3,linestyles='solid')
        manual_locations = [(554743,64375), (517235,81424)]
        clbls            = plt.clabel(Ct,fmt='%i', fontsize=4,manual=manual_locations)

    # 3.7. Plot the desired countour variable and add some plot resources.
    if plot_slp==True:
        x_slp, y_slp = m(to_np(lons_wrf), to_np(lats_wrf))
        Ct_slp       = ax.contour(gaussian_filter(x_slp,0),gaussian_filter(y_slp,0),gaussian_filter(to_np(slp),0),clevs_slp,colors='white',latlon=True,linewidths=0.9,linestyles='solid')
        clbls_slp    = plt.clabel(Ct_slp,fmt='%i',inline=1, fontsize=5)

    if plot_var==True:
        if contourf_var=='1':
            h1  = m.contourf(lon_var, lat_var, var, clevs,latlon=True,cmap=cmap,extend="both")  
            if project=='1':
                cax = fig.add_axes([0.37, 0.03, 0.27, 0.025])   
            if project=='2':
                cax = fig.add_axes([0.37, 0.196, 0.27, 0.025])       
            if project=='3':
                cax = fig.add_axes([0.37, 0.08, 0.27, 0.025])                  
            cb  = fig.colorbar(h1, cax=cax, orientation="horizontal",panchor=(0.5,0.5),shrink=0.3,ticks=ticks)
            cb.set_label(r'Sea Surface Temperature [$^\circ\!$C]', fontsize=7, color='0.2',labelpad=-0.5)
            cb.ax.tick_params(labelsize=7, length=2, color='0.2', labelcolor='0.2',direction='in') 
            cb.set_ticks(ticks)

        if contourf_var=='2':
            h1  = m.contourf(lon_var, lat_var, var, clevs,latlon=True,cmap=cmap)  
            cax = fig.add_axes([0.37, 0.19, 0.27, 0.025])     
            cb  = fig.colorbar(h1, cax=cax, orientation="horizontal",panchor=(0.5,0.5),shrink=0.3,ticks=ticks)
            cb.set_label(r'Sea Surface Salinity [PSU]', fontsize=5, color='0.2',labelpad=-0.5)
            cb.ax.tick_params(labelsize=5, length=2, color='0.2', labelcolor='0.2',direction='in') 
            cb.set_ticks(ticks)
        if contourf_var=='3':
            h1  = ax.contourf(lon_rain, lat_rain, var, clevs,latlon=True,cmap=cmap)  
            cax = fig.add_axes([0.37, 0.19, 0.27, 0.025])     
            cb  = fig.colorbar(h1, cax=cax, orientation="horizontal",panchor=(0.5,0.5),shrink=0.3,ticks=ticks)
            cb.set_label(r'Accumulated precipitation [mm]', fontsize=5, color='0.2',labelpad=-0.5)
            cb.ax.tick_params(labelsize=5, length=2, color='0.2', labelcolor='0.2',direction='in') 
            cb.set_ticks(ticks)  
        if contourf_var=='4':
            h1  = m.contourf(lon_var, lat_var, var, clevs,latlon=True,cmap=cmap)  
            cax = fig.add_axes([0.37, 0.19, 0.27, 0.025])     
            cb  = fig.colorbar(h1, cax=cax, orientation="horizontal",panchor=(0.5,0.5),shrink=0.3,ticks=ticks)
            cb.set_label(r'Heat Fluxes [W m.s⁻²]', fontsize=5, color='0.2',labelpad=-0.5)
            cb.ax.tick_params(labelsize=5, length=2, color='0.2', labelcolor='0.2',direction='in') 
            cb.set_ticks(ticks)
        if contourf_var=='5':
            h1  = m.contourf(lon_var, lat_var, var, clevs,latlon=True,cmap=cmap,extend="both")  
            cax = fig.add_axes([0.37, 0.11, 0.27, 0.025]) 
            cb  = fig.colorbar(h1, cax=cax, orientation="horizontal",panchor=(0.5,0.5),shrink=0.3,ticks=ticks)
            cb.set_label(r'Fractional sea ice cover [%]', fontsize=12, color='0.2',labelpad=-0.5)
            cb.ax.tick_params(labelsize=12, length=2, color='0.2', labelcolor='0.2',direction='in') 
            cb.set_ticks(ticks)  

    # 3.8. Save figures.
    if contourf_var=='1':
        try:
            os.makedirs("sst")
        except FileExistsError:
            pass
        plt.savefig('./sst/temp_{0:03d}.png'.format(i), transparent=True, bbox_inches = 'tight', pad_inches=0, dpi=250)
        plt.clf()
    elif contourf_var=='2':
        try:
            os.makedirs("sss")
        except FileExistsError:
            pass
        plt.savefig('./sss/salt_{0:03d}.png'.format(i), transparent=False, bbox_inches = 'tight', pad_inches=0, dpi=250)
        plt.clf()   
    elif contourf_var=='3':
        try:
            os.makedirs("precipitation")
        except FileExistsError:
            pass        
        plt.savefig('./precipitation/prec_{0:03d}.png'.format(i), transparent=False, bbox_inches = 'tight', pad_inches=0, dpi=250)
        plt.clf()
    elif contourf_var=='4':
        try:
            os.makedirs("heat_fluxes")
        except FileExistsError:
            pass
        plt.savefig('./heat_fluxes/fluxes_{0:03d}.png'.format(i), transparent=False, bbox_inches = 'tight', pad_inches=0, dpi=250)
        plt.clf()
    elif contourf_var=='5':
        try:
            os.makedirs("ice")
        except FileExistsError:
            pass
        plt.savefig('./ice/ice_{0:03d}.png'.format(i), transparent=False, bbox_inches = 'tight', pad_inches=0, dpi=250)
        plt.clf()   
    bar.next()
bar.finish()
# 4. Create mp4 file from figures.
if create_video==True:
    cwd = os.getcwd()
    if contourf_var=='1':
        exists = os.path.isfile('./sst/wrf_roms_temp.mp4')
        if exists==True:
            os.system("rm -rf ./sst/wrf_roms_temp.mp4")
            os.system("ffmpeg -r 10 -pattern_type glob -i '"+cwd+"/sst/*.png' -c:v libx264 -pix_fmt yuv420p -crf 15 ./sst/wrf_roms_temp.mp4")
        else:
            os.system("ffmpeg -r 10 -pattern_type glob -i '"+cwd+"/sst/*.png' -c:v libx264 -pix_fmt yuv420p -crf 15 ./sst/wrf_roms_temp.mp4")
            print(bg.da_cyan+'Wish to delete the .png files? (1) Yes or (2) No.'+bg.rs)
        removefiles = input()
        if removefiles=='1':
            os.system("rm -rf ./sst/*.png")
        if removefiles=='2':
            pass
    if contourf_var=='2':
        exists = os.path.isfile('./sss/wrf_roms_salt.mp4')
        if exists==True:
            os.system("rm -rf ./salinity/wrf_roms_salt.mp4")
            os.system("ffmpeg -r 10 -pattern_type glob -i '"+cwd+"/sss/*.png' -c:v libx264 -pix_fmt yuv420p -crf 15 ./sss/wrf_roms_salt.mp4")
        else:
            os.system("ffmpeg -r 10 -pattern_type glob -i '"+cwd+"/sss/*.png' -c:v libx264 -pix_fmt yuv420p -crf 15 ./sss/wrf_roms_salt.mp4")
            print(bg.da_cyan+'Wish to delete the .png files? (1) Yes or (2) No.'+bg.rs)
        removefiles = input()
        if removefiles=='1':
            os.system("rm -rf ./sss/*.png")
        if removefiles=='2':
            pass
    if contourf_var=='3':
        exists = os.path.isfile('./precipitation/wrf_roms_prec.mp4')
        if exists==True:
            os.system("rm -rf ./precipitation/wrf_roms_prec.mp4")
            os.system("ffmpeg -r 10 -pattern_type glob -i '"+cwd+"/precipitation/*.png' -c:v libx264 -pix_fmt yuv420p -crf 15 ./precipitation/wrf_roms_prec.mp4")
        else:
            os.system("ffmpeg -r 10 -pattern_type glob -i '"+cwd+"/precipitation/*.png' -c:v libx264 -pix_fmt yuv420p -crf 15 ./precipitation/wrf_roms_prec.mp4")
            print(bg.da_cyan+'Wish to delete the .png files? (1) Yes or (2) No.'+bg.rs)
        removefiles = input()
        if removefiles=='1':
            os.system("rm -rf ./precipitation/*.png")
        if removefiles=='2':
            pass
    if contourf_var=='4':
        exists = os.path.isfile('./heat_fluxes/wrf_roms_flux.mp4')
        if exists==True:
            os.system("rm -rf ./heat_fluxes/wrf_roms_flux.mp4")
            os.system("ffmpeg -r 10 -pattern_type glob -i '"+cwd+"/heat_fluxes/*.png' -c:v libx264 -pix_fmt yuv420p -crf 15 ./heat_fluxes/wrf_roms_flux.mp4")
        else:
            os.system("ffmpeg -r 10 -pattern_type glob -i '"+cwd+"/heat_fluxes/*.png' -c:v libx264 -pix_fmt yuv420p -crf 15 ./heat_fluxes/wrf_roms_flux.mp4")
            print(bg.da_cyan+'Wish to delete the .png files? (1) Yes or (2) No.'+bg.rs)
        removefiles = input()
        if removefiles=='1':
            os.system("rm -rf ./heat_fluxes/*.png")
        if removefiles=='2':
            pass
    if contourf_var=='5':
        exists = os.path.isfile('./ice/wrf_roms_flux.mp4')
        if exists==True:
            os.system("rm -rf ./ice/wrf_roms_flux.mp4")
            os.system("ffmpeg -r 10 -pattern_type glob -i '"+cwd+"/ice/*.png' -c:v libx264 -pix_fmt yuv420p -crf 15 ./ice/wrf_roms_flux.mp4")
        else:
            os.system("ffmpeg -r 10 -pattern_type glob -i '"+cwd+"/ice/*.png' -c:v libx264 -pix_fmt yuv420p -crf 15 ./ice/wrf_roms_flux.mp4")
            print(bg.da_cyan+'Wish to delete the .png files? (1) Yes or (2) No.'+bg.rs)
        removefiles = input()
        if removefiles=='1':
            os.system("rm -rf ./ice/*.png")
        if removefiles=='2':
            pass
else:
    pass
