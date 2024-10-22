#!/home/uesleisutil/anaconda3/bin/python python
# -*- coding: utf-8 -*-

"""
File name:      roms_wrf_horizontal.py
Author:         Ueslei Adriano Sutil
Email:          uesleisutil1@gmail.com
Created:        27 February 2019
Last modified:  13 August 2019
Version:        3.5
Python:         3.7.1

Create horizontal plots from ROMS (his) and WRF-ARW outputs, displaying:
    - Sea surface temperature (Contourf; °C);
    - Sea surface salinity (Contourf; PSU);
    - Heat fluxes (Contourf; W m-2);
    - 200 and 1000 meters bathymetry (Contour; m);
    - Wind vectors at 10 m (Vector; m s-1);
    - Ocean current at surface (Vector; m s-1);
    - Sea level pressure (Contourf; hPa);
    - Accumulated total precipitation (Contourf; mm);
    - Sea ice thickness (Contourf; %);
    - Convective Avaliable Potential Energy (Contourf; J kg-1).

CAPE = Defined as the accumulated negative buoyant energy from the parcel startin point to the 
       level of free convection. The word "parcel" refers to a 500 meter deep parcel, with the
       actual temperature and moisture averaged over that depth.
bbox = [lon_min,lon_max,lat_min,lat_max]
"""

import numpy                 as     np
import matplotlib.pyplot     as     plt
from   mpl_toolkits.basemap  import Basemap
import netCDF4
from   roms_libs             import *
import cmocean
from   OceanLab.utils        import download_bathy
import os
from   wrf                   import to_np, getvar, latlon_coords, smooth2d, extract_times
from   scipy.ndimage.filters import gaussian_filter
from   sty                   import bg, rs
from   datetime              import datetime
import pandas                as     pd
from   progress.bar          import IncrementalBar
matplotlib.use('Agg')

# 2. Customizations.

print(bg.da_cyan+'Which project? (1) SC_2008, (2) ATLEQ or (3) Antartic.'+bg.rs)
project = input()
if project=='1':
    roms_file = '/media/ueslei/Ueslei/INPE/PCI/Projetos/SC_2008/Outputs/normal/roms.nc'
    wrf_file  = '/media/ueslei/Ueslei/INPE/PCI/Projetos/SC_2008/Outputs/normal/wrf.nc'
    bbox          = [-53, -40, -32, -23]
    initloop      = 144
    zlev          = -1 # Last sigma layer corresponds to surface.
    plot_var      = True 
    plot_bathy    = True
    plot_currents = False
    plot_wind     = True
    plot_slp      = False
    create_video  = False
    ppt_fig       = False # If True, it will generate a figure with transparent background.
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
    plot_wind     = True
    plot_slp      = False
    create_video  = True
    ppt_fig       = False
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
    ppt_fig       = False
    nc_roms       = netCDF4.Dataset(roms_file)
    nc_wrf        = netCDF4.Dataset(wrf_file)
    atemp         = nc_roms.variables['temp']
    ntimes        = len(atemp)

if project =='1' or project=='2':
    print(bg.da_cyan+'Which contour? (1) Sea Surface Temperature, (2) Salinity, (3) Precipitation, (4) Heat Fluxes.'+bg.rs)
    contourf_var  = input()
if project =='3':
    print(bg.da_cyan+'Which contour? (1) Sea Surface Temperature, (2) Salinity, (3) Precipitation, (4) Heat Fluxes or (5) Fraction of cell covered by ice.'+bg.rs)
    contourf_var  = input() 

# 3. Start looping through time
range_loop = [i for i in range(168,ntimes,1)]
bar        = IncrementalBar(bg.da_cyan+'Creating figures:'+bg.rs, max=len(range_loop))
for i in range_loop:
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
        slp         = smooth2d(slp,120)
        lats_wrf, lons_wrf = latlon_coords(slp)
        clevs_slp   = np.arange(900,1030,1)
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
        if  contourf_var=='1':
            lon_rho     = nc_roms.variables['lon_rho'][:]
            lat_rho     = nc_roms.variables['lat_rho'][:]
            i0,i1,j0,j1 = bbox2ij(lon_rho,lat_rho,bbox)
            lon_var     = lon_rho[j0:j1, i0:i1]
            lat_var     = lat_rho[j0:j1, i0:i1]
            var         = nc_roms.variables['temp'][i, zlev,  j0:j1, i0:i1]
            clevs       = np.arange(16,25.1,0.01)
            ticks       = np.arange(min(clevs),max(clevs),1)  
            cmap        = cmocean.cm.thermal   
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
            #rainnc0  = nc_wrf.variables['RAINNC'][i-1,:,:]
            #rainc0   = nc_wrf.variables['RAINC'][i-1,:,:]            
            #rainnc1  = nc_wrf.variables['RAINNC'][i,:,:]
            #rainc1   = nc_wrf.variables['RAINC'][i,:,:]
            #rainnc   = rainnc1-rainnc0
            #rainc    = rainc1-rainc0
            rainnc   = nc_wrf.variables['RAINNC'][i,:,:]
            rainc    = nc_wrf.variables['RAINC'][i,:,:]                
            var      = rainnc+rainc
            var1     = getvar(nc_wrf, "pw", i)
            lat_wrf, lon_wrf = latlon_coords(var1)      
            clevs    = np.arange(0,400,5)
            ticks    = np.array([0,50,100,150,200,250,300,350,400])  
            cmap     = cmocean.cm.rain    
        if contourf_var=='4': 
            latent   = nc_wrf.variables['LH'][i,:,:]
            sensible = nc_wrf.variables['HFX'][i,:,:]                
            var      = latent+sensible
            var1     = getvar(nc_wrf, "pw", i)
            lat_wrf, lon_wrf = latlon_coords(var1)      
            clevs    = np.arange(0,401,5)
            ticks    = np.array([0,50,100,150,200,250,300,350,400])  
            cmap     = cmocean.cm.thermal    
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
        plt.xlabel('Longitude'u' [\N{DEGREE SIGN}]',labelpad=18,size=10)
        plt.ylabel('Latitude'u' [\N{DEGREE SIGN}]',labelpad=33,size=10)
    ax   = fig.add_subplot(111)
    if project=='1' or project=='2':
        plt.title(timestr, fontsize=11)
    if project=='3':
        plt.title(timestr, fontsize=10, pad=25)
        
    # 3.3. Add coastline, continents and lat/lon.
    if project=='1':
        m.drawparallels(np.arange(-90.,120.,1.), linewidth=0.00, color='black', labels=[1,0,0,1],labelstyle="N/S",fontsize=10)
        m.drawmeridians(np.arange(-180.,180.,2.), linewidth=0.00,color='black', labels=[1,0,0,1],labelstyle="N/S",fontsize=10)
    if project=='2':
        m.drawparallels(np.arange(-90.,120.,5), linewidth=0.00, color='black', labels=[1,0,0,1],labelstyle="N/S",fontsize=6)
        m.drawmeridians(np.arange(-180.,180.,5), linewidth=0.00,color='black', labels=[1,0,0,1],labelstyle="N/S",fontsize=6)
    if project=='3':
        m.drawmeridians(np.arange(0,360,10),labels=[1,1,1,0],linewidth=0.1,fontsize=8.75,labelstyle="N/S")
        m.drawparallels(np.arange(-90,90,5),linewidth=0.1,fontsize=8.75,labelstyle="N/S")
    if contourf_var=='3':
        m.drawcountries(color = '#000000')
        m.drawcoastlines(color = '#000000')
    else:
        m.fillcontinents(color = '#ffffff')
        m.drawcountries(color = '#000000',linewidth=0.5)
        m.drawcoastlines(color = '#000000',linewidth=0.5)
        m.drawstates(color = '#000000',linewidth=0.5)       

    # 3.4. Plot Current vector.
    if plot_currents==True:
        x_rho, y_rho = m(lon,lat)
        if project=='1':
            nsub  = 20
            scale = 0.065
        if project=='2':
            nsub  = 20
            scale = 0.03
        if project=='3':
            nsub  = 20
            scale = 0.07          
        C = ax.quiver(x_rho[::nsub,::nsub],y_rho[::nsub,::nsub],u[::nsub,::nsub],v[::nsub,::nsub],alpha=0.5,scale=1/scale, zorder=1e35, width=0.0025,color='black',pivot='middle')
        if project=='1':
            qk = ax.quiverkey(C, .215, -0.15, 0.5, ' Sea Surface Current\n 0.5 m.s⁻¹ ', coordinates='axes',color='#444444',labelsep=0.05, labelcolor='black',alpha=1,fontproperties={'size': '9'})
        if project=='2':
            qk = ax.quiverkey(C, .22, -0.25, 0.5, ' Sea Surface Current\n 0.5 m.s⁻¹ ', coordinates='axes',color='black',labelsep=0.05, labelcolor='black',alpha=0.4,fontproperties={'size': '5'})
        if project=='3':
            qk = ax.quiverkey(C, .2, -0.11, 0.5, ' Sea Surface Current\n 0.5 m.s⁻¹ ', coordinates='axes',color='black',labelsep=0.05, labelcolor='black',alpha=0.4,fontproperties={'size': '6'})
    
    # 3.5. Plot wind speed at 10 meters.
    if plot_wind==True:
        x, y    = m(to_np(lons_wrf), to_np(lats_wrf))
        if project=='1':
            spacing = 7
            scale   = 0.03
        if project=='2':
            spacing = 40
            scale   = 0.03
        if project=='3':
            spacing = 25
            scale   = 0.02
        W = ax.quiver(x[::spacing,::spacing], y[::spacing,::spacing], to_np(uvmet10[0,::spacing, ::spacing]),to_np(uvmet10[1,::spacing, ::spacing]),pivot='middle',scale=8/scale, zorder=1e35, width=0.005,color='gray',headlength=3, headaxislength=2.8 )
        if project=='1':
            wk = ax.quiverkey(W, 0.76, -0.15, 10, 'Wind Vector at 10 m\n 10 m.s⁻¹ ', coordinates='axes',color='black',labelsep=0.05, alpha=0.5,labelcolor='black',fontproperties={'size': '9'})
        if project=='2':
            wk = ax.quiverkey(W, .75, -.25, 10, ' Wind Vector at 10 m\n 10 m.s⁻¹ ', coordinates='axes',color='#444444',labelsep=0.05, labelcolor='black',fontproperties={'size': '5'})
        if project=='3':
            wk = ax.quiverkey(W, .75, -.25, 10, ' Wind Vector at 10 m\n 10 m.s⁻¹ ', coordinates='axes',color='#444444',labelsep=0.05, labelcolor='black',fontproperties={'size': '5'})
        
    # 3.6. Plot the bathymetry then apply a gaussian filter to smooth it.
    if plot_bathy==True:
        bathy_levels     = [-1000,-200]
        bLON,bLAT,BAT    = download_bathy(lnd=bbox[0],lnu=-39,ltd=bbox[2],ltu=bbox[3])
        Ct               = m.contour(gaussian_filter(bLON,3),gaussian_filter(bLAT,3),gaussian_filter(BAT,3),bathy_levels,colors='black',latlon=True,linewidths=0.3,linestyles='solid')
        manual_locations = [(411291,138582), (304620,57927)]
        clbls            = plt.clabel(Ct,fmt='%i', fontsize=9,manual=manual_locations)

    # 3.7. Plot the desired countour variable and add some plot resources.
    if plot_slp==True:
        x_slp, y_slp = m(to_np(lons_wrf), to_np(lats_wrf))
        Ct_slp       = ax.contour(gaussian_filter(x_slp,10),gaussian_filter(y_slp,10),gaussian_filter(to_np(slp),10),clevs_slp,colors='white',latlon=True,linewidths=1,linestyles='solid')
        clbls_slp    = plt.clabel(Ct_slp,fmt='%i',inline=1, fontsize=9)

    if plot_var==True:
        if contourf_var=='1':
            h1  = m.contourf(lon_var, lat_var, var, clevs,latlon=True,cmap=cmap,extend="both")  
            if project=='1':
                cax = fig.add_axes([0.37, 0.01, 0.27, 0.025])   
            if project=='2':
                cax = fig.add_axes([0.37, 0.196, 0.27, 0.025])       
            if project=='3':
                cax = fig.add_axes([0.37, 0.08, 0.27, 0.025])                  
            cb  = fig.colorbar(h1, cax=cax, orientation="horizontal",panchor=(0.5,0.5),shrink=0.3,ticks=ticks)
            cb.set_label(r'Sea Surface Temperature [$^\circ\!$C]', fontsize=9, color='0.2',labelpad=0)
            cb.ax.tick_params(labelsize=9, length=2, color='0.2', labelcolor='0.2',direction='in') 
            cb.set_ticks(ticks)
        if contourf_var=='2':
            h1  = m.contourf(lon_var, lat_var, var, clevs,latlon=True,cmap=cmap)  
            cax = fig.add_axes([0.37, 0.01, 0.27, 0.025])   
            cb  = fig.colorbar(h1, cax=cax, orientation="horizontal",panchor=(0.5,0.5),shrink=0.3,ticks=ticks)
            cb.set_label(r'Sea Surface Salinity [PSU]', fontsize=9, color='0.2',labelpad=-0)
            cb.ax.tick_params(labelsize=5, length=2, color='0.2', labelcolor='0.2',direction='in') 
            cb.set_ticks(ticks)
        if contourf_var=='3':
            x1, y1 = m(to_np(lon_wrf), to_np(lat_wrf))
            h1  = ax.contourf(x1, y1, to_np(var), clevs,cmap=cmap,latlon=True,extend="both")  
            cax = fig.add_axes([0.37, 0.01, 0.27, 0.025])     
            cb  = fig.colorbar(h1, cax=cax, orientation="horizontal",panchor=(0.5,0.5),shrink=0.3,ticks=ticks)
            cb.set_label(r'Accumulated Precipitation [mm]', fontsize=9, color='0.2',labelpad=0)
            cb.ax.tick_params(labelsize=9, length=2, color='0.2', labelcolor='0.2',direction='in') 
            cb.set_ticks(ticks)
        if contourf_var=='4':
            x1, y1 = m(to_np(lon_wrf), to_np(lat_wrf))
            h1  = ax.contourf(x1, y1, to_np(var), clevs,cmap=cmap,latlon=True,extend="both")  
            cax = fig.add_axes([0.37, 0.01, 0.27, 0.025])     
            cb  = fig.colorbar(h1, cax=cax, orientation="horizontal",panchor=(0.5,0.5),shrink=0.3,ticks=ticks)
            cb.set_label(r'Accumulated Precipitation [W m.s⁻²]', fontsize=9, color='0.2',labelpad=0)
            cb.ax.tick_params(labelsize=9, length=2, color='0.2', labelcolor='0.2',direction='in') 
            cb.set_ticks(ticks)            
        if contourf_var=='5':
            h1  = m.contourf(lon_var, lat_var, var, clevs,latlon=True,cmap=cmap,extend="both")  
            cax = fig.add_axes([0.37, 0.01, 0.27, 0.025])   
            cb  = fig.colorbar(h1, cax=cax, orientation="horizontal",panchor=(0.5,0.5),shrink=0.3,ticks=ticks)
            cb.set_label(r'Fractional sea ice cover [%]', fontsize=9, color='0.2',labelpad=-0)
            cb.ax.tick_params(labelsize=9, length=2, color='0.2', labelcolor='0.2',direction='in') 
            cb.set_ticks(ticks)  

    # 3.8. Save figures.
    if contourf_var=='1':
        try:
            os.makedirs("sst")
        except FileExistsError:
            pass
        if ppt_fig==True:
            plt.savefig('./sst/temp_{0:03d}.png'.format(i), transparent=True, bbox_inches = 'tight', pad_inches=0, dpi=250)
        if ppt_fig==False:
            plt.savefig('./sst/temp_{0:03d}.png'.format(i), transparent=False, bbox_inches = 'tight', pad_inches=0, dpi=250)        
        plt.clf()
    elif contourf_var=='2':
        try:
            os.makedirs("sss")
        except FileExistsError:
            pass
        if ppt_fig==True:
            plt.savefig('./sss/salt_{0:03d}.png'.format(i), transparent=True, bbox_inches = 'tight', pad_inches=0, dpi=250)
        if ppt_fig==False:
            plt.savefig('./sss/salt_{0:03d}.png'.format(i), transparent=False, bbox_inches = 'tight', pad_inches=0, dpi=250)            
        plt.clf()   
    elif contourf_var=='3':
        try:
            os.makedirs("precipitation")
        except FileExistsError:
            pass
        if ppt_fig==True:
            plt.savefig('./precipitation/prec_{0:03d}.png'.format(i), transparent=True, bbox_inches = 'tight', pad_inches=0, dpi=250)
        if ppt_fig==False:
            plt.savefig('./precipitation/prec_{0:03d}.png'.format(i), transparent=False, bbox_inches = 'tight', pad_inches=0, dpi=250)           
        plt.clf()
    elif contourf_var=='4':
        try:
            os.makedirs("heat_fluxes")
        except FileExistsError:
            pass
        if ppt_fig==True:
            plt.savefig('./heat_fluxes/fluxes_{0:03d}.png'.format(i), transparent=True, bbox_inches = 'tight', pad_inches=0, dpi=250)
        if ppt_fig==False:
            plt.savefig('./heat_fluxes/fluxes_{0:03d}.png'.format(i), transparent=False, bbox_inches = 'tight', pad_inches=0, dpi=250)            
        plt.clf()
    elif contourf_var=='5':
        try:
            os.makedirs("ice")
        except FileExistsError:
            pass
        if ppt_fig==True:
            plt.savefig('./ice/ice_{0:03d}.png'.format(i), transparent=True, bbox_inches = 'tight', pad_inches=0, dpi=250)
        if ppt_fig==False:
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
