#!/home/uesleisutil/anaconda3/bin/python python
# -*- coding: utf-8 -*-

"""
File name:      wrf_cape.py
Author:         Ueslei Adriano Sutil
Email:          uesleisutil1@gmail.com
Created:        27 February 2019
Last modified:  13 August 2019
Version:        3.5
Python:         3.7.1

Creates horizontal Convective Avaliable Potential Energy plot from and WRF-ARW output.

bbox = [lon_min,lon_max,lat_min,lat_max]
"""

import numpy                 as     np
import matplotlib.pyplot     as     plt
from   mpl_toolkits.basemap  import Basemap
import netCDF4
import os
from   wrf                   import to_np, getvar, latlon_coords, extract_times, cape_2d
from   scipy.ndimage.filters import gaussian_filter
import pandas                as     pd
from   progress.bar          import IncrementalBar
from   roms_libs             import *
import cmocean
plt.switch_backend('Agg')
import warnings
warnings.filterwarnings("ignore")
 
wrf_file      = '/media/ueslei/Ueslei/INPE/PCI/Projetos/SC_2008/Outputs/normal/wrf.nc'
lonbounds     = [-57,-40] 
latbounds     = [-22,-33]
bbox          = [-57,-40,-34,-22.2]
plot_var      = True 
plot_wind     = False
create_video  = True
ppt_fig       = False # If True, it generates a figure with transparent background.
nc_wrf        = netCDF4.Dataset(wrf_file)
atemp         = nc_wrf.variables['Times']
ntimes        = len(atemp)
foo           = getvar(nc_wrf, "LH", 0)
lat1, lon1    = latlon_coords(foo, as_np=True)
lat2          = lat1[:,0]
lon2          = lon1[0,:]
lon           = lon2[lonbounds[0]:lonbounds[1]]
lat           = lat2[latbounds[1]:latbounds[0]]


range_loop = [i for i in range(168,ntimes,1)]
bar        = IncrementalBar(' ', max=len(range_loop))
for i in range_loop:   
    timestr1  = extract_times(nc_wrf,timeidx=None,meta=False,do_xtime=False)
    timestr11 = timestr1[i]
    timestr   = pd.to_datetime(timestr11, format="%b %d %Y %H:%M")

    if plot_wind==True:
        uvmet10 = getvar(nc_wrf, "uvmet10", i, units="m s-1")
    if plot_var==True:
        pres    = getvar(nc_wrf, "pressure", i)
        tk      = getvar(nc_wrf, "tk", i)
        qv      = getvar(nc_wrf, "QVAPOR", i)        
        geopt   = getvar(nc_wrf, "geopt", i)
        terrain = getvar(nc_wrf, "ter", i, units="m")
        psfc    = getvar(nc_wrf, "PSFC", i)
        psfc    = psfc/100
        var     = cape_2d(pres,tk,qv,geopt,terrain,psfc, ter_follow=True, meta=False)
        cape    = var[0,:,:]
        cin     = var[1,:,:]
        var     = cape
        
        clevs   = np.arange(0,5005,10)
        ticks   = np.array([0,1000,2000,3000,4000,5000])  
        cmap    = cmocean.cm.thermal    
        
    # 3.2. Create a figure.
    m   = Basemap(projection='merc',llcrnrlat=latbounds[1],urcrnrlat=latbounds[0],llcrnrlon=lonbounds[0],urcrnrlon=lonbounds[1], lat_ts=30,resolution='i')
    fig = plt.figure(1,figsize=(10,8))
    plt.xlabel('Longitude'u' [\N{DEGREE SIGN}]',labelpad=18,size=10)
    plt.ylabel('Latitude'u' [\N{DEGREE SIGN}]',labelpad=33,size=10)
    ax  = fig.add_subplot(111)
    plt.title(timestr, fontsize=12)

    # 3.3. Add coastline, continents and lat/lon.
    m.drawparallels(np.arange(-90.,120.,1.), linewidth=0.00, color='black', labels=[1,0,0,1],labelstyle="N/S",fontsize=10)
    m.drawmeridians(np.arange(-180.,180.,2.), linewidth=0.00,color='black', labels=[1,0,0,1],labelstyle="N/S",fontsize=10)
    #m.fillcontinents(color = '#000000')
    m.drawcountries(color = 'k',linewidth=0.5)
    m.drawcoastlines(color = 'k',linewidth=0.5)

    # 3.5. Plot wind speed at 10 meters.
    if plot_wind == True:
        x, y    = m(to_np(lon1), to_np(lat1))
        spacing = 35
        scale   = 0.03
        W       = ax.quiver(x[::spacing,::spacing], y[::spacing,::spacing], to_np(uvmet10[0,::spacing, ::spacing]),to_np(uvmet10[1,::spacing, ::spacing]),pivot='middle',scale=8/scale, zorder=1e35, width=0.005,color='gray',headlength=3, headaxislength=2.8 )
        wk      = ax.quiverkey(W, 0.76, -0.14, 10, ' Wind Vector at 10 m\n 10 m.s⁻¹ ', coordinates='axes',color='#444444',labelsep=0.05, labelcolor='black',fontproperties={'size': '9'})

    if plot_var==True:
        h1  = m.contourf(lon1, lat1, var, clevs,latlon=True,cmap=cmap,extend="both")  
        cax = fig.add_axes([0.37, 0.05, 0.27, 0.025])                   
        cb  = fig.colorbar(h1, cax=cax, orientation="horizontal",panchor=(0.5,0.5),shrink=0.3,ticks=ticks)
        cb.set_label(r'Convective Available Potential Energy [J kg-1]', fontsize=9, color='0.2',labelpad=0)
        cb.ax.tick_params(labelsize=9, length=2, color='0.2', labelcolor='0.2',direction='in') 
        cb.set_ticks(ticks)

    if plot_var==True:
        try:
            os.makedirs("wrf_out")
        except FileExistsError:
            pass
        if ppt_fig==True:
            plt.savefig('./wrf_out/cape_{0:03d}.png'.format(i), transparent=True, bbox_inches = 'tight', pad_inches=0, dpi=250)
        if ppt_fig==False:
            plt.savefig('./wrf_out/cape_{0:03d}.png'.format(i), transparent=False, bbox_inches = 'tight', pad_inches=0, dpi=250)                 
        plt.clf()

    bar.next()
bar.finish()
if create_video==True:
    try:
        os.makedirs("wrf_out")
    except FileExistsError:
        pass
    cwd = os.getcwd()
    exists = os.path.isfile('./wrf_out/wrf_cape.mp4')
    if exists==True:
        os.system("rm -rf ./wrf_out/wrf_cape.mp4")
        os.system("ffmpeg -r 10 -pattern_type glob -i '"+cwd+"/wrf_out/*.png' -c:v libx264 -pix_fmt yuv420p -crf 15 ./wrf_out/wrf_cape.mp4")
    else:
        os.system("ffmpeg -r 10 -pattern_type glob -i '"+cwd+"/wrf_out/*.png' -c:v libx264 -pix_fmt yuv420p -crf 15 ./wrf_out/wrf_cape.mp4")
        print('Wish to delete the .png files? (1) Yes or (2) No.')
    removefiles = input()
    if removefiles=='1':
        os.system("rm -rf ./wrf_out/*.png")
    if removefiles=='2':
        pass