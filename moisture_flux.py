#!/home/uesleisutil/anaconda3/bin/python python
# -*- coding: utf-8 -*-

"""
File name:      wrf_vimfc.py
Author:         Ueslei Adriano Sutil
Email:          uesleisutil1@gmail.com
Created:        26 August 2019
Last modified:  26 August 2019
Version:        1.0
Python:         3.7.1

Calculate Vertical Integrated Moisture Flux Convergence from WRF output then
plot as contour.

Banacos, P. C. and Schultz, D. M. The Use of Moisture Flux Convergence in Forecasting 
Convective Initiation: Historical and Poeration Perspectives. Weather and Forecasting. 
v. 20, p. 351-366, 2005.

van Zomeren, J.; van Delden, A. Vertically integrated moisture flux convergence
as a predictor of thunderstorms. Atmospheric Research, v. 83, p. 435-445, 2007.


UNDER CONSTRUCTION!!!!!!!!
"""

# Python preable.
import numpy                as     np
from   wrf                  import getvar, interplevel, latlon_coords
import netCDF4
import matplotlib.pyplot    as     plt
from   mpl_toolkits.basemap import Basemap
import cmocean
from   progress.bar          import IncrementalBar
import os

# Set file directory and lon/lat.
wrf_file      = '/media/ueslei/Ueslei/INPE/PCI/Projetos/SC_2008/Outputs/normal/wrf.nc'
nc_wrf        = netCDF4.Dataset(wrf_file)
bbox          = [-52.5, -45.5, -30.5, -25.5]
clevs         = np.arange(-5,5,0.5)
ticks         = np.arange(min(clevs),max(clevs),2)  
cmap          = cmocean.cm.balance   
create_video  = False
ppt_fig       = False # If True, it will generate a figure with transparent background.

# Open variables.
foo        = nc_wrf.variables['LH'][:,0,0]
ntimes     = len(foo)
range_loop = [i for i in range(0,ntimes,1)]
bar        = IncrementalBar(max=len(range_loop))

for i in range_loop:
    ua         = getvar(nc_wrf, "ua",timeidx=i,units="m s-1")
    va         = getvar(nc_wrf, "va",timeidx=i,units="m s-1")
    p          = getvar(nc_wrf,"pres",timeidx=i,units="hPa")
    qvapor     = nc_wrf.variables['QVAPOR'][i,:,:,:]

    # Interpolate wind and specific humidity to 1000 and 700 hPa levels.
    # First layer;
    u1000 = interplevel(ua, p, 1000.0)
    lats_wrf, lons_wrf = latlon_coords(u1000, as_np=True)
    u1000 = u1000.values
    u925  = interplevel(ua, p, 925.0)
    u925  = u925.values
    v1000 = interplevel(va, p, 1000.0)
    v1000 = v1000.values
    v925  = interplevel(va, p, 925.0)
    v925  = v925.values
    q1000 = interplevel(qvapor, p, 1000.0)
    q1000 = q1000.values
    q925  = interplevel(qvapor, p, 925.0)
    q925  = q925.values
    h1000 = q1000/(1+q1000)
    h925  = q925/(1+q925)
    
    # Second layer.
    u850  = interplevel(ua, p, 850.0)
    u850  = u850.values
    v850  = interplevel(va, p, 850.0)
    v850  = v850.values
    q850  = interplevel(qvapor, p, 850.0)
    q850  = q850.values
    h850  = q850/(1+850)
    
    # Third layer.      
    u700  = interplevel(ua, p, 700.0)
    u700  = u700.values
    v700  = interplevel(va, p, 700.0)
    v700  = v700.values
    q700  = interplevel(qvapor, p, 700.0)
    q700  = q700.values
    h700  = q700/(1+q700)
    
    # Calculate VIMFC.
    # First layer.
    q1           = (q1000+q925)/2
    u1           = (u1000+u925)/2
    v1           = (v1000+v925)/2
    dp1          = 1000-925
    vimfc1       = ((q1*u1+q1*v1)*dp1)/9.81
    nan1         = np.isnan(vimfc1)
    vimfc1[nan1] = 0
    
    # Second layer.  
    q2           = (q925+q850)/2
    u2           = (u925+u850)/2
    v2           = (v925+v850)/2
    dp2          = 925-800
    vimfc2       = ((q2*u2+q2*v2)*dp2)/9.81
    nan2         = np.isnan(vimfc2)
    vimfc2[nan2] = 0
    
    # Third layer.   
    q3           = (q850+q700)/2
    u3           = (u850+u700)/2
    v3           = (v850+v700)/2
    dp3          = 850-700
    vimfc3       = ((q3*u3+q3*v3)*dp3)/9.81
    nan3         = np.isnan(vimfc3)
    vimfc3[nan3] = 0    
    
    # Sum layers.
    vimfc  = ((vimfc1 + vimfc2 + vimfc3))
    
    # Start plotting.
    m    = Basemap(projection='merc',llcrnrlat=bbox[2],urcrnrlat=bbox[3],llcrnrlon=bbox[0],urcrnrlon=bbox[1], lat_ts=30,resolution='i')
    fig  = plt.figure(1,figsize=(10,8))

    # Plor lat and lon labels.
    plt.xlabel('Longitude'u' [\N{DEGREE SIGN}]',labelpad=18,size=10)
    plt.ylabel('Latitude'u' [\N{DEGREE SIGN}]',labelpad=33,size=10)
    ax   = fig.add_subplot(111)

    # Draw parallels and meridians.
    m.drawparallels(np.arange(-90.,120.,1.), linewidth=0.00, color='black', labels=[1,0,0,1],labelstyle="N/S",fontsize=10)
    m.drawmeridians(np.arange(-180.,180.,1.), linewidth=0.00,color='black', labels=[1,0,0,1],labelstyle="N/S",fontsize=10)
    m.drawcountries(color = '#000000')
    m.drawcoastlines(color = '#000000')
    m.drawstates(color = '#000000')   

    # Plot variable.
    h1  = m.contourf(lons_wrf, lats_wrf, vimfc, clevs,latlon=True,cmap=cmap,extend="both")       

    # Set colorbar options.
    cax = fig.add_axes([0.37, 0.01, 0.27, 0.025])     
    cb  = fig.colorbar(h1, cax=cax, orientation="horizontal",panchor=(0.5,0.5),shrink=0.3,ticks=ticks)
    cb.set_label(r'Vertically Integrated Moisture Flux Convergence [10-5 Kg ms-2]', fontsize=11, color='0.2',labelpad=-0.5)
    cb.ax.tick_params(labelsize=10, length=2, color='0.2', labelcolor='0.2',direction='in') 
    cb.set_ticks(ticks)

    # Save figure
    if ppt_fig==True:
        plt.savefig('./vimfc_{0:03d}.png'.format(i), transparent=True, bbox_inches = 'tight', pad_inches=0, dpi=250)
    if ppt_fig==False:
        plt.savefig('./vimfc_{0:03d}.png'.format(i), transparent=False, bbox_inches = 'tight', pad_inches=0, dpi=250)
    plt.clf() 
    bar.next()
bar.finish()

# Create video, if True.
if create_video==True:  
    os.system("ffmpeg -r 10 -pattern_type glob -i '/*.png' -c:v libx264 -pix_fmt yuv420p -crf 15 ./wrf_flux.mp4")  
if create_video==False:  
    pass
