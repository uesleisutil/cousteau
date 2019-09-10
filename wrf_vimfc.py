"""
File name:      wrf_vimfc.py
Author:         Ueslei Adriano Sutil
Email:          uesleisutil1@gmail.com
Created:        26 August 2019
<<<<<<< HEAD
Last modified:  02 September 2019
Version:        1.2
=======
Last modified:  26 August 2019
Version:        1.0
>>>>>>> 26292d69271599df1c602f8e8c2e1d9669c544c6
Python:         3.7.1

Plot Vertical Integrated Moisture Flux Convergence from WRF output.

First, I created a netCDF file using the NCL script wrf_vimfc.ncl, then 
I used the original WRF output and the netCDF file to create the plots.

References:
Banacos, P. C. and Schultz, D. M. The Use of Moisture Flux Convergence in Forecasting 
Convective Initiation: Historical and Poeration Perspectives. Weather and Forecasting. 
v. 20, p. 351-366, 2005.

van Zomeren, J.; van Delden, A. Vertically integrated moisture flux convergence
as a predictor of thunderstorms. Atmospheric Research, v. 83, p. 435-445, 2007.
<<<<<<< HEAD
=======

>>>>>>> 26292d69271599df1c602f8e8c2e1d9669c544c6
"""

# Python preable.
import numpy                as     np
from   wrf                  import getvar,interplevel,latlon_coords,to_np,extract_times
import matplotlib.pyplot    as     plt
<<<<<<< HEAD
import matplotlib
import matplotlib.cm
=======
>>>>>>> 26292d69271599df1c602f8e8c2e1d9669c544c6
from   mpl_toolkits.basemap import Basemap
from   progress.bar         import IncrementalBar
import pandas               as     pd
import netCDF4
import os
import cmocean
<<<<<<< HEAD
import sys
import warnings
if not sys.warnoptions:
    warnings.simplefilter("ignore")
matplotlib.use('Agg')
=======
>>>>>>> 26292d69271599df1c602f8e8c2e1d9669c544c6

# Set file directory and lon/lat.
wrf_file      = '/media/ueslei/Ueslei/INPE/PCI/Projetos/SC_2008/Outputs/normal/vimfc.nc'
wrf_file2     = '/media/ueslei/Ueslei/INPE/PCI/Projetos/SC_2008/Outputs/normal/wrf.nc'
nc_wrf        = netCDF4.Dataset(wrf_file)
nc_wrf2       = netCDF4.Dataset(wrf_file2)
bbox          = [-52.5, -45.5, -30.5, -25.5]
<<<<<<< HEAD
clevs         = np.arange(-0.10,0.105,0.005)
ticks         = np.arange(min(clevs),max(clevs),0.05)  
cmap          = matplotlib.cm.bwr
=======
clevs         = np.arange(-0.1,0.105,0.005)
ticks         = np.arange(min(clevs),max(clevs),0.05)  
cmap          = cmocean.cm.balance   
>>>>>>> 26292d69271599df1c602f8e8c2e1d9669c544c6
create_video  = True
ppt_fig       = False # If True, generate a figure with transparent background.

# Open variables.
foo        = getvar(nc_wrf2, "td2",units="degC")
foo2       = nc_wrf2.variables['LH'][:,0,0]
ntimes     = len(foo2)
range_loop = [i for i in range(0,ntimes,1)]
bar        = IncrementalBar(max=len(range_loop))

for i in range_loop:
    timestr1  = extract_times(nc_wrf2,timeidx=None,meta=False,do_xtime=False)
    timestr11 = timestr1[i]
    timestr   = pd.to_datetime(timestr11, format="%b %d %Y %H:%M")

    # Load variables.
    uvmet   = getvar(nc_wrf2, "uvmet",timeidx=i,units="m s-1")
    ua      = uvmet[0,:,:,:]
    va      = uvmet[1,:,:,:]   
    uvmet10 = getvar(nc_wrf2, "uvmet10",timeidx=i,units="m s-1")
    p       = getvar(nc_wrf2,"pres",timeidx=i,units="hPa")
    u850    = interplevel(ua, p, 850.0)
    v850    = interplevel(va, p, 850.0)
    vimfc   = nc_wrf.variables['vimfc'][i,:,:]
    lats_wrf, lons_wrf = latlon_coords(foo, as_np=True)

    # Basemap options.
    m   = Basemap(projection='merc',llcrnrlat=bbox[2],urcrnrlat=bbox[3],llcrnrlon=bbox[0],urcrnrlon=bbox[1],lat_ts=30,resolution='i')
    fig = plt.figure(1,figsize=(10,8))
    ax  = fig.add_subplot(111)
    plt.xlabel('Longitude'u' [\N{DEGREE SIGN}]',labelpad=18,size=10)
    plt.ylabel('Latitude'u' [\N{DEGREE SIGN}]',labelpad=33,size=10)
    ax  = fig.add_subplot(111)
    
    # Draw parallels and meridians.
    m.drawparallels(np.arange(-90.,120.,1.), linewidth=0.00, color='black', labels=[1,0,0,1],labelstyle="N/S",fontsize=10)
    m.drawmeridians(np.arange(-180.,180.,1.), linewidth=0.00,color='black', labels=[1,0,0,1],labelstyle="N/S",fontsize=10)
    m.drawcountries(color = '#000000')
    m.drawcoastlines(color = '#000000')
    m.drawstates(color = '#000000')   

    # Plot winds as windbarbs.
    spacing = 45
    scale   = 0.03
    x, y    = m(to_np(lons_wrf), to_np(lats_wrf))
    W       = ax.quiver(x[::spacing,::spacing], y[::spacing,::spacing], uvmet10[0,::spacing, ::spacing],uvmet10[1,::spacing, ::spacing], pivot='middle',scale=8/scale, zorder=1e35, width=0.005,color='black',headlength=3,headaxislength=2.8 )
    wk      = ax.quiverkey(W, 0.76, -0.13, 10, ' Wind Vector at 10 m\n 10 m.s⁻¹ ', coordinates='axes',color='#444444',labelsep=0.05, labelcolor='black',fontproperties={'size': '9'})
    K       = ax.quiver(x[::spacing,::spacing], y[::spacing,::spacing],u850[::spacing, ::spacing],v850[::spacing, ::spacing],pivot='middle',scale=8/scale, zorder=1e35, width=0.005,color='blue',headlength=3,headaxislength=2.8 )
    Skk     = ax.quiverkey(K, 0.20, -0.13, 10, ' Wind Vector at 850 hPa\n 10 m.s⁻¹ ',coordinates='axes',color='blue',labelsep=0.05, labelcolor='black',fontproperties={'size': '9'})

    # Plot contourf variable.
<<<<<<< HEAD
    h1  = m.contourf(lons_wrf, lats_wrf, vimfc, clevs,latlon=True,cmap=cmap,extend="both",spacing="proportional",midpoint=0)       
=======
    h1  = m.contourf(lons_wrf, lats_wrf, vimfc, clevs,latlon=True,cmap=cmap,extend="both")       
>>>>>>> 26292d69271599df1c602f8e8c2e1d9669c544c6
    cax = fig.add_axes([0.37, 0.01, 0.27, 0.025])     
    cb  = fig.colorbar(h1, cax=cax, orientation="horizontal",panchor=(0.5,0.5),shrink=0.3,ticks=ticks)
    cb.set_label(r'Vertically Integrated Moisture Flux Convergence [10⁻⁵ Kg m⁻²s⁻¹]', fontsize=9, color='0.2',labelpad=-0.5)
    cb.ax.tick_params(labelsize=8, length=2, color='0.2', labelcolor='0.2',direction='in') 
    cb.set_ticks(ticks)
    plt.title(timestr, fontsize=10, pad=500)
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
    os.system("ffmpeg -r 10 -pattern_type glob -i '*.png' -c:v libx264 -pix_fmt yuv420p -crf 15 ./wrf_flux.mp4")  
if create_video==False:  
    pass
