import numpy as np
from   numpy import ma
import matplotlib.pyplot as plt
import matplotlib.patheffects as PathEffects
import matplotlib.cm as cm
from   mpl_toolkits.basemap import Basemap
import netCDF4
from   roms_libs import *
import cmocean
from   OceanLab.utils import download_bathy
import subprocess, os
from   wrf import to_np, getvar, latlon_coords, smooth2d
from scipy.ndimage.filters import gaussian_filter

wrf_file      = '/media/ueslei/Ueslei/INPE/PCI/Antartic/Outputs/wrf.nc'
initloop      = 0
zlev          = -1 # Last sigma layer corresponds to surface.
nc_wrf        = netCDF4.Dataset(wrf_file)

for i in range(0,1):
    slp         = getvar(nc_wrf, "slp", i, units="hPa")
    slp         = smooth2d(slp,30)
    lats_wrf, lons_wrf = latlon_coords(slp)
    clevs_slp   = np.arange(980,1030,1)
    var   = nc_wrf.variables['SST'][i,:,:]
    lat_var = nc_wrf.variables['XLAT'][i,:,0]
    lon_var = nc_wrf.variables['XLONG'][i,0,:]
    clevs    = np.arange(0,501,5)
    ticks    = np.array([0,100,200,300,400,500])  
    cmap     = matplotlib.pyplot.jet()  #cmocean.cm.rain    
    m = Basemap(projection='spstere',boundinglat=-45,lon_0=90,resolution='l')   
    fig  = plt.figure(1,figsize=(10,8))
    plt.xlabel('Longitude'u' [\N{DEGREE SIGN}]',labelpad=11,size=6)
    plt.ylabel('Latitude'u' [\N{DEGREE SIGN}]',labelpad=17,size=6)
    ax   = fig.add_subplot(111)
    m.drawcoastlines(color = '#000000',linewidth=0.5)
    x_slp, y_slp = m(to_np(lons_wrf), to_np(lats_wrf))
    Ct_slp       = ax.contour(x_slp,y_slp,to_np(slp),clevs_slp,colors='white',latlon=True,linewidths=0.9,linestyles='solid')
    clbls_slp    = plt.clabel(Ct_slp,fmt='%i',inline=1, fontsize=5)
    plt.savefig('./antartic.png', transparent=False, bbox_inches = 'tight', pad_inches=0, dpi=250)
    h1  = m.contourf(to_np(lon_var), to_np(lat_var), to_np(var), clevs,latlon=True,cmap=cmap)  