"""
File name:      roms_bathymetry.py
Author:         Ueslei Adriano Sutil
Email:          uesleisutil1@gmail.com
Created:        03 September 2017
Last modified:  25 July 2019
Version:        2.0
Python:         3.7.1

Creates horizontal bathymetry plot from ROMS output file, in meters.

bbox = [lon_min,lon_max,lat_min,lat_max]
"""

# Importe as bibliotecas.
import numpy                  as np
import matplotlib.pyplot      as plt
import matplotlib.patheffects as PathEffects
import matplotlib.cm          as cm
from   mpl_toolkits.basemap   import Basemap
import netCDF4
from   roms_libs              import *
import cmocean

# Defina os limites do mapa.

roms_file   = '/media/ueslei/Ueslei/INPE/PCI/Projetos/Antartic/Outputs/antartic_grd.nc'
nc_roms     = netCDF4.Dataset(roms_file)
bbox        = [-70., 0, -80., -50.]
lon_rho     = nc_roms.variables['lon_rho'][:]
lat_rho     = nc_roms.variables['lat_rho'][:]
i0,i1,j0,j1 = bbox2ij(lon_rho,lat_rho,bbox)
lon_var     = lon_rho[j0:j1, i0:i1]
lat_var     = lat_rho[j0:j1, i0:i1]
var         = nc_roms.variables['hraw'][0,j0:j1, i0:i1]
mask        = 1 - nc_roms.variables['mask_rho'][j0:j1, i0:i1]
clevs       = np.arange(-5000,200,50)
ticks       = np.arange(min(clevs),max(clevs),1000)  
cmap        = cmocean.cm.deep_r


m           = Basemap(projection='merc',llcrnrlat=bbox[2],urcrnrlat=bbox[3],llcrnrlon=bbox[0],urcrnrlon=bbox[1], lat_ts=30,resolution='i')
m.fillcontinents(color='#BDA973') 
m.drawcoastlines(color='black',linewidth=0.5)

plt.xlabel('Longitude'u' [\N{DEGREE SIGN}]',labelpad=11,size=6)
plt.ylabel('Latitude'u' [\N{DEGREE SIGN}]',labelpad=19,size=6)
fig  = plt.figure(1,figsize=(10,8))
ax   = fig.add_subplot(111)
m.drawparallels(np.arange(-90.,120.,5), linewidth=0.00, color='black', labels=[1,0,0,1],labelstyle="N/S",fontsize=6)
m.drawmeridians(np.arange(-180.,180.,10), linewidth=0.00,color='black', labels=[1,0,0,1],labelstyle="N/S",fontsize=6)

h1  = m.contourf(lon_var, lat_var, var, clevs,latlon=True,cmap=cmap) 
cax = fig.add_axes([0.37, 0.02, 0.27, 0.025]) 
cb  = fig.colorbar(h1, cax=cax, orientation="horizontal",panchor=(0.5,0.5),shrink=0.3,ticks=ticks)
cb.set_label(r'Bathymetry [m]', fontsize=5, color='0.2',labelpad=-0.5)
cb.ax.tick_params(labelsize=5, length=2, color='0.2', labelcolor='0.2',direction='in') 
cb.set_ticks(ticks)
plt.savefig('./batimetria.png', transparent=True, bbox_inches = 'tight', pad_inches=0, dpi=250)