"""
File name:      etopo1_relief.py
Author:         Ueslei Adriano Sutil
Email:          uesleisutil1@gmail.com
Created:        03 September 2017
Last modified:  30 October 2020
Version:        3.0
Python:         3.8.3

Creates horizontal relief plot from ETOPO1 data.

bbox = [lon_min,lon_max,lat_min,lat_max]
"""


# Library import.
from   scipy.ndimage.filters import gaussian_filter
from   netCDF4               import Dataset
from   mpl_toolkits.basemap  import Basemap
import numpy                 as     np
import matplotlib.pyplot     as     plt
import cmocean
import matplotlib.colors     as     colors

def create_topo(lon_1=-49,lon_2=-33,lat_1=-34,lat_2=-20):
    """Plot ETOPO1 data as a bathymetry contour.
    """
    # Open data.
    etopodata = Dataset('/media/ueslei/Ueslei_HD/Scripts/Python/Auxiliares/OceanLab/ETOPO1.nc')

    # Read data.
    etopo_lon = etopodata.variables['x'][:]
    etopo_lat = etopodata.variables['y'][:]

    # Create subset.
    x_data = (etopo_lon>lon_1)&(etopo_lon<lon_2)
    y_data = (etopo_lat>lat_1)&(etopo_lat<lat_2)
    topo = etopodata.variables['z'][y_data,:][:,x_data]

    # Make grid.
    topo_lon,topo_lat = np.meshgrid(etopo_lon[x_data],etopo_lat[y_data])

    return topo_lon,topo_lat,topo

def create_topo2(lon_1=-49,lon_2=-33,lat_1=-34,lat_2=-20):
    """Plot ETOPO1 data as a bathymetry contour.
    """
    # Open data.
    etopodata = Dataset('/media/ueslei/Ueslei_HD/Scripts/Python/Auxiliares/OceanLab/ETOPO1_SC.nc')

    # Read data.
    etopo_lon = etopodata.variables['x'][:]
    etopo_lat = etopodata.variables['y'][:]

    # Create subset.
    x_data = (etopo_lon>lon_1)&(etopo_lon<lon_2)
    y_data = (etopo_lat>lat_1)&(etopo_lat<lat_2)
    topo2 = etopodata.variables['z'][y_data,:][:,x_data]

    # Make grid.
    topo_lon2,topo_lat2 = np.meshgrid(etopo_lon[x_data],etopo_lat[y_data])

    return topo_lon2,topo_lat2,topo2

def create_bathy(lon_1=-49,lon_2=-33,lat_1=-34,lat_2=-20):
    """Plot ETOPO1 data as a bathymetry contour.
    """
    # Open data.
    etopodata = Dataset('/media/ueslei/Ueslei_HD/Scripts/Python/Auxiliares/OceanLab/ETOPO1.nc')

    # Read data.
    etopo_lon = etopodata.variables['x'][:]
    etopo_lat = etopodata.variables['y'][:]

    # Create subset.
    x_data = (etopo_lon>lon_1)&(etopo_lon<lon_2)
    y_data = (etopo_lat>lat_1)&(etopo_lat<lat_2)
    bathy_topo = etopodata.variables['z'][y_data,:][:,x_data]

    # Make grid.
    bathy_lon,bathy_lat = np.meshgrid(etopo_lon[x_data],etopo_lat[y_data])

    return bathy_lon,bathy_lat,bathy_topo

bbox = [-53,-44,-31,-23] 
topo_lon,topo_lat,topo = create_topo(lon_1=bbox[0],lon_2=bbox[1],lat_1=bbox[2],lat_2=bbox[3])
	
m = Basemap(projection='merc',llcrnrlat=bbox[2],urcrnrlat=bbox[3],llcrnrlon=bbox[0],urcrnrlon=bbox[1], lat_ts=30,resolution='i')
fig  = plt.figure(1,figsize=(10,8))  
plt.xlabel('Longitude'u' [\N{DEGREE SIGN}]',labelpad=18,size=9)
plt.ylabel('Latitude'u' [\N{DEGREE SIGN}]',labelpad=33,size=9)
ax   = fig.add_subplot(111) 
m.drawparallels(np.arange(-90.,120.,1.), linewidth=0.00, color='black', labels=[1,0,0,1],labelstyle="N/S",fontsize=10)
m.drawmeridians(np.arange(-180.,180.,2.), linewidth=0.00,color='black', labels=[1,0,0,1],labelstyle="N/S",fontsize=10)

m.drawcountries(color = '#000000',linewidth=0.5)
m.drawcoastlines(color = '#000000',linewidth=0.5)
m.drawstates(color = '#000000',linewidth=0.5)  


bathy_levels     = [-1000,-200]        
bathy_lon,bathy_lat,bathy_topo = create_bathy(lon_1=bbox[0],lon_2=bbox[1],lat_1=bbox[2],lat_2=bbox[3])
Ct               = m.contour(gaussian_filter(bathy_lon,2),gaussian_filter(bathy_lat,2),gaussian_filter(bathy_topo,2),bathy_levels,colors='black',latlon=True,linewidths=0.7,linestyles='solid')
manual_locations = [(368129,16411),(409187,20034)]       
clbls            = plt.clabel(Ct,fmt='%i', fontsize=9,manual=manual_locations,colors="black")

topo_levels       = [300,500]        
topo_lon2,topo_lat2,topo2 = create_topo2(lon_1=bbox[0],lon_2=bbox[1],lat_1=bbox[2],lat_2=bbox[3])
Ct2               = m.contour(gaussian_filter(topo_lon2,2),gaussian_filter(topo_lat2,2),gaussian_filter(topo2,2),topo_levels,colors='white',latlon=True,linewidths=0.7,linestyles='solid')
manual_locations2 = [(358462,301659),(292944,405689)]       
clbls2            = plt.clabel(Ct2,fmt='%i', fontsize=9,manual=manual_locations2,colors="white")

ax.text(397297,487414, 'Blumenau',color='red',fontsize=9, bbox=dict(fill=True, edgecolor='black', linewidth=0,alpha=0.3))
plt.plot(387827,464829, color='red', marker='o', linestyle='dashed',linewidth=2, markersize=6)

ax.text(420550,419981, 'Major Gercino',color='red',fontsize=9, bbox=dict(fill=True, edgecolor='black', linewidth=0,alpha=0.3))
plt.plot(412055,397620, color='red', marker='o', linestyle='dashed',linewidth=2, markersize=6)


clevs  = np.arange(-3000,1820,5)
ticks  = np.arange(min(clevs),max(clevs),1000)  
cmap   = cmocean.cm.topo
x1, y1 = m(topo_lon, topo_lat)
h1     = ax.contourf(x1, y1, topo, clevs,cmap=cmap,latlon=True,vmin=-3000, vmax=3000,extend="both") 
cax    = fig.add_axes([0.37, 0.039, 0.27, 0.025]) 
cb     = fig.colorbar(h1, cax=cax, orientation="horizontal",panchor=(0.5,0.5),shrink=0.3,ticks=ticks)   
cb.set_label(r'Relief [m]', fontsize=9, color='0.2',labelpad=0) 
cb.ax.tick_params(labelsize=9, length=2, color='0.2', labelcolor='0.2',direction='in') 
cb.set_ticks(ticks) 
plt.show()
plt.savefig('./topo.png', transparent=False, bbox_inches = 'tight', pad_inches=0, dpi=100)    