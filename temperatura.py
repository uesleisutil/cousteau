from __future__ import unicode_literals
import numpy as np
import mpl_toolkits
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import os, sys, string
from datetime import datetime, timedelta
from cmocean import cm
import netCDF4
from OceanLab.utils import load_pickle
import matplotlib.patheffects as PathEffects

_author_   = 'Ueslei Adriano Sutil'
_email_    = 'ueslei@outlook.com'
_created_  = datetime(2017, 01, 31)
_modified_ = datetime(2017, 03, 14)
_version_  = "0.1.5"
_status_   = "Development"


# Abra o arquivo e carregue as variáveis.
nc = '/home/uesleisutil/Documentos/PCI/Anita/WRF/wrf.nc'
file = netCDF4.Dataset(nc)
lat  = file.variables['latitude'][:]
lon  = file.variables['longitude'][:]
data = file.variables['thetao'][0,0,:,:]

# Converta o dado de short para float.
data = np.array(data,'f')

# Plot o mapa usando o Basemap
#m = Basemap(llcrnrlat=(lat.min()+10),urcrnrlat=(lat.max()-10),llcrnrlon=(np.nanmin(lon)+15),urcrnrlon=(np.nanmax(lon)-20),rsphere=(6378137.00,6356752.3142),resolution='l',area_thresh=1000.,projection='lcc', lat_1=lat.max(),lon_0=lon.min())

m = Basemap(llcrnrlon=(lon.min()+10), llcrnrlat=(lat.min()+10), urcrnrlon=(lon.max()-10), urcrnrlat=(lat.max()),resolution='h',projection='merc')

#Adicione a linha de costa.
m.drawcoastlines(linewidth=1)

# Adicione os continentes.
m.fillcontinents(color='#BDA973')

# Adicione os paralelos e meridianos .
m.drawparallels(np.arange(-90.,120.,10.), labels=[1,0,0,1],labelstyle="+/-")
m.drawmeridians(np.arange(-180.,180.,10.),labels=[1,0,0,1],labelstyle="+/-")

# Adicione países
m.drawcountries(color = '#000000')

# Adicione uma escala ao mapa.
m.drawmapscale(-35,-45,lon.min(),lat.max(),1000, barstyle='fancy', units='km', labelstyle='simple',fillcolor1='w');

# Adicione alguns argumentos para o plot.
kw = dict(levels=range(5,30,1),cmap=cm.amp,latlon=True)
lon,lat = np.meshgrid(lon, lat)
Cf = m.contourf(lon,lat,data,**kw);

# Adicione a colorbar.
cbar = plt.colorbar(Cf);

# Adicione a batimetria.
BDATA = load_pickle('/home/uesleisutil/Documentos/python_scripts/batimetria')
bLON,bLAT,BAT = BDATA['LON'],BDATA['LAT'],BDATA['BAT']

# Adicione argumentos para a batimetria
Ct = m.contour(bLON,bLAT,BAT,[-1000,-200],colors='k',latlon=True)
clbls = plt.clabel(Ct,fmt='%i')

# Adicione um título.
plt.title('March 2004')

# Adicione um nome para o eixo da colorbar.
cbar.ax.set_ylabel('Sea Surface Temperature [C]') #precisamos acessar o objeto ax dentro de cbar

# Salve o arquivo
plt.savefig('/home/uesleisutil/Documentos/python_scripts/temperatura.png',dpi=500)
