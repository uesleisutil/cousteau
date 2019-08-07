from __future__ import unicode_literals
import numpy as np
import mpl_toolkits
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import os, sys, string
from datetime import datetime, timedelta
from cmocean import cm
import netCDF4
from  wrf import to_np, getvar, CoordPair

_author_   = 'Ueslei Adriano Sutil'
_email_    = 'ueslei@outlook.com'

# Abra o arquivo e carregue as variáveis.
nc   = '/media/ueslei/Ueslei/INPE/PCI/SC_2008/Outputs/normal/wrf_normal.nc'
file = netCDF4.Dataset(nc)

# Defina algumas opções
titulo = 'Marine Atmospheric Boundary Layer Stability \n Tropical Storm Anita'
ylabel = '[C]'
salvar = './mabl.png'

# Carregue as variáveis
lat  = file.variables['XLAT'][0,:,0]
lon  = file.variables['XLONG'][0,0,:]
sst  = file.variables['SST'][0,:,:]
sst  = sst-273.15
psfc = file.variables['PSFC'][0,:,:]
psfc = psfc/100
t2   = file.variables['T2'][0,:,:]
t2   = t2-273.15
stab = t2-sst

m = Basemap(llcrnrlon=(lon.min()+10), llcrnrlat=lat.min(), urcrnrlon=lon.max(), urcrnrlat=lat.max(),resolution='h',projection='merc')

#Adicione a linha de costa.
m.drawcoastlines(linewidth=1)

# Adicione os continentes.
m.fillcontinents(color='#000000')

# Adicione os paralelos e meridianos .
m.drawparallels(np.arange(-90.,120.,10.), linewidth=0.1, dashes = [1,5], color='black', labels=[1,0,0,1],labelstyle="+/-")
m.drawmeridians(np.arange(-180.,180.,10.), linewidth=0.1,dashes = [1,5], color='black', labels=[1,0,0,1],labelstyle="+/-")

# Adicione países
m.drawcountries(color = '#ffffff')

# Adicione uma escala ao mapa.
m.drawmapscale(-25,-50,lon.min(),lat.max(),1500, barstyle='fancy', units='km', labelstyle='simple',fillcolor1='w')

# Defina para plotar uma variável.
kw = dict(levels=range(-5,5,1),cmap=cm.balance,latlon=True)
lon,lat = np.meshgrid(lat, lat)
Cf = m.contourf(lon,lat,stab,**kw)
cbar = plt.colorbar(Cf)
cbar.ax.set_ylabel(ylabel) #precisamos acessar o objeto ax dentro de cbar

Cb = m.contour(stab,20,linewidths=2.0, colors='darkgreen', linestyles='-')
# Adicione um título.
plt.title(titulo)

# Salve o arquivo
plt.show()
plt.savefig(salvar,dpi=500)
