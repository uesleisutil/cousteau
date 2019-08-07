from OceanLab.utils import download_bathy
from OceanLab.utils import save_pickle
import netCDF4
import numpy as np # módulo de funções matemáticas
from mpl_toolkits.basemap import Basemap # módulo de mapas
import matplotlib.pyplot as plt #módulo de criação de imagens

longitude_min,longitude_max = -90.5, 90.0 #aqui aumentamos o tamanho
latitude_min,latitude_max   = -90.0, 90.0 #do mapa

bLON,bLAT,BAT = download_bathy(lnd=longitude_min,lnu=longitude_max,
                               ltd=latitude_min,ltu=latitude_max )

BDATA = dict(LON=bLON,LAT=bLAT,BAT=BAT)
save_pickle(BDATA,'/home/uesleisutil/Documentos/Libs/Batimetria')