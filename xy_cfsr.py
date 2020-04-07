#!/home/uesleisutil/anaconda3/bin/python python
# -*- coding: utf-8 -*-

"""
File name:      xy_cfsr.py
Author:         Ueslei Adriano Sutil
Email:          uesleisutil1@gmail.com
Created:        07 April 2020
Last modified:  07 April 2020
Version:        1.0
Python:         3.7.4

Create XY graph from CFSR reanalysis.
"""

# Library import.
print("Phase 1: Starting program.")
import numpy                as np
import matplotlib.pyplot    as plt
import matplotlib.dates
from   netCDF4              import num2date
import datetime
import matplotlib
import netCDF4
import matplotlib.ticker as ticker

matplotlib.use('Agg')

# Import file and set bounding box.
print("Phase 2: Importing files.")
cfsr_file    = '/media/ueslei/Ueslei/INPE/Embarques/OP/OP38/2_Fase/Dados/CFSR/cfsr_oct2019.nc'
lonbounds    = [306,308] 
latbounds    = [-45,-43]

# Adjust variables inside the bouding box. 
print("Phase 3: Importing variables and cropping to boubing box.")
nc_cfsr      = netCDF4.Dataset(cfsr_file)
lon_cfsr     = nc_cfsr.variables['lon'][:]
lat_cfsr     = nc_cfsr.variables['lat'][:]
latli        = np.argmin(np.abs(lat_cfsr-latbounds[1]))
latui        = np.argmin(np.abs(lat_cfsr-latbounds[0])) 
lonli        = np.argmin(np.abs(lon_cfsr-lonbounds[0]))
lonui        = np.argmin(np.abs(lon_cfsr-lonbounds[1]))
lon_cfsr     = nc_cfsr.variables['lon'][lonli:lonui]
lat_cfsr     = nc_cfsr.variables['lat'][latli:latui]
lh_cfsr      = nc_cfsr.variables['LHTFL_L1'][:,latli:latui,lonli:lonui]             
sh_cfsr      = nc_cfsr.variables['SHTFL_L1'][:,latli:latui,lonli:lonui] 

# Load time and convert timestamp.
print("Phase 4: Load time.")
time         = nc_cfsr.variables['time'][:]
time_unit    = nc_cfsr.variables['time'].units 
date_time    = nc_cfsr.variables['time'].calendar 
dt           = netCDF4.num2date(time,time_unit,date_time)
dt           = [date_obj.strftime('%d/%m %Hh') for date_obj in dt]

# Calculate desired variables and average the bounding box area.
print("Phase 5: Calculate desired variables.")
th_cfsr      = lh_cfsr+sh_cfsr
th_cfsr_1    = np.average(th_cfsr,axis=1) # Average over lat.
th_cfsr_xy   = np.average(th_cfsr_1,axis=1) # Average over lon.
lh_cfsr_1    = np.average(lh_cfsr,axis=1) 
lh_cfsr_xy   = np.average(lh_cfsr_1,axis=1) 
sh_cfsr_1    = np.average(sh_cfsr,axis=1) 
sh_cfsr_xy   = np.average(sh_cfsr_1,axis=1) 

# Create figures subplots and adjustments.
print("Phase 6: Create figures subplots and adjustments.")
fig, axs      = plt.subplots(1, 1)
fig.subplots_adjust(wspace=0, hspace=0.12)
y_major_ticks = np.arange(0, 600.1, 50)

# Plot figure.
print("Phase 7: Plot figure.")
axs.plot(dt,th_cfsr_xy,label='Total Heat Flux',linestyle='-',marker='',linewidth=2,color='black')    
axs.plot(dt,lh_cfsr_xy,label='Latent Heat Flux',linestyle='-',marker='',linewidth=2,color='crimson')
axs.plot(dt,sh_cfsr_xy,label='Sensible Heat Flux',linestyle='-',marker='',linewidth=2,color='C0')
axs.legend(loc='upper right',fontsize=10)
axs.xaxis.set_major_formatter(ticker.FixedFormatter(dt))
spacing = 3
visible = axs.xaxis.get_ticklabels()[::spacing]
for label in axs.xaxis.get_ticklabels():
    if label not in visible:
        label.set_visible(False)
axs.tick_params(which='both',  width=2)
axs.tick_params(which='major', length=7)
axs.tick_params(which='minor', length=4)
axs.set_yticks(y_major_ticks)
axs.set_ylabel('Heat Fluxes \n [W m-2]',fontsize=12, labelpad=10)
axs.set_title(' ')

# Save figure.
print("Phase 8: Save figure.")
fig.set_size_inches(13, 8)
fig.autofmt_xdate()

fig.savefig('xy_cfsr_era5.png',transparent=False, bbox_inches = 'tight', pad_inches=0, dpi=250)
