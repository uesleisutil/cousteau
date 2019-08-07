#!/home/ueslei/anaconda3/bin/python python
# -*- coding: utf-8 -*-

"""
File name:      roms_wrf_horizontal.py
Author:         Ueslei Adriano Sutil
Email:          uesleisutil1@gmail.com
Created:        02 August 2017
Last modified:  25 April 2019
Version:        2.0
Python:         3.7.1

Creates vertical plots from ROMS (his) and WRF-ARW outputs.

bbox = [lon_min,lon_max,lat_min,lat_max]
"""

import matplotlib.pyplot   as plt1
import matplotlib.pyplot   as plt
import cmocean
from   wrf                 import to_np, getvar, CoordPair, vertcross
import matplotlib.gridspec as gridspec
from   netCDF4             import Dataset,num2date
from   pyroms              import vgrid
import numpy               as np

# 1. WRF options.
wrf_file    = "/media/ueslei/Ueslei/INPE/PCI/SC_2008/Outputs/normal/wrf.nc"
wrf_start   = CoordPair(lat=0.008, lon=-49.61)
wrf_end     = CoordPair(lat=0.008, lon=1.1)
wrf_levels  = np.arange(0,601.,1.) # Height (in meters).
wrf_h_spac  = 100 # Height spacing in figure (in meters).
wrf_time    = 1

# 2. ROMS options.
roms_file   = '/media/ueslei/Ueslei/INPE/PCI/SC_2008/Outputs/normal/roms.nc'
roms_time   = 1
roms_lat    = 183
roms_lon    = 150
	#TODO: Change ROMS inds to lat lon
roms_depth  = -600 # (in meters)

# 3. Plotting options.
col_label   = 'Temperature [°C]'
clevs       = np.arange(5,30,0.01) # Min, max and spacing (In °C).
cxlatlon    = 'Longitude [°]'
fig_dir     = './'
fig_name    = 'wrf_roms_cross_section.png'
dpi_opt     = 300
cmap1       = cmocean.cm.thermal # Collor pallete.

# No need to change the above.

# WRF module
# Open the NetCDF file and extract variables.
ncfile = Dataset(wrf_file)
z      = getvar(ncfile, "z")

if wrf_time>1:
	temp = getvar(ncfile, "tc", timeidx=wrf_time)
	ua   = getvar(ncfile, "ua", timeidx=wrf_time)
	va   = getvar(ncfile, "va", timeidx=wrf_time)
elif wrf_time<=1:
	temp = getvar(ncfile, "tc")
	ua   = getvar(ncfile, "ua")
	va   = getvar(ncfile, "va")

# Compute the vertical cross-section interpolation.
temp_cross = vertcross(temp, z, wrfin=ncfile, start_point=wrf_start, end_point=wrf_end,
					   latlon=True, meta=True,levels=wrf_levels)

# Create the figure
fig = plt1.figure(1,figsize=(12,6))
ax  = fig.add_subplot(2,1,1)
ax  = plt1.axes()
gs  = gridspec.GridSpec(2, 1, hspace=0)
ax.set_subplotspec(gs[0:1])
fig.tight_layout()
temp_contours = ax.contourf(to_np(temp_cross), cmap=cmap1,levels=clevs)

# Set the x-ticks to use latitude and longitude labels.
coord_pairs = to_np(temp_cross.coords["xy_loc"])
x_ticks     = np.arange(coord_pairs.shape[0])
x_labels    = [pair.latlon_str(fmt="{:.2f}, {:.2f}") for pair in to_np(coord_pairs)]

# Disable xticks
plt1.setp( ax.get_xticklabels(), visible=False)
ax.xaxis.grid(False, which='minor')

# Set the y-ticks to be height.
vert_vals = to_np(temp_cross.coords["vertical"])
v_ticks   = np.arange(vert_vals.shape[0])

ax.set_yticks(v_ticks[::wrf_h_spac])
ax.set_yticklabels(vert_vals[::wrf_h_spac])

# Set the x-axis and  y-axis labels
ax.set_ylabel("WRF Altitude [m]", fontsize=12)

# ROMS module
# Set variables
data   = Dataset(roms_file)
cs_r   = data.variables[u'Cs_r'][:]
hc     = data.variables[u'hc'][:]
s_rho  = data.variables[u's_rho'][:]
N      = len(s_rho)
rho_pt = data.variables[u'lon_rho'][260,100:300]
Vtrans = data.variables[u'Vtransform'][:]
zeta   = data.variables[u'zeta'][roms_time,:,:]
depth  = data.variables[u'h'][:,:]
lat    = rho_pt.repeat(N)
lats   = lat.reshape(rho_pt.shape[0],N)
var0   = data.variables[u'temp'][roms_time,:,260,100:300]
var    = var0.transpose()
time   = data.variables['ocean_time']


# Calculate the depth (in meters)
[z]=vgrid.z_r(depth,hc,N,s_rho,cs_r,zeta,Vtrans)
zz=z[:,260,100:300]
zt=zz.transpose()

# Setting the initial figure
fig1 = plt.figure(1,figsize=(12,6))
ax1 = fig1.add_subplot(gs[1])
fig1.tight_layout()
ax1.set_ylabel('ROMS Depth [m]',fontsize=12)
ax1.set_xlabel(cxlatlon, fontsize=12)
ax1.set_ylim(roms_depth,zt.max())

# Plotting values
vare    = ax1.contourf(lats,zt,var,clevs,shading='flat',cmap=cmap1)
cb      = plt.colorbar(vare,ax=[ax,ax1],label='[C]')
#cb      = plt.colorbar(vare,shrink=1.0,orientation='horizontal',fraction=0.046, pad=0.04)
cb.set_label(col_label,fontsize=12)
plt.savefig(fig_dir+fig_name,dpi=dpi_opt,bbox_inches='tight')
