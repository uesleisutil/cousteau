#!/home/ueslei/anaconda3/bin/python python
# -*- coding: utf-8 -*-

"""
File name:      roms_wrf_horizontal.py
Author:         Ueslei Adriano Sutil
Email:          uesleisutil1@gmail.com
Created:        02 August 2017
Last modified:  30 September 2019
Version:        2.5
Python:         3.7.1

Creates vertical plots from ROMS (his) and WRF-ARW outputs.

WARNING: Since ROMS does not plot continent data, Use coordinates 
         only inside the oceanm otherwise it will show a generical
		 error. 

bbox = [lon_min,lon_max,lat_min,lat_max]

"""

import matplotlib.pyplot   as plt1
import matplotlib.pyplot   as plt
import cmocean
from   wrf                 import to_np, getvar, CoordPair, vertcross
import matplotlib.gridspec as gridspec
from   netCDF4             import Dataset
from   pyroms              import vgrid
import numpy               as np
from   roms_libs           import bbox2ij
plt.switch_backend('Agg')
plt.rcParams.update({'font.size': 9})
# 1. WRF options.
wrf_file    = '/media/ueslei/Ueslei/INPE/PCI/Projetos/SC_2008/Outputs/normal/antigo/wrf.nc'
wrf_start   = CoordPair(lat=-5, lon=-20)
wrf_end     = CoordPair(lat=5, lon=-20)
wrf_levels  = np.arange(0,601.,1.) # Height (in meters).
wrf_h_spac  = 100 # Height spacing in figure (in meters).
wrf_time    = 1

# 2. ROMS options.
roms_file   = '/media/ueslei/Ueslei/INPE/PCI/Projetos/SC_2008/Outputs/normal/antigo/roms.nc'
roms_time   = 1
roms_depth  = -1500 # (in meters)

# 3. Plotting options.
col_label   = 'Temperature [°C]'
clevs       = np.arange(10,30,0.01) # Min, max and spacing (In °C).
cxlatlon    = 'Latitude [°]'
fig_dir     = './'
fig_name    = 'wrf_roms_cross_section.png'
dpi_opt     = 300
cmap1       = plt.jet()  
plt.title('Cross-section at 20 °W',fontsize=10)

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
ax.set_ylabel("WRF Altitude [m]", fontsize=9)

# ROMS module
# Set variables
data   = Dataset(roms_file)
cs_r   = data.variables[u's_rho'][:]
hc     = 25 # hmin
s_rho  = data.variables[u's_rho'][:]
N      = len(s_rho)
rho_pt = data.variables[u'lat_rho'][121:243,480]
Vtrans = 2 # Vtransform
zeta   = data.variables[u'zeta'][:,:]
depth  = data.variables[u'h'][:,:]
lat    = rho_pt.repeat(N)
lats   = lat.reshape(rho_pt.shape[0],N)
var0   = data.variables[u'temp'][0,:,121:243,480]
var    = var0.transpose()
time   = data.variables['ocean_time']

# Calculate the depth (in meters)
[z]=vgrid.z_r(depth,hc,N,s_rho,cs_r,zeta,Vtrans)
zz=z[:,121:243,480]
zt=zz.transpose()

# Setting the initial figure
fig1 = plt.figure(1,figsize=(12,6))
ax1  = fig1.add_subplot(gs[1])
fig1.tight_layout()
ax1.set_ylabel('ROMS Depth [m]',fontsize=9)
ax1.set_xlabel(cxlatlon, fontsize=9)
ax1.set_ylim(roms_depth,zt.max())

# Plotting values
vare    = ax1.contourf(lats,zt,var,clevs,shading='flat',cmap=cmap1,extend="both")
cb      = plt.colorbar(vare,ax=[ax,ax1],label='[C]')
#cb      = plt.colorbar(vare,shrink=1.0,orientation='horizontal',fraction=0.046, pad=0.04)
cb.set_label(col_label, fontsize=9, color='0.2',labelpad=5)
cb.ax.tick_params(labelsize=9, length=2, color='0.2', labelcolor='0.2',direction='in') 
plt.savefig(fig_dir+fig_name,dpi=dpi_opt,bbox_inches='tight')