### wrf_roms_cross_section.py ###################
# Author:        Ueslei Adriano Sutil           #
# Created:       25 aug 2017                    #
# Last modified: 29 aug 2017                    #
#                                               #
# About:         Create a cross-section from    #
#                WRF and ROMS outputs.          #
#################################################

import numpy               as np
import matplotlib.pyplot   as plt1
from   matplotlib.cm       import get_cmap
import cartopy.crs         as crs
from   cartopy.feature     import NaturalEarthFeature
from   cmocean             import cm
from   wrf                 import to_np, getvar, CoordPair, vertcross
import matplotlib.gridspec as gridspec
from   netCDF4             import Dataset,num2date
from   pyroms              import vgrid
# from   octant              import tools
import numpy               as np
import matplotlib.pyplot   as pltexitr

####### WRF options #######
# WRF input
filename_wrf = "/media/ueslei/Ueslei/INPE/2014/Outputs/WRF/wrf_I_t01.nc"
filename_wr = "/media/ueslei/Ueslei/INPE/2014/Outputs/WR/wr_I_t01.nc"
filename_wrs = "/media/ueslei/Ueslei/INPE/2014/Outputs/WRS/wrs_I_t01.nc"

# Create the start point and end point for the cross section
start_point = CoordPair(lat=-41, lon=-60)
end_point   = CoordPair(lat=-41, lon=-30)
# WRF max height (meters)
levels      = np.arange(0,601.,1.)
# WRF grid height spacing
wrf_spacing = 100
# WRF timestep
timeidx     = 80

####### ROMS options #######
# ROMS input
path        = '/media/ueslei/Ueslei/INPE/2014/Outputs/WRS/'
name        = 'cbm12_ocean_his.nc'


# ROMS time
ntime       = 40
# ROMS grid index
#ind_lat     = 125:420
ind_lon     = 250
# ROMS max depth
dep_max     = -600

###### Plot options ######
cmap1       = cm.balance
col_label   = ' '
clevs       = np.arange(-3,3,0.005) # (min,max,spacing)
cxlatlon    = 'Longitude'
fig_dir     = './'
fig_name    = 'temp2.png'
dpi_opt     = 300

################### Don't change above ###################

# WRF module
# Open the NetCDF file
ncfile_wrs = Dataset(filename_wrs)
ncfile_wr  = Dataset(filename_wr)
ncfile_wrf = Dataset(filename_wrf)

# Extract the model height and wind speed
z_wrs    = getvar(ncfile_wrs, "z")
temp_wrs = getvar(ncfile_wrs, "tc", timeidx=timeidx)


z_wr    = getvar(ncfile_wr, "z")
temp_wr = getvar(ncfile_wr, "tc", timeidx=timeidx)


z_wrf    = getvar(ncfile_wrf, "z")
temp_wrf = getvar(ncfile_wrf, "tc", timeidx=timeidx)


# Compute the vertical cross-section interpolation.  Also, include the lat/lon
# points along the cross-section.
temp_cross_wrs = vertcross(temp_wrs, z_wrs, wrfin=ncfile_wrs, start_point=start_point, end_point=end_point,
                       latlon=True, meta=True,levels=levels)

temp_cross_wr = vertcross(temp_wr, z_wr, wrfin=ncfile_wr, start_point=start_point, end_point=end_point,
                       latlon=True, meta=True,levels=levels)

temp_cross_wrf = vertcross(temp_wrf, z_wrf, wrfin=ncfile_wrf, start_point=start_point, end_point=end_point,
                       latlon=True, meta=True,levels=levels)

# Differences
#temp_cross_wrf_wrs = temp_cross_wrf-temp_cross_wrs
#temp_cross_wrf_wr  = temp_cross_wrf-temp_cross_wr
temp_cross_wrs_wr  = temp_cross_wrf-temp_cross_wrs


# Create the figure
fig = plt1.figure(1,figsize=(12,6))
ax  = fig.add_subplot(2,1,1)
ax  = plt1.axes()
gs = gridspec.GridSpec(2, 1, hspace=0)
ax.set_subplotspec(gs[0:1])
fig.tight_layout()

# Make the contour plot
temp_contours = ax.contourf(to_np(temp_cross_wrs_wr), cmap=cmap1,levels=clevs)

# Set the x-ticks to use latitude and longitude labels.
coord_pairs = to_np(temp_cross_wrs_wr.coords["xy_loc"])
x_ticks     = np.arange(coord_pairs.shape[0])
x_labels    = [pair.latlon_str(fmt="{:.2f}, {:.2f}") for pair in to_np(coord_pairs)]

# Disable xticks
plt1.setp( ax.get_xticklabels(), visible=False)
ax.xaxis.grid(False, which='minor')

# Set the y-ticks to be height.
vert_vals = to_np(temp_cross_wrs_wr.coords["vertical"])
v_ticks   = np.arange(vert_vals.shape[0])

ax.set_yticks(v_ticks[::wrf_spacing])
ax.set_yticklabels(vert_vals[::wrf_spacing])

# Set the x-axis and  y-axis labels
ax.set_ylabel("Altitude (m)", fontsize=12)

# ROMS module
# Set variables
data   = Dataset(path+name)
cs_r   = data.variables[u'Cs_r'][:]
hc     = data.variables[u'hc'][:]
s_rho  = data.variables[u's_rho'][:]
N      = len(s_rho)
rho_pt = data.variables[u'lon_rho'][260,120:485]
Vtrans = data.variables[u'Vtransform'][:]
zeta   = data.variables[u'zeta'][ntime,:,:]
depth  = data.variables[u'h'][:,:]
lat    = rho_pt.repeat(N)
lats   = lat.reshape(rho_pt.shape[0],N)
var0   = data.variables[u'temp'][ntime,:,260,120:485]
var    = var0.transpose()
time   = data.variables['ocean_time']


# Calculate the depth (in meters)
[z]=vgrid.z_r(depth,hc,N,s_rho,cs_r,zeta,Vtrans)
zz=z[:,260,120:485]
zt=zz.transpose()

# Setting the initial figure
fig1 = plt.figure(1,figsize=(12,6))

ax1 = fig1.add_subplot(gs[1])
fig1.tight_layout()
ax1.set_ylabel('Depth (m)',fontsize=12)
ax1.set_xlabel(cxlatlon, fontsize=12)
ax1.set_ylim(dep_max,zt.max())

# Plotting values
vare    = ax1.contourf(lats,zt,var,clevs,shading='flat',cmap=cmap1)
cb      = plt.colorbar(vare,ax=[ax,ax1])

#cb      = plt.colorbar(vare,shrink=1.0,orientation='horizontal',fraction=0.046, pad=0.04)
cb.set_label(col_label,fontsize=12)
plt.savefig(fig_dir+fig_name,dpi=dpi_opt,bbox_inches='tight')
