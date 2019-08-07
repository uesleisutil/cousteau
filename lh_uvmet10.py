'''
Make LH, SLP and UVMET10 plots with WRF-Python tool from WRF-ARW output using Cartopy, Matplotlib and Cmocean colormap.
'''

# Import the libraries.
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import cartopy.crs as crs
from cartopy.feature import NaturalEarthFeature
import cartopy.feature as cfeature
from wrf import to_np, getvar, smooth2d, get_cartopy, cartopy_xlim, cartopy_ylim, latlon_coords, ALL_TIMES
from cmocean import cm
import numpy as np
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import os, sys, string
from datetime import datetime, timedelta

# Author and script informations.
_author_   = 'Ueslei Adriano Sutil'
_email_    = 'ueslei@outlook.com'
_created_  = datetime(2017, 07, 20)
_modified_ = datetime(2017, 07, 20)
_version_  = "0.1"
_status_   = "Development"

# Open the NetCDF file.
ncfile = Dataset("/home/uesleisutil/Documentos/INPE/PCI/2014/Outputs/wrs_I_t02.nc")
ncfile2 = Dataset("/home/uesleisutil/Documentos/INPE/PCI/2014/Outputs/wrs_I_t02_LH.nc")

# Get the data.
slp     = getvar(ncfile, "slp", timeidx=78)
uvmet10 = getvar(ncfile, "uvmet10", units="km h-1", timeidx=78)
u10     = uvmet10[1,:,:]
v10     = uvmet10[0,:,:]
lh      = getvar(ncfile2, "LH", timeidx=78)

# Smooth the sea level pressure since it tends to be noisy near the mountains.
smooth_slp = smooth2d(slp, 10)

# Get the latitude and longitude points.
lats, lons = latlon_coords(slp)

# Get the cartopy mapping object.
cart_proj = get_cartopy(slp)

# Create a figure.
fig = plt.figure(figsize=(12,9))

# Set the GeoAxes to the projection used by WRF.
ax = plt.axes(projection=cart_proj)

# Download and add the states and coastlines.
states = NaturalEarthFeature(category='cultural', scale='50m', facecolor='none',name='admin_0_countries')
ax.coastlines('50m', linewidth=1.7)
ax.add_feature(cfeature.BORDERS,linewidth=1.7)

# Make the contour outlines and filled contours for the Latent Heat Flux.
clevs2 = np.arange(-100,650.,5.)
c2 = plt.contourf(to_np(lons), to_np(lats), to_np(lh), clevs2, transform=crs.PlateCarree(),cmap=cm.thermal,zorder=1)
plt.colorbar(ax=ax, shrink=.52,orientation='horizontal',aspect=20,fraction=0.026, pad=0.04)

# Plot the wind barbs for the U10 and V10 from UVMET10.
c3 = plt.barbs(to_np(lons[::75,::75]), to_np(lats[::75,::75]), to_np(u10[::75, ::75]), to_np(v10[::75, ::75]),color='black',transform=crs.PlateCarree(), length=6, zorder=3)

# Make the contour outlines and filled contours for the smoothed Sea Level Pressure.
clevs1 = np.arange(980,1040.,2.)
c1 = ax.contour(lons, lats, to_np(smooth_slp), levels=clevs1, colors="white",transform=crs.PlateCarree(), linewidths=1.3,zorder=2)
plt.clabel(c1,inline=1,inline_spacing=1,fontsize=8,fmt='%1.0f',colors='w')

# Set the map limits.  Not really necessary, but used for demonstration.
ax.set_xlim(cartopy_xlim(smooth_slp))
ax.set_ylim(cartopy_ylim(smooth_slp))

ax.set_xlabel('Latitude (degrees)')
ax.set_ylabel('Longitude (degrees)')

# Add the gridlines.
g1 = ax.gridlines(color="gray", linestyle="--",draw_labels=True,linewidth=0.2)
g1.xlabels_top = False
g1.ylabels_right= False

# Save.
plt.title("COA_WRS (ROMS+WRF+SWAN)")
plt.savefig('/home/uesleisutil/Documentos/Scripts/Python/slp_200.png',dpi=200,bbox_inches='tight')
