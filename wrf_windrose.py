"""
File name:      wrf_windrose.py
Author:         Ueslei Adriano Sutil
Email:          uesleisutil1@gmail.com
Created:        24 October 2019
Last modified:  24 October 2019
Version:        1.0
Python:         3.7.1

Given one latitude/longitude point, extract all WRF-ARW wind speed and direction at 
10 meters height rotated to Earth coordinates as .txt file, then plot its wind rose
plot
"""
# Import libraries.
from   windrose      import WindroseAxes
import pandas        as pd
import numpy         as np
from   matplotlib    import pyplot as plt
import matplotlib.cm as cm
from   wrf           import getvar, ll_to_xy, extract_times
from   progress.bar  import IncrementalBar
import csv
import netCDF4

# Open file and choose desired lat/lon.
wrf_file = '/media/ueslei/Ueslei/INPE/PCI/Projetos/SC_2008/Outputs/warm_sst_100/wrf.nc'
nc_wrf   = netCDF4.Dataset(wrf_file)
lat      = -26.555
lon      = -49.355

# Set time, lenght of the WRF file and the variables to stock the wind speed and direction data.
timestr           = extract_times(nc_wrf,timeidx=None,meta=False,do_xtime=False)
timestr           = pd.to_datetime(timestr, format="%b %d %Y")
foo               = nc_wrf.variables['LH']
ntimes            = len(foo)
range_loop        = [i for i in range(168,ntimes,1)]
bar               = IncrementalBar('', max=len(range_loop))
wspd10_earth_list = np.zeros([len(range_loop)])
wdir10_earth_list = np.zeros([len(range_loop)])

# Start looping through time and stock data inside a list.
for i in range(0, len(range_loop),1):
    latlon       = ll_to_xy(nc_wrf,lat,lon,timeidx=i,squeeze=True,meta=False,as_int=True,stagger='m')
    uvmet10      = getvar(nc_wrf, "uvmet10_wspd_wdir",  timeidx=i,  units="m s-1", meta=False)
    wspd10_earth = uvmet10[0,latlon[0],latlon[1]]
    wdir10_earth = uvmet10[1,latlon[0],latlon[1]]

    wspd10_earth_list[i] = round(wspd10_earth,1)
    wdir10_earth_list[i] = round(wdir10_earth,1)

    bar.next()
bar.finish()     

# Create a variable to save data as table.
rows = zip(timestr, wspd10_earth_list,wdir10_earth_list)

# Save the variables in a CSV file
with open('./wrf_windrose.csv', "w") as f:
    writer = csv.writer(f)
    for row in rows:
        writer.writerow(row)

# Read the CSV file.
csv_file      = pd.read_csv('./wrf_windrose.csv', keep_default_na=True, delimiter=',', header=None,names=['TIME','WSPD','WDIR'])
wspd10_earth  = csv_file['WSPD']
wspd10_earth  = wspd10_earth.values
wdir10_earth  = csv_file['WDIR']
wdir10_earth  = wdir10_earth.values

# Create wind rose figure and save.
fig  = plt.figure(1,figsize=(10,8)) 
ax   = fig.add_subplot(111)
ax   = WindroseAxes.from_ax()
plot = ax.bar(wdir10_earth, wspd10_earth, normed=True, opening=0.8, edgecolor='white')
ax.set_legend()
plt.savefig('wrf_windrose1.png',transparent=False, bbox_inches = 'tight', pad_inches=0, dpi=250)
