#!/home/uesleisutil/anaconda3/bin/python python
# -*- coding: utf-8 -*-

"""
File name:      extract_data.py
Author:         Ueslei Adriano Sutil
Email:          uesleisutil1@gmail.com
Created:        15 August 2019
Last modified:  17 August 2019
Version:        1.0
Python:         3.7.1

Extract data from WRF output and save in ASCII format:
    - Air temperature at 2 meters (°C);
    - Wind velocity at 10 meters (m.s⁻¹);
    - Sea level pressure (hPa);
    - Relative humidity at 2 meters (%).
"""

# Import variables.
from   wrf import getvar,ll_to_xy,extract_times
import netCDF4
import csv
from   progress.bar import IncrementalBar
import numpy as np
from   datetime              import datetime
import pandas as pd

# Select WRF latitude/longitude and time-step.
lat      = -25.53
lon      = -48.51
inittime = 0

# Insert WRF file directory.
wrf_dir    = '/media/ueslei/Ueslei/INPE/PCI/Projetos/SC_2008/Outputs/normal/wrf.nc'
nc_file    = netCDF4.Dataset(wrf_dir)

# Create a variable in order to use the IncrementalBar.
timestr    = extract_times(nc_file,timeidx=None,meta=False,do_xtime=False)
timestr    = pd.to_datetime(timestr, format="%b %d %Y")
range_loop = len(timestr)
bar        = IncrementalBar('', max=len(timestr))

# Create a variable with zeros to list the data throught loop.
tc_list    = np.zeros([range_loop])
rh_list    = np.zeros([range_loop])
slp_list   = np.zeros([range_loop])
uvmet_list = np.zeros([range_loop])

# Starting looping throught time.
for i in range(inittime,range_loop,1):
    # Open WRF variables.
    rh         = getvar(nc_file,'rh2',timeidx=i,method='cat',squeeze=True,meta=False)
    uvmet      = getvar(nc_file,'uvmet10_wspd_wdir',timeidx=i,method='cat',squeeze=True,meta=False,units="m s-1")
    slp        = getvar(nc_file,'slp',timeidx=i,method='cat',squeeze=True,meta=False,units="hPa")
    tc         = nc_file.variables['T2'][i,:,:]
    
    # Calculate the XY location throught the specified latitude and longitude coordinates.
    latlon     = ll_to_xy(nc_file,lat,lon,timeidx=i,squeeze=True,meta=False,as_int=True,stagger='m')
    
    # Store the data from the desired location.
    rh            = rh[latlon[0],latlon[1]]
    uvmet         = uvmet[0,latlon[0],latlon[1]]
    slp           = slp[latlon[0],latlon[1]]
    tc            = tc[latlon[0],latlon[1]]
    tc            = tc-273.15
    tc_list[i]    = round(tc,1)
    rh_list[i]    = round(rh,0)
    slp_list[i]   = round(slp,1)
    uvmet_list[i] = round(uvmet,1)
    
    # Update the IncrementalBar.
    bar.next()
bar.finish()     

# Store the lists inside one variable.
rows = zip(timestr, tc_list,rh_list,slp_list,uvmet_list)

# Save the variables in a CSV file
with open('./output.csv', "w") as f:
    writer = csv.writer(f)
    for row in rows:
        writer.writerow(row)
        