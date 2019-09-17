#!/home/uesleisutil/anaconda3/bin/python python
# -*- coding: utf-8 -*-

"""
File name:      xy_sst_plot.py
Author:         Ueslei Adriano Sutil
Email:          uesleisutil1@gmail.com
Created:        10 July 2019
Last modified:  17 September 2019
Version:        1.5
Python:         3.6.8

Creates XY graph from Global Ocean Physics Reanalisys (GLORYS12V1; von Schuckmann et al. 2016),
Simple Ocean Data Assimilation 3 (SODA3; Carton et al., 2018), Optimum Interpolation Sea Surface
Temperature V2 High Resolution Dataset (OISST; Reynolds et al., 2007) or 

GLORYS12V1:
Overview: https://www.mercator-ocean.fr/wp-content/uploads/2017/06/FS-GLORYS2V4_EN-2.pdf
von Schuckmann, K., and Coauthors, 2016: The Copernicus marine environment monitoring service ocean state report.
J. Oper. Oceanogr., 9 (Suppl.), S235–S320, https://doi.org/10.1080/1755876X.2016.1273446.

SODA 3.12.2:
Carton, J.A., G.A. Chepurin, and L. Chen, 2018: SODA3: a new ocean climate reanalysis, J. Clim., 31, 6967-6983,
DOI:10.1175/JCLI-D-18-0149.1. 

OISST V2 HRD:
Reynolds, Richard W., Thomas M. Smith, Chunying Liu, Dudley B. Chelton, Kenneth S. Casey, Michael G. Schlax, 2007: 
Daily High-Resolution-Blended Analyses for Sea Surface Temperature. J. Climate, 20, 5473-5496. 

!!! IMPORTANT !!!
To calculate the anomaly, you need to process Glorys data using CDO. What I do:
cdo -b F64 -ymonsub glorys_d03_1993_2015_montlhy.nc  -ymonmean glorys_d03_1993_2015_montlhy.nc glorys_d03_1993_2015_montlhy_anon.nc

To calculate the anomaly, you need to process SODA data using CDO. What I do:
cdo sellonlatbox,-56.5993,-42.42,-21.8274,-34.3007 soda_1980_2016_montlhy.nc soda_d02_1980_2016_montlhy.nc
cdo sellevel,5.03355 soda_d02_1980_2016_montlhy.nc soda_sst_d02_1980_2016_montlhy.nc
cdo -b F64 -ymonsub soda_sst_d02_1980_2016_montlhy.nc  -ymonmean soda_sst_d02_1980_2016_montlhy.nc soda_sst_d02_1980_2016_anom.nc

To calculate the anomaly, you need to process ERA5 data using CDO. What I do:
cdo -b F64 -ymonsub era5.nc  -ymonmean era5.nc era5_ssta.nc

"""

import xarray as xr
import os
import pandas as pd
import numpy as np
from matplotlib.ticker import AutoMinorLocator
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from sklearn.linear_model import LinearRegression
from roms_libs import *
from netCDF4 import Dataset
matplotlib.use('Agg')

# Open files.
print('Choose dataset: (1) GLORYS2V4, (2) SODA3, (3) OISST or (4) ERA5:')
dataset    = input()
if dataset == '1':
    file_dir   = os.path.expanduser('/media/ueslei/Ueslei/INPE/PCI/Projetos/SC_2008/Dados/Glorys2V4/')
    file_name  = os.path.join(file_dir,'glorys_d01_1993_2015_monthly.nc')
    file_name2 = os.path.join(file_dir,'glorys_d01_1993_2015_monthly_anom.nc')
if dataset == '2':
    file_dir   = os.path.expanduser('/media/ueslei/Ueslei/INPE/PCI/Projetos/SC_2008/Dados/SODA/')
    file_name  = os.path.join(file_dir,'soda_sst_d01_1980_2016_monthly.nc')
    file_name2 = os.path.join(file_dir,'soda_sst_d01_1980_2016_monthly_anom.nc')
if dataset == '3':
    file_dir   = os.path.expanduser ('/media/ueslei/Ueslei/INPE/PCI/Projetos/SC_2008/Dados/OISST/')
    file_name  = os.path.join(file_dir,'oisst_d01_monmean.nc')
    file_name2 = os.path.join(file_dir,'oisst_ssta_d01_monmean.nc')    
    lonbounds  = [-65,-35]
    latbounds  = [-40,-15]
if dataset == '4':
    file_dir   = os.path.expanduser ('/media/ueslei/Ueslei/INPE/PCI/Projetos/SC_2008/Dados/ERA5/')
    file_name  = os.path.join(file_dir,'era5_montlhy.nc')
    file_name2 = os.path.join(file_dir,'era5_ssta_montlhy.nc')    
    lonbounds  = [-65,-35]
    latbounds  = [-40,-15]    

# Open variables.
if dataset =='1':
    sst_data  = xr.open_dataset(file_name)
    ssta_data = xr.open_dataset(file_name2)
    sst       = sst_data.temperature.data-273.15
    ssta      = ssta_data.temperature.data
if dataset =='2':
    sst_data  = xr.open_dataset(file_name)
    ssta_data = xr.open_dataset(file_name2)
    sst       = sst_data.temp.data
    ssta      = ssta_data.temp.data   
if dataset =='3':
    sst_data  = netCDF4.Dataset(file_name)
    ssta_data =netCDF4.Dataset(file_name2)
    lats      = sst_data.variables['lat'][:]
    lons      = sst_data.variables['lon'][:]-360
    latli     = np.argmin(np.abs(lats-latbounds[0]))
    latui     = np.argmin(np.abs(lats-latbounds[1])) 
    lonli     = np.argmin(np.abs(lons-lonbounds[0]))
    lonui     = np.argmin(np.abs(lons-lonbounds[1])) 
    sst       = sst_data.variables['sst']
    sst       = sst[:,latli:latui,lonli:lonui] 
    ssta      = ssta_data.variables['anom']
    ssta      = ssta[:,latli:latui,lonli:lonui]
if dataset =='4':
    sst_data  = netCDF4.Dataset(file_name)
    ssta_data =netCDF4.Dataset(file_name2)
    lats      = sst_data.variables['latitude'][:]*-1
    lons      = sst_data.variables['longitude'][:]-180
    latli     = np.argmin(np.abs(lats-latbounds[0]))
    latui     = np.argmin(np.abs(lats-latbounds[1])) 
    lonli     = np.argmin(np.abs(lons-lonbounds[0]))
    lonui     = np.argmin(np.abs(lons-lonbounds[1])) 
    sst       = sst_data.variables['sst']
    sst       = sst[:,latli:latui,lonli:lonui]-273.15
    ssta      = ssta_data.variables['sst']
    ssta      = ssta[:,latli:latui,lonli:lonui]
    
# Calculate means.
if dataset == '1' or dataset == '2':
    sst_mean  = np.nanmean(np.nanmean(sst,axis=(1,2)),axis=1,dtype='float64')
    ssta_mean = np.nanmean(np.nanmean(ssta,axis=(1,2)),axis=1,dtype='float64')
if dataset == '3' or dataset == '4':
    sst_mean   = np.nanmean(sst,axis=(1,2),dtype='float64')
    ssta_mean  = np.nanmean(ssta,axis=(1,2),dtype='float64')

# Create X-axis time range.
if dataset == '1':
    dates = pd.date_range(start='19930101', end='20151216',freq='MS')
if dataset == '2':
    dates = pd.date_range(start='19800101', end='20101216',freq='MS') 
if dataset == '3':
    dates = pd.date_range(start='19810901', end='20121231',freq='MS')   
if dataset == '4':
    dates = pd.date_range(start='19790101', end='20121201',freq='MS')   

# Create DataFrame using Pandas.
sst_df  = pd.DataFrame({'SST':sst_mean},index=dates,columns=['SST'])
ssta_df = pd.DataFrame({'SSTA':ssta_mean},index=dates,columns=['SSTA','Trend line'])

# Create figures subplots.
fig = plt.figure()
fig.subplots_adjust(wspace=0, hspace=0.15)
ax1 = fig.add_subplot(211)

# Y-axis range for both plots.
if dataset == '1':
    y_major_ticks  = np.arange(14, 25, 1)
    y_major_ticks2 = np.arange(-2, 2, 0.5)
if dataset == '2':
    y_major_ticks  = np.arange(14, 25, 1)
    y_major_ticks2 = np.arange(-2, 2, 0.5)
if dataset == '3':
    y_major_ticks  = np.arange(14, 25, 1)
    y_major_ticks2 = np.arange(-2, 2, 0.5)
if dataset == '4':
    y_major_ticks  = np.arange(19, 35, 2)
    y_major_ticks2 = np.arange(-2, 2, 0.5)

# Options for the first subplot.
ax1.set_ylabel('Sea Surface Temperature [$^\circ\!$C] \n',fontsize=9)
ax1.set_xlabel('Years')
ax1.set_yticks(y_major_ticks)
ax1.xaxis.set_minor_locator(AutoMinorLocator())
ax1.yaxis.set_minor_locator(AutoMinorLocator())
ax1.tick_params(which='both',  width=2)
ax1.tick_params(which='major', length=7)
ax1.tick_params(which='minor', length=4)
ax1.xaxis.set_visible(False)
if dataset == '1':
    ax1.set_title('Monthly Sea Surface Temperature Mean [$^\circ\!$C] from GLORYS2V4',fontsize=10)
if dataset == '2':
    ax1.set_title('Monthly Sea Surface Temperature Mean [$^\circ\!$C] from SODA3',fontsize=10)
if dataset == '3':
    ax1.set_title('Monthly Sea Surface Temperature Mean [$^\circ\!$C] from OISST V2 HRD',fontsize=10)
if dataset == '4':
    ax1.set_title('Monthly Sea Surface Temperature Mean [$^\circ\!$C] from ERA5',fontsize=10)

# Options for the second subplot.
ax2   = fig.add_subplot(212)
ax2.set_ylabel('Sea Surface Temperature Anomaly [$^\circ\!$C]', fontsize=9)
ax2.set_xlabel('Year', color='k',fontsize=9)
ax2.set_yticks(y_major_ticks2)
ax2.yaxis.set_minor_locator(AutoMinorLocator())
ax2.tick_params(which='both',  width=2)
ax2.tick_params(which='major', length=7)
ax2.tick_params(which='minor', length=4)
ax2.xaxis.set_major_locator(mdates.MonthLocator(bymonth =(1), bymonthday=15))
ax2.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
if dataset == '1':
    ax2.set_title('Monthly Sea Surface Temperature Anomaly [$^\circ\!$C] from GLORYS2V4',fontsize=10)
if dataset == '2':
    ax2.set_title('Monthly Sea Surface Temperature Anomaly [$^\circ\!$C] from SODA3',fontsize=10)
if dataset == '3':
    ax2.set_title('Monthly Sea Surface Temperature Anomaly [$^\circ\!$C] from OISST V2 HRD',fontsize=10)    
if dataset == '4':
    ax2.set_title('Monthly Sea Surface Temperature Anomaly [$^\circ\!$C] from ERA5',fontsize=10)  
colors = ['crimson','C0']

# Add a trend line.
x       = np.arange(dates.size)
fit     = np.polyfit(x, ssta_df['SSTA'], 1)
fit_fn  = np.poly1d(fit)
ssta_tl = pd.DataFrame({'Trend line':fit_fn(x)},index=dates,columns=['Trend line'])

# Plot the figures.
sst_df.plot.line(ax=ax1,color='k',linestyle='-',linewidth=2,legend=False,x_compat=True)      
ssta_tl.plot.line(ax=ax2,linestyle='-',linewidth=1.7,legend=False,rot=0,color='Black')
ssta_df.plot.line(ax=ax2,linestyle='-',linewidth=0.7,legend=False,rot=0,color='Black')
ax2.fill_between(ssta_df.index, 0, ssta_df['SSTA'], where=ssta_df['SSTA'] >= 0, facecolor='crimson', interpolate=True,alpha=1)
ax2.fill_between(ssta_df.index, 0, ssta_df['SSTA'], where=ssta_df['SSTA'] <= 0, facecolor='C0', interpolate=True,alpha=1)

# Insert reference line and text to the figure.
ax1.axvline(pd.to_datetime('2008-11-01'), color='k', linestyle=':', lw=1)
ax2.axhline(y=0,linestyle='-',color='k',linewidth=1)
ax2.axvline(pd.to_datetime('2008-11-01'), color='k', linestyle=':', lw=1)

if dataset=='1':
    ax1.text(pd.to_datetime('2008-11-01'),14,'November 2008: 18.7 $^\circ\!$C',color='Black',fontsize=9)
    ax2.text(pd.to_datetime('2008-11-01'),-1.08,'November 2008: 0.3 $^\circ\!$C',color='crimson',fontsize=9)
if dataset=='2':
    ax1.text(pd.to_datetime('2003-11-01'),14,'November 2008: 18.9 $^\circ\!$C',color='k',fontsize=9)
    ax2.text(pd.to_datetime('2003-11-01'),-1.5,'November 2008: 0.35 $^\circ\!$C',color='crimson',fontsize=9)
if dataset=='3':
    ax1.text(pd.to_datetime('2003-11-01'),14,'November 2008: 18.9 $^\circ\!$C',color='k',fontsize=9)
    ax2.text(pd.to_datetime('2003-11-01'),-1.5,'November 2008: 0.35 $^\circ\!$C',color='crimson',fontsize=9)
if dataset=='4':
    ax1.text(pd.to_datetime('2003-06-01'),18.3,'November 2008: 24.7 $^\circ\!$C',color='k',fontsize=9)
    ax2.text(pd.to_datetime('2003-06-01'),-1.1,'November 2008: 0.73 $^\circ\!$C',color='crimson',fontsize=9)

# Save figure.
fig.set_size_inches(13, 7)
if dataset=='1':
    fig.savefig('xy_glorys_plot.png',transparent=False, bbox_inches = 'tight', pad_inches=0, dpi=250)
if dataset=='2':
    fig.savefig('xy_soda_plot.png',transparent=False, bbox_inches = 'tight', pad_inches=0, dpi=250)
if dataset=='3':
    fig.savefig('xy_oisst_plot.png',transparent=False, bbox_inches = 'tight', pad_inches=0, dpi=250)
if dataset=='4':
    fig.savefig('xy_era5_plot.png',transparent=False, bbox_inches = 'tight', pad_inches=0, dpi=250)