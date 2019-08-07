#!/home/uesleisutil/anaconda3/bin/python python
# -*- coding: utf-8 -*-

"""
File name:      xy_csv_plot.py
Author:         Ueslei Adriano Sutil
Email:          uesleisutil1@gmail.com
Created:        11 July 2019
Last modified:  12 July 2019
Version:        1.0
Python:         3.7.1

Creates XY graph from CSV file.
"""

import pandas as pd
import numpy as np
from matplotlib.ticker import AutoMinorLocator
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
# Open CSV file and import variables.
df_blu         = pd.read_csv('/media/ueslei/Ueslei/INPE/PCI/Projetos/SC_2008/Dados/In_situ_data/blumenau/blumenau.csv', keep_default_na=True, delimiter=',', header=None,names=['Year','Total', 'DayMax'])
time_blu       = df_blu['Year']
time_blu       = time_blu.values
max_rain_blu   = df_blu['DayMax']
max_rain_blu   = max_rain_blu.values
total_rain_blu = df_blu['Total']
total_rain_blu = total_rain_blu.values

df_maj         = pd.read_csv('/media/ueslei/Ueslei/INPE/PCI/Projetos/SC_2008/Dados/In_situ_data/major/major.csv', keep_default_na=True, delimiter=',', header=None,names=['Year','Total', 'DayMax'])
time_maj       = df_maj['Year']
time_maj       = time_maj.values
max_rain_maj   = df_maj['DayMax']
max_rain_maj   = max_rain_maj.values
total_rain_maj = df_maj['Total']
total_rain_maj = total_rain_maj.values

# Create figures subplots.
fig, axs = plt.subplots(2, 1)
fig.subplots_adjust(wspace=0, hspace=0.12)
# Y-axis range for both plots.
y_major_ticks  = np.arange(0, 1001, 100)
y_major_ticks2 = np.arange(0, 1001, 100)
x_major_ticks  = np.arange(1945, 2011, 5)

# Plot the figures.
axs[0].plot(time_blu,max_rain_blu,label='November total precipitation (mm)',linestyle='-',marker='.',linewidth=2,color='C0')      
axs[0].plot(time_blu,total_rain_blu,label='November daily maximum precipitation (mm)',linestyle='-',marker='.',linewidth=2,color='crimson')
axs[0].legend(loc='upper left',fontsize=8)
axs[0].yaxis.set_minor_locator(AutoMinorLocator())
axs[0].tick_params(which='both',  width=2)
axs[0].tick_params(which='major', length=7)
axs[0].tick_params(which='minor', length=4)
axs[0].xaxis.set_visible(False)
axs[0].set_yticks(y_major_ticks)
axs[0].set_ylabel('Precipitation [mm]',fontsize=9)
axs[0].set_title('November precipitation history: Blumenau')
axs[0].axvline(2008, color='k', linestyle=':', lw=1)
axs[0].text(2002,980,'1001.2 mm',color='C0')
axs[0].text(2002,250,'250.9 mm',color='crimson')


axs[1].plot(time_maj,max_rain_maj,label='November total precipitation (mm)',linestyle='-',marker='.',linewidth=2,color='C0')      
axs[1].plot(time_maj,total_rain_maj,label='November daily maximum precipitation (mm)',linestyle='-',marker='.',linewidth=2,color='crimson')
axs[1].set_xlabel('Years')
axs[1].set_xticks(x_major_ticks)
axs[1].legend(loc='upper left',fontsize=8)
axs[1].xaxis.set_minor_locator(AutoMinorLocator())
axs[1].yaxis.set_minor_locator(AutoMinorLocator())
axs[1].tick_params(which='both',  width=2)
axs[1].tick_params(which='major', length=7)
axs[1].tick_params(which='minor', length=4)
axs[1].set_yticks(y_major_ticks)
axs[1].set_ylabel('Precipitation [mm]',fontsize=9)
axs[1].set_title('November precipitation history: Major Gercino')
axs[1].axvline(2008, color='k', linestyle=':', lw=1)
axs[1].text(2002,860,'860.6.2 mm',color='C0')
axs[1].text(2002,300,'300.5 mm',color='crimson')


# Save figure.
fig.set_size_inches(13, 7)
fig.savefig('xy_prec_plot.png',transparent=False, bbox_inches = 'tight', pad_inches=0, dpi=250)