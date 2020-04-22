"""
File name:      xy_lcamd.py
Author:         Ueslei Adriano Sutil
Email:          uesleisutil1@gmail.com
Created:        07 Apr 2020
Last modified:  22 Apr 2020
Version:        2.0

Creates XY graph from CSV file.
"""

import pandas as pd
import numpy as np
from matplotlib.ticker   import AutoMinorLocator
import matplotlib.pyplot as plt
import matplotlib.dates  as mdates
import matplotlib
import netCDF4
import scipy.stats
matplotlib.use('Agg')
import warnings
warnings.filterwarnings("ignore")

# MAtplotlib resources.
plt.rcParams.update({
    'font.family':'sans-serif',
    'font.sans-serif':['Liberation Sans'],
    'font.size':'12',
    })

# Open CSV file and import variables.
csv_import = pd.read_csv('/media/ueslei/Ueslei/INPE/Doutorado/LCAMD/v02/12_04_2020/data.txt',
                          keep_default_na=True, delimiter=',', header=None,
                          names=['Temp_WD','Pressure_BME'])

temp_wd    = csv_import['Temp_WD']
temp_wd    = temp_wd.values
pres_bme   = csv_import['Pressure_BME']
pres_bme   = pres_bme.values

csv_import = pd.read_csv('/media/ueslei/Ueslei/INPE/Doutorado/LCAMD/v02/12_04_2020/inmet_taubate.csv',
                          keep_default_na=True, delimiter=',', header=None,
                          names=['Temp','Pressure'])

temp_inmet = csv_import['Temp']
temp_inmet = temp_inmet.values
pres_inmet = csv_import['Pressure']
pres_inmet = pres_inmet.values-5

# Calculate correlation, trend and RMSE
# Variable 1.
stats = scipy.stats.spearmanr(temp_wd,temp_inmet, axis=0, nan_policy='propagate')
cor1  = str(stats[0])
pval1 = stats[1]
cor1  = cor1[0:5]
z1    = np.polyfit(temp_wd,temp_inmet, 1)
p1    = np.poly1d(z1)
rmse1 = np.sqrt(np.mean((temp_wd-temp_inmet)**2))
rmse1 = str(rmse1)
rmse1 = rmse1[0:4]

bias1 = temp_wd-temp_inmet

# Variable 2.
stats = scipy.stats.spearmanr(pres_bme,pres_inmet, axis=0, nan_policy='propagate')
cor2  = str(stats[0])
pval2 = stats[1]
cor2  = cor2[0:5]
z2    = np.polyfit(pres_bme,pres_inmet, 1)
p2    = np.poly1d(z2)

rmse2 = np.sqrt(np.mean((pres_bme-pres_inmet)**2))
rmse2 = str(rmse2)
rmse2 = rmse2[0:4]

bias2 = pres_bme-pres_inmet

# Rearrange time.
idx   = pd.date_range('2020-04-12 12:00', '2020-04-19 18:45', freq = '1H')
hours = mdates.HourLocator(interval = 6)
h_fmt = mdates.DateFormatter('%d/%M %Hh')

# Create figures subplots.
axs = plt.subplots(4, 1, constrained_layout=True)
fig = plt.figure()

# Y-axis range for both plots.
y_major_ticks  = np.arange(15, 32, 3)
y_major_ticks2 = np.arange(938, 953.1, 3)
y_major_ticks3 = np.arange(0, 1, 2)
y_major_ticks4 = np.arange(-4, 5, 0.5)

# Plot the figures.
# Plot 1.
ax1 = plt.subplot2grid((3,3), (0,0),rowspan=1,colspan=2)
ax1.plot(idx,temp_wd,linestyle='-',label='LCAMD', marker='',linewidth=2,color='C0')  
ax1.plot(idx,temp_inmet,linestyle='-',label='INMET',marker='',linewidth=2,color='crimson') 
#fill_between(idx, temp_wd, temp_inmet,alpha=0.3,color='gray')
ax1.legend(loc='upper left',fontsize=10)
ax1.xaxis.set_minor_locator(AutoMinorLocator())
ax1.yaxis.set_minor_locator(AutoMinorLocator())
ax1.tick_params(which='both',  width=2)
ax1.tick_params(which='major', length=7)
ax1.tick_params(which='minor', length=4)
ax1.set_yticks(y_major_ticks)
ax1.set_ylabel('Temperatura do Ar \n [°C]',fontsize=12, labelpad=10)
ax1.set_title('Teste #02 LCAMD v.2 \n 12/04/2020 12:00 a 19/04/2020 18:45')
for label in ax1.get_xmajorticklabels():
    label.set_rotation(30)
    label.set_horizontalalignment("right")
ax1.twinx()
plt.ylabel("Viés da Temperatura do Ar \n [°C]",color="k",rotation=270,labelpad=40),
plt.bar(idx,bias1,alpha=0.3, color='steelblue',width=0.035,align='center') 
plt.axhline(y=0,xmin=0.0454, xmax=0.954,linestyle='-',label='Viés',alpha=0.3,color='steelblue',linewidth=1)
plt.legend(loc='upper right',fontsize=10)

# Plot 2.
ax2 = plt.subplot2grid((3,3), (1,0), rowspan=1,colspan=2) 
ax2.plot(idx,pres_bme,label='LCAMD',linestyle='-',marker='',linewidth=2,color='C0',alpha=1)    
ax2.plot(idx,pres_inmet,label='INMET',linestyle='-',marker='',linewidth=2,color='crimson')
#ax2.fill_between(idx, pres_bme, pres_inmet,alpha=0.3,color='gray')
ax2.legend(loc='upper left',fontsize=10)
ax2.xaxis.set_minor_locator(AutoMinorLocator())
ax2.yaxis.set_minor_locator(AutoMinorLocator())
ax2.tick_params(which='both',  width=2)
ax2.tick_params(which='major', length=7)
ax2.tick_params(which='minor', length=4)
ax2.set_yticks(y_major_ticks2)
ax2.set_ylabel('Pressão Atmosférica \n [hPa]',fontsize=12, labelpad=10)
for label in ax2.get_xmajorticklabels():
    label.set_rotation(30)
    label.set_horizontalalignment("right")
ax2.twinx()
plt.ylabel("Viés da Pressão Atmosférica \n [hPa]",color="k",rotation=270,labelpad=25),
plt.bar(idx,bias2,alpha=0.3, color='steelblue',width=0.035,align='center') 
plt.axhline(y=0,xmin=0.0454, xmax=0.954,linestyle='-',label='Viés',alpha=0.3,color='steelblue',linewidth=1)
plt.legend(loc='upper right',fontsize=10)

# Plot 3.
ax3 = plt.subplot2grid((3,3), (2,0)) 
ax3.scatter(temp_wd,temp_inmet,color='C0')   
ax3.plot(temp_wd,p1(temp_wd),linestyle='-',color='crimson',linewidth=2)
ax3.xaxis.set_minor_locator(AutoMinorLocator())
ax3.yaxis.set_minor_locator(AutoMinorLocator())
ax3.tick_params(which='both',  width=2)
ax3.tick_params(which='major', length=7)
ax3.tick_params(which='minor', length=4)
ax3.set_yticks(y_major_ticks)
ax3.set_xticks(y_major_ticks)
ax3.set_ylabel('Temperatura do Ar \n INMET [°C]',fontsize=12, labelpad=10)
ax3.set_xlabel('Temperatura do Ar \n LCAMD [°C]',fontsize=12, labelpad=10)
ax3.text(26.4,14.2,'r² = '+cor1,color='k',fontsize=12)
if pval1 < 0.001:
    ax3.text(26.5,15.2,'p < 0.001',color='k',fontsize=12)
else:
    ax3.text(26.6,15.2,'p > 0.001',color='k',fontsize=12)  
ax3.text(25.9,16.2,'RMSE: '+rmse1,color='k',fontsize=12)

# plot 4.
ax4 = plt.subplot2grid((3,3), (2,1),)     
ax4.scatter(pres_bme,pres_inmet,color='C0')  
ax4.plot(pres_bme,p2(pres_bme),linestyle='-',color='crimson',linewidth=2)
ax4.xaxis.set_minor_locator(AutoMinorLocator())
ax4.yaxis.set_minor_locator(AutoMinorLocator())
ax4.tick_params(which='both',  width=2)
ax4.tick_params(which='major', length=7)
ax4.tick_params(which='minor', length=4)
ax4.set_yticks(y_major_ticks2)
ax4.set_xticks(y_major_ticks2)
ax4.set_ylabel('Pressão Atmosférica \n INMET [hPa]',fontsize=12, labelpad=10)
ax4.set_xlabel('Pressão Atmosférica \n LCAMD [hPa]',fontsize=12, labelpad=10)
ax4.text(948.8,936.9,'r² = '+cor2,color='k',fontsize=12)
if pval2 < 0.001:
    ax4.text(948.9,937.9,'p < 0.001',color='k',fontsize=12)
else:
    ax4.text(948.9,937.9,'p > 0.001',color='k',fontsize=12)
ax4.text(948.3,938.9,'RMSE: '+rmse2,color='k',fontsize=12)


# Set final fig resources and then save figure.
fig.set_size_inches(13, 13)
fig.tight_layout()
fig.subplots_adjust(hspace = 0.3)
fig.savefig('xy_lcamd.png',transparent=False, bbox_inches = 'tight', pad_inches=0, dpi=250)
