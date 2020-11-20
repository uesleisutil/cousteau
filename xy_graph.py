"""
File name:      xy_graph.py
Author:         Ueslei Adriano Sutil
Email:          uesleisutil1@gmail.com
Created:        19 November 2020
Last modified:  20 November 2020
"""

# Import libraries
import pandas as pd
import numpy as np
from matplotlib.ticker import AutoMinorLocator
import matplotlib.pyplot as plt

# Open CSV file and import variables.
file_blumenau = pd.read_csv('/media/ueslei/Ueslei_HD/INPE/PCI/Projetos/SC_2008/Dados/Tabelas/Blumenau.csv', 
                        keep_default_na=True,delimiter=',', header=None,names=['Simu','SST', 'Umid','Prec'])
file_major = pd.read_csv('/media/ueslei/Ueslei_HD/INPE/PCI/Projetos/SC_2008/Dados/Tabelas/Major.csv', 
                        keep_default_na=True,delimiter=',', header=None,names=['Simu','SST', 'Umid','Prec'])

# Import CSV data.
sst = file_blumenau['SST']
sst = sst.values

umid_blumenau = file_blumenau['Umid']
umid_blumenau = umid_blumenau.values
umid_blumenau_err = 2*np.std(umid_blumenau)

prec_blumenau = file_blumenau['Prec']
prec_blumenau = prec_blumenau.values


umid_major = file_major['Umid']
umid_major = umid_major.values
prec_major = file_major['Prec']
prec_major = prec_major.values

# Create figures plots.
fig, axs = plt.subplots(2, 2)
fig.subplots_adjust(wspace=0.25, hspace=0.25)

# Set X and Y ticks.
x_major_ticks = np.arange(21.5, 23.5, 0.5)
y_major_ticks_1 = np.arange(11, 12.1, 0.2)
y_major_ticks_2 = np.arange(20, 95, 10)

# Plot the figures. 
axs[0,0].plot(sst, umid_blumenau,marker='.', linewidth=1,markersize=12, color='C10')   
axs[0,0].plot(sst, umid_blumenau,linestyle='-', linewidth=1, color='black')    

axs[0,1].plot(sst, umid_major,marker='.', linewidth=1,markersize=12, color='C10')   
axs[0,1].plot(sst, umid_major,linestyle='-', linewidth=1, color='black')    

axs[1,0].plot(sst, prec_blumenau,marker='.', linewidth=1,markersize=12, color='C2')   
axs[1,0].plot(sst, prec_blumenau,linestyle='-', linewidth=1, color='black')    

axs[1,1].plot(sst, prec_major,marker='.', linewidth=1,markersize=12, color='C2')   
axs[1,1].plot(sst, prec_major,linestyle='-', linewidth=1, color='black')    

# Create nice axis.
axs[0,0].yaxis.set_minor_locator(AutoMinorLocator())
axs[0,0].tick_params(which='both',  width=2)
axs[0,0].tick_params(which='major', length=6)
axs[0,0].tick_params(which='minor', length=0)

axs[0,1].yaxis.set_minor_locator(AutoMinorLocator())
axs[0,1].tick_params(which='both',  width=2)
axs[0,1].tick_params(which='major', length=6)
axs[0,1].tick_params(which='minor', length=0)

axs[1,0].yaxis.set_minor_locator(AutoMinorLocator())
axs[1,0].tick_params(which='both',  width=2)
axs[1,0].tick_params(which='major', length=6)
axs[1,0].tick_params(which='minor', length=0)

axs[1,1].yaxis.set_minor_locator(AutoMinorLocator())
axs[1,1].tick_params(which='both',  width=2)
axs[1,1].tick_params(which='major', length=6)
axs[1,1].tick_params(which='minor', length=0)

# Import ticks.
axs[0,0].set_xticks(x_major_ticks)
axs[0,0].set_yticks(y_major_ticks_1)

axs[0,1].set_xticks(x_major_ticks)
axs[0,1].set_yticks(y_major_ticks_1)

axs[1,0].set_xticks(x_major_ticks)
axs[1,0].set_yticks(y_major_ticks_2)

axs[1,1].set_xticks(x_major_ticks)
axs[1,1].set_yticks(y_major_ticks_2)

# Set titles
axs[0,0].set_title('Blumenau',fontsize=9)
axs[0,1].set_title('Major Gercino',fontsize=9)

# Set axis labels.
axs[0,0].set_ylabel('Specific Humidity [g/kg]',fontsize=9)
axs[0,0].set_xlabel('Sea Surface Temperature [ᵒC] ',fontsize=9)
axs[0,0].tick_params(axis='x', labelsize=9)
axs[0,0].tick_params(axis='y', labelsize=9)

axs[0,1].set_ylabel('Specific Humidity [g/kg]',fontsize=9)
axs[0,1].set_xlabel('Sea Surface Temperature [ᵒC] ',fontsize=9)
axs[0,1].tick_params(axis='x', labelsize=9)
axs[0,1].tick_params(axis='y', labelsize=9)

axs[1,0].set_ylabel('Precipitation [mm]',fontsize=9)
axs[1,0].set_xlabel('Sea Surface Temperature [ᵒC] ',fontsize=9)
axs[1,0].tick_params(axis='x', labelsize=9)
axs[1,0].tick_params(axis='y', labelsize=9)

axs[1,1].set_ylabel('Precipitation [mm]',fontsize=9)
axs[1,1].set_xlabel('Sea Surface Temperature [ᵒC] ',fontsize=9)
axs[1,1].tick_params(axis='x', labelsize=9)
axs[1,1].tick_params(axis='y', labelsize=9)

# Create ref line.
axs[0,0].axvline(22.444, color='crimson', linestyle='--', lw=1)
axs[0,0].text(22.47,11,'Normal',color='crimson',fontsize=8)

axs[0,1].axvline(22.444, color='crimson', linestyle='--', lw=1)
axs[0,1].text(22.47,11,'Normal',color='crimson',fontsize=8)

axs[1,0].axvline(22.444, color='crimson', linestyle='--', lw=1)
axs[1,0].text(22.47,20,'Normal',color='crimson',fontsize=8)

axs[1,1].axvline(22.444, color='crimson', linestyle='--', lw=1)
axs[1,1].text(22.47,20,'Normal',color='crimson',fontsize=8)

# Save figure.
fig.set_size_inches(10,10)
fig.savefig('xy_graph_csv.png', bbox_inches = 'tight', pad_inches=0, dpi=125)
