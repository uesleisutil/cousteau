#!/home/uesleisutil/anaconda3/bin/python python
# -*- coding: utf-8 -*-

"""
File name:      taylor_target_diagrams.py
Author:         Ueslei Adriano Sutil
Email:          uesleisutil1@gmail.com
Created:        16 August 2019
Last modified:  16 August 2019
Version:        1.0
Python:         3.7.1

Extract data from ASCII files, save as pickle then create Taylor and Target Diagrams.
"""

import matplotlib.pyplot as plt
from matplotlib import rcParams
import skill_metrics as sm
from sys import version_info
import pandas as pd
import numpy as np
import matplotlib.gridspec as gridspec

# Open CSV file and import variables.
data1     = pd.read_csv('/media/ueslei/Ueslei/INPE/PCI/Projetos/SC_2008/Dados/In_situ_data/INMET/Paranagua_83844.csv', keep_default_na=True, delimiter=',', header=None,names=['T_EM','UR_EM','SLP_EM','VVel_EM','T_WRF','UR_WRF','SLP_WRF', 'VVel_WRF'])
vvel_wrf1 = data1['VVel_WRF']
vvel_wrf1 = vvel_wrf1.values
slp_wrf1  = data1['SLP_WRF']
slp_wrf1  = slp_wrf1.values
ur_wrf1   = data1['UR_WRF']
ur_wrf1   = ur_wrf1.values
t_wrf1    = data1['T_WRF']
t_wrf1    = t_wrf1.values
vvel_em1  = data1['VVel_EM']
vvel_em1  = vvel_em1.values
slp_em1   = data1['SLP_EM']
slp_em1   = slp_em1.values
ur_em1    = data1['UR_EM']
ur_em1    = ur_em1.values
t_em1     = data1['T_EM']
t_em1     = t_em1.values

# Start estatistics.
taylor_stats1 = sm.taylor_statistics(t_wrf1,t_em1,'data')
taylor_stats2 = sm.taylor_statistics(ur_wrf1,ur_em1,'data')
taylor_stats3 = sm.taylor_statistics(slp_wrf1,slp_em1,'data')
taylor_stats4 = sm.taylor_statistics(vvel_wrf1,vvel_em1,'data')

target_stats1 = sm.target_statistics(t_wrf1,t_em1,'data')
target_stats2 = sm.target_statistics(ur_wrf1,ur_em1,'data')
target_stats3 = sm.target_statistics(slp_wrf1,slp_em1,'data')
target_stats4 = sm.target_statistics(vvel_wrf1,vvel_em1,'data')

# Store estatistics in arrays.
crmsd1 = np.array([taylor_stats1['crmsd'][1], taylor_stats2['crmsd'][1], taylor_stats3['crmsd'][1], taylor_stats4['crmsd'][1]])
ccoef1 = np.array([taylor_stats1['ccoef'][1], taylor_stats2['ccoef'][1], taylor_stats3['ccoef'][1], taylor_stats4['ccoef'][1]])
sdev1  = np.array([taylor_stats1['sdev'][1], taylor_stats2['sdev'][1], taylor_stats3['sdev'][1], taylor_stats4['sdev'][1]])

bias   = np.array([target_stats1['bias'], target_stats2['bias'], target_stats3['bias'], target_stats4['bias']])
crmsd  = np.array([target_stats1['crmsd'], target_stats2['crmsd'], target_stats3['crmsd'],target_stats4['crmsd']])
rmsd   = np.array([target_stats1['rmsd'], target_stats2['rmsd'], target_stats3['rmsd'],target_stats4['rmsd']])


label = ['Air temperature at 2 meters', 'Relative humidity at 2 meters', 'Sea level pressure', 'Wind Velocity at 10 meters']


fig = plt.figure(1,figsize=(15, 15))
fig.set_figwidth(15) # optional setting the width of the image
fig.add_subplot(221)
sm.taylor_diagram(sdev1,crmsd1,ccoef1, numberpanels=2, alpha = 0.0, markerLabel = label, markerSize=8,markerColor='b',markerLegend='on', colOBS = 'crimson', colRMS='C0', colCOR='k',markerobs = '*', styleOBS='-',styleCOR='-', styleRMS=':', titleCOR='on',titleSTD='on',titleRMS='on',widthOBS=0.7,widthCOR=0.2,widthRMS=1,widthSTD=0.7,tickSTD=[2,4,6,0,-2,-4,-6],tickCOR=[1, .99, .95, .90, .8, .6, .3, -.3, -.6, -.8, -.9, -.95, -.99, -1])
fig.add_subplot(233)
sm.target_diagram(bias,crmsd,rmsd,markerLabel=label, markersize=8,ticks=np.arange(20,20,2),circleLineWidth=0.8,markerLabelColor = 'b',markerColor = 'b', markerLegend = 'on')
fig.subplots_adjust(left=0, bottom=0, right=1, top=1,wspace=0, hspace=-1)
fig.tight_layout()
#plt.subplot(fig1,fig2)
plt.savefig('taylor2.png', transparent=False, bbox_inches = 'tight', pad_inches=0, dpi=250)


# Write plot to file
#plt.savefig('taylor2.png', transparent=False, bbox_inches = 'tight', pad_inches=0, dpi=250)


