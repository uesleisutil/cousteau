#!/home/uesleisutil/anaconda3/bin/python python
# -*- coding: utf-8 -*-

"""
File name:      taylor_target_diagrams.py
Author:         Ueslei Adriano Sutil
Email:          uesleisutil1@gmail.com
Created:        16 August 2019
Last modified:  16 August 2019
Version:        2.0
Python:         3.7.1

Extract data from ASCII files, save as pickle then create Taylor and Target Diagrams.
"""

import matplotlib.pyplot as plt
import skill_metrics as sm
import pandas as pd
import numpy as np

print('Choose data: (1) Paranaguá, (2) Indaial, (3) Florianópolis or (4) São Joaquim')
choose_city = input()

# Open CSV file and import variables.
if choose_city == '1':
    title  = 'Paranaguá'
    latlon = '-25.53 °S and -48.51 °W'
    data   = pd.read_csv('/media/ueslei/Ueslei/INPE/PCI/Projetos/SC_2008/Dados/In_situ_data/INMET/Paranagua_83844.csv', keep_default_na=True, delimiter=',', header=None,names=['T_EM','UR_EM','SLP_EM','VVel_EM','T_WRF','UR_WRF','SLP_WRF', 'VVel_WRF'])
if choose_city == '2':
    title  = 'Indaial'
    latlon = '-26.90 °S and -49.21 °W'
    data   = pd.read_csv('/media/ueslei/Ueslei/INPE/PCI/Projetos/SC_2008/Dados/In_situ_data/INMET/Indaial_83872.csv', keep_default_na=True, delimiter=',', header=None,names=['T_EM','UR_EM','SLP_EM','VVel_EM','T_WRF','UR_WRF','SLP_WRF', 'VVel_WRF'])
if choose_city == '3':
    title  = 'Florianópolis'
    latlon = '-27.58 °S and -48.56 °W'
    data   = pd.read_csv('/media/ueslei/Ueslei/INPE/PCI/Projetos/SC_2008/Dados/In_situ_data/INMET/Floripa_83897.csv', keep_default_na=True, delimiter=',', header=None,names=['T_EM','UR_EM','SLP_EM','VVel_EM','T_WRF','UR_WRF','SLP_WRF', 'VVel_WRF'])
if choose_city == '4':
    title  ='São Joaquim'
    latlon = '-28.30 °S and -49.93 °W'
    data   = pd.read_csv('/media/ueslei/Ueslei/INPE/PCI/Projetos/SC_2008/Dados/In_situ_data/INMET/SaoJoaquim_83920.csv', keep_default_na=True, delimiter=',', header=None,names=['T_EM','UR_EM','SLP_EM','VVel_EM','T_WRF','UR_WRF','SLP_WRF', 'VVel_WRF'])

print('Choose variable: (1) Temperature, (2) Relative Humidity, (3) Sea Level Pressure or (4) Wind Velocity:')
choose_var  = input()
if choose_var == '1':
    var = 'Temperature [°C]'
if choose_var == '2':
    var = 'Relative Humidity [%]'
if choose_var == '3':
    var = 'Sea Level Pressure [hPa]'
if choose_var == '4':
    var = 'Wind Velocity [m s-1]'   

# First Data.
if choose_city == '1' and choose_var=='1':
    t_insitu     = data['T_EM']
    t_insitu     = t_insitu.values 
    t_model      = data['T_WRF']
    t_model      = t_model.values
    taylor_stats = sm.taylor_statistics(t_model,t_insitu)
    target_stats = sm.target_statistics(t_model,t_insitu)
if choose_city == '1' and choose_var=='2':
    ur_insitu    = data['UR_EM']
    ur_insitu    = ur_insitu.values
    ur_model     = data['UR_WRF']
    ur_model     = ur_model.values
    taylor_stats = sm.taylor_statistics(ur_model,ur_insitu)
    target_stats = sm.target_statistics(ur_model,ur_insitu)
if choose_city == '1' and choose_var=='3':
    slp_insitu   = data['SLP_EM']
    slp_insitu   = slp_insitu.values
    slp_model    = data['SLP_WRF']
    slp_model    = slp_model.values
    taylor_stats = sm.taylor_statistics(slp_model,slp_insitu)
    target_stats = sm.target_statistics(slp_model,slp_insitu)  
if choose_city == '1' and choose_var=='4':
    vvel_insitu  = data['VVel_EM']
    vvel_insitu  = vvel_insitu.values
    vvel_model   = data['VVel_WRF']
    vvel_model   = vvel_model.values
    taylor_stats = sm.taylor_statistics(vvel_model,vvel_insitu)
    target_stats = sm.target_statistics(vvel_model,vvel_insitu)  

# Second data.
if choose_city == '2' and choose_var=='1':
    t_insitu     = data['T_EM']
    t_insitu     = t_insitu.values 
    t_model      = data['T_WRF']
    t_model      = t_model.values
    taylor_stats = sm.taylor_statistics(t_model,t_insitu)
    target_stats = sm.target_statistics(t_model,t_insitu)
if choose_city == '2' and choose_var=='2':
    ur_insitu    = data['UR_EM']
    ur_insitu    = ur_insitu.values
    ur_model     = data['UR_WRF']
    ur_model     = ur_model.values
    taylor_stats = sm.taylor_statistics(ur_model,ur_insitu)
    target_stats = sm.target_statistics(ur_model,ur_insitu)
if choose_city == '2' and choose_var=='3':
    slp_insitu   = data['SLP_EM']
    slp_insitu   = slp_insitu.values
    slp_model    = data['SLP_WRF']
    slp_model    = slp_model.values
    taylor_stats = sm.taylor_statistics(slp_model,slp_insitu)
    target_stats = sm.target_statistics(slp_model,slp_insitu)  
if choose_city == '2' and choose_var=='4':
    vvel_insitu  = data['VVel_EM']
    vvel_insitu  = vvel_insitu.values
    vvel_model   = data['VVel_WRF']
    vvel_model   = vvel_model.values
    taylor_stats = sm.taylor_statistics(vvel_model,vvel_insitu)
    target_stats = sm.target_statistics(vvel_model,vvel_insitu)  

# Third data.
if choose_city == '3' and choose_var=='1':
    t_insitu     = data['T_EM']
    t_insitu     = t_insitu.values 
    t_model      = data['T_WRF']
    t_model      = t_model.values
    taylor_stats = sm.taylor_statistics(t_model,t_insitu)
    target_stats = sm.target_statistics(t_model,t_insitu)
if choose_city == '3' and choose_var=='2':
    ur_insitu    = data['UR_EM']
    ur_insitu    = ur_insitu.values
    ur_model     = data['UR_WRF']
    ur_model     = ur_model.values
    taylor_stats = sm.taylor_statistics(ur_model,ur_insitu)
    target_stats = sm.target_statistics(ur_model,ur_insitu)
if choose_city == '3' and choose_var=='3':
    slp_insitu   = data['SLP_EM']
    slp_insitu   = slp_insitu.values
    slp_model    = data['SLP_WRF']
    slp_model    = slp_model.values
    taylor_stats = sm.taylor_statistics(slp_model,slp_insitu)
    target_stats = sm.target_statistics(slp_model,slp_insitu)  
if choose_city == '3' and choose_var=='4':
    vvel_insitu  = data['VVel_EM']
    vvel_insitu  = vvel_insitu.values
    vvel_model   = data['VVel_WRF']
    vvel_model   = vvel_model.values
    taylor_stats = sm.taylor_statistics(vvel_model,vvel_insitu)
    target_stats = sm.target_statistics(vvel_model,vvel_insitu)  

# Fourth data.
if choose_city == '4' and choose_var=='1':
    t_insitu     = data['T_EM']
    t_insitu     = t_insitu.values 
    t_model      = data['T_WRF']
    t_model      = t_model.values
    taylor_stats = sm.taylor_statistics(t_model,t_insitu)
    target_stats = sm.target_statistics(t_model,t_insitu)
if choose_city == '4' and choose_var=='2':
    ur_insitu    = data['UR_EM']
    ur_insitu    = ur_insitu.values
    ur_model     = data['UR_WRF']
    ur_model     = ur_model.values
    taylor_stats = sm.taylor_statistics(ur_model,ur_insitu)
    target_stats = sm.target_statistics(ur_model,ur_insitu)
if choose_city == '4' and choose_var=='3':
    slp_insitu   = data['SLP_EM']
    slp_insitu   = slp_insitu.values
    slp_model    = data['SLP_WRF']
    slp_model    = slp_model.values
    taylor_stats = sm.taylor_statistics(slp_model,slp_insitu)
    target_stats = sm.target_statistics(slp_model,slp_insitu)  
if choose_city == '4' and choose_var=='4':
    vvel_insitu  = data['VVel_EM']
    vvel_insitu  = vvel_insitu.values
    vvel_model   = data['VVel_WRF']
    vvel_model   = vvel_model.values
    taylor_stats = sm.taylor_statistics(vvel_model,vvel_insitu)
    target_stats = sm.target_statistics(vvel_model,vvel_insitu)  

# Taylor Diagram arrays.
sdev  = np.array([taylor_stats['sdev'][0],taylor_stats['sdev'][1]])
crmsd = np.array([taylor_stats['crmsd'][1],taylor_stats['crmsd'][1]])
ccoef = np.array([taylor_stats['ccoef'][1],taylor_stats['ccoef'][1]])

# Target diagram arrays.
bias  = np.array([target_stats['bias'], target_stats['bias']])
crmsd = np.array([target_stats['crmsd'], target_stats['crmsd']])
rmsd  = np.array([target_stats['rmsd'], target_stats['rmsd']])


# Plot Taylor diagram.
fig = plt.figure(1,figsize=(15, 11))
fig.set_figwidth(15)
fig.add_subplot(221)
plt.title(title +': ' + var +'\n'+ latlon, pad=60)
sm.taylor_diagram(sdev,crmsd,ccoef, numberpanels=2, alpha = 0.0, checkstats = 'on', markerSize=12,markerColor='r',markerLegend='off', colOBS = 'crimson', colRMS='C0', colCOR='k',markerobs = '*', styleOBS='-',styleCOR='-', styleRMS=':', titleCOR='on',titleSTD='on',titleRMS='on',widthOBS=0.7,widthCOR=0.2,widthRMS=1,widthSTD=0.7,tickSTD=[1,2,3,0,-1,-2,-3],tickCOR=[1, .99, .95, .90, .8, .6, .3, -.3, -.6, -.8, -.9, -.95, -.99, -1])

# Plot Target Diagram.
fig.add_subplot(223)
sm.target_diagram(bias,crmsd,rmsd,circles=[1,2,3,4,5], markersize=9,ticks=np.arange(-5,6,1),equalAxes='on',circleLineWidth=0.7,markerLabelColor = 'b',markerColor = 'r', markerLegend = 'off')
fig.subplots_adjust(left=0, bottom=0, right=1, top=1,wspace=0, hspace=1)
fig.tight_layout()

# Save figure.
if choose_city == '1' and choose_var=='1':
    plt.savefig('taylor_target_paranagua_temp.png', transparent=False, bbox_inches ='tight', pad_inches=0, dpi=250)
if choose_city == '1' and choose_var=='2':
    plt.savefig('taylor_target_paranagua_ur.png', transparent=False, bbox_inches ='tight', pad_inches=0, dpi=250)
if choose_city == '1' and choose_var=='3':
    plt.savefig('taylor_target_paranagua_slp.png', transparent=False, bbox_inches ='tight', pad_inches=0, dpi=250)
if choose_city == '1' and choose_var=='4':
    plt.savefig('taylor_target_paranagua_wind.png', transparent=False, bbox_inches ='tight', pad_inches=0, dpi=250)

if choose_city == '2' and choose_var=='1':
    plt.savefig('taylor_target_indaial_temp.png', transparent=False, bbox_inches ='tight', pad_inches=0, dpi=250)
if choose_city == '2' and choose_var=='2':
    plt.savefig('taylor_target_indaial_ur.png', transparent=False, bbox_inches ='tight', pad_inches=0, dpi=250)
if choose_city == '2' and choose_var=='3':
    plt.savefig('taylor_target_indaial_slp.png', transparent=False, bbox_inches ='tight', pad_inches=0, dpi=250)
if choose_city == '2' and choose_var=='4':
    plt.savefig('taylor_target_indaial_wind.png', transparent=False, bbox_inches ='tight', pad_inches=0, dpi=250)

if choose_city == '3' and choose_var=='1':
    plt.savefig('taylor_target_floripa_temp.png', transparent=False, bbox_inches ='tight', pad_inches=0, dpi=250)
if choose_city == '3' and choose_var=='2':
    plt.savefig('taylor_target_floripa_ur.png', transparent=False, bbox_inches ='tight', pad_inches=0, dpi=250)
if choose_city == '3' and choose_var=='3':
    plt.savefig('taylor_target_floripa_slp.png', transparent=False, bbox_inches ='tight', pad_inches=0, dpi=250)
if choose_city == '3' and choose_var=='4':
    plt.savefig('taylor_target_floripa_wind.png', transparent=False, bbox_inches ='tight', pad_inches=0, dpi=250)

if choose_city == '4' and choose_var=='1':
    plt.savefig('taylor_target_saojoaquim_temp.png', transparent=False, bbox_inches ='tight', pad_inches=0, dpi=250)
if choose_city == '4' and choose_var=='2':
    plt.savefig('taylor_target_saojoaquim_ur.png', transparent=False, bbox_inches ='tight', pad_inches=0, dpi=250)
if choose_city == '4' and choose_var=='3':
    plt.savefig('taylor_target_saojoaquim_slp.png', transparent=False, bbox_inches ='tight', pad_inches=0, dpi=250)
if choose_city == '4' and choose_var=='4':
    plt.savefig('taylor_target_saojoaquim_wind.png', transparent=False, bbox_inches ='tight', pad_inches=0, dpi=250)
















