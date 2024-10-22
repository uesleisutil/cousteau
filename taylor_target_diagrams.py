#!/home/uesleisutil/anaconda3/bin/python python
# -*- coding: utf-8 -*-

"""
File name:      taylor_target_diagrams.py
Author:         Ueslei Adriano Sutil
Email:          uesleisutil1@gmail.com
Created:        16 August 2019
Last modified:  20 August 2019
Version:        2.1
Python:         3.7.1

Extract data from ASCII files then create Taylor and Target Diagrams.
"""

import matplotlib.pyplot as plt
import skill_metrics as sm
import pandas as pd
import numpy as np

make_taylor = False
make_target = True
print('Choose plots: (1) Taylor and Target, (2) Taylor diagram or (3) Target diagram:')
choose_plot = input()

print('Choose data: (1) Paranaguá, (2) Indaial, (3) Florianópolis or (4) São Joaquim')
choose_city = input()
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
    var = 'Temperature'
    t_insitu     = data['T_EM']
    t_insitu     = t_insitu.values 
    t_model      = data['T_WRF']
    t_model      = t_model.values
    taylor_stats = sm.taylor_statistics(t_model,t_insitu)
    target_stats = sm.target_statistics(t_model,t_insitu)
if choose_var == '2':
    var          = 'Relative Humidity'
    ur_insitu    = data['UR_EM']
    ur_insitu    = ur_insitu.values
    ur_model     = data['UR_WRF']
    ur_model     = ur_model.values
    taylor_stats = sm.taylor_statistics(ur_model,ur_insitu)
    target_stats = sm.target_statistics(ur_model,ur_insitu)
if choose_var == '3':
    var          = 'Sea Level Pressure'
    slp_insitu   = data['SLP_EM']
    slp_insitu   = slp_insitu.values
    slp_model    = data['SLP_WRF']
    slp_model    = slp_model.values
    taylor_stats = sm.taylor_statistics(slp_model,slp_insitu)
    target_stats = sm.target_statistics(slp_model,slp_insitu)  
if choose_var == '4':
    var          = 'Wind Velocity'   
    vvel_insitu  = data['VVel_EM']
    vvel_insitu  = vvel_insitu.values
    vvel_model   = data['VVel_WRF']
    vvel_model   = vvel_model.values
    taylor_stats = sm.taylor_statistics(vvel_model,vvel_insitu)
    target_stats = sm.target_statistics(vvel_model,vvel_insitu)  

sdev  = np.array([taylor_stats['sdev'][0],taylor_stats['sdev'][1]])
crmsd = np.array([taylor_stats['crmsd'][1],taylor_stats['crmsd'][1]])
ccoef = np.array([taylor_stats['ccoef'][1],taylor_stats['ccoef'][1]])

bias  = np.array([target_stats['bias'], target_stats['bias']])
crmsd = np.array([target_stats['crmsd'], target_stats['crmsd']])
rmsd  = np.array([target_stats['rmsd'], target_stats['rmsd']])

if choose_plot == '1':
    fig = plt.figure(1,figsize=(15, 11))
    fig.set_figwidth(15)
    fig.add_subplot(221)
    plt.title(title +': ' + var ++ latlon, pad=60)
    sm.taylor_diagram(sdev,crmsd,ccoef, numberpanels=2, alpha = 0.0, checkstats = 'on', markerSize=24,markerColor='r',markerLegend='off', colOBS = 'crimson', colRMS='C0', colCOR='k',markerobs = '*', styleOBS='-',styleCOR='-', styleRMS=':', titleCOR='on',titleSTD='on',titleRMS='on',widthOBS=0.7,widthCOR=0.2,widthRMS=1,widthSTD=0.7,tickSTD=[1,2,3,0,-1,-2,-3],tickCOR=[1, .99, .95, .90, .8, .6, .3, -.3, -.6, -.8, -.9, -.95, -.99, -1])
    
    fig.add_subplot(223)
    sm.target_diagram(bias,crmsd,rmsd,circles=[1,2,3,4,5], markersize=9,ticks=np.arange(-5,6,1),equalAxes='on',circleLineWidth=0.7,markerLabelColor = 'b',markerColor = 'r', markerLegend = 'off')
    fig.subplots_adjust(left=0, bottom=0, right=1, top=1,wspace=0, hspace=1)
    fig.tight_layout()
    
    plt.savefig('TaylorTarget_'+title+'_'+var+'.png', transparent=False, bbox_inches ='tight', pad_inches=0, dpi=250)

if choose_plot == '2':
    plt.title(title+': '+  latlon + '\n' + var, pad=20, fontsize=7)
    #sm.taylor_diagram(sdev,crmsd,ccoef,numberpanels=2, alpha = 1.0, checkstats = 'on', markerSize=12,markerColor='r',markerLegend='off', colOBS = 'olive', colRMS='C0', colCOR='k',markerobs = '.', styleOBS='-',styleCOR='-', styleRMS=':', titleCOR='on',titleSTD='on',titleRMS='on',widthOBS=0.7,widthCOR=0.2,widthRMS=1,widthSTD=0.7,tickSTD=[2,4,6,8,10,0,-2,-4,-6,-8,-10],tickRMS=[2.5,5,7.5,10,12.5],tickCOR=[1, .99, .95, .90, .8, .6, .3, -.3, -.6, -.8, -.9, -.95, -.99, -1])
    sm.taylor_diagram(sdev,crmsd,ccoef, titlermsdangle=155, markerSize=12, markerColor='r',
                      markerLegend='off', colOBS = 'C0', colRMS='C0', 
                      colCOR='k', markerobs= '.', styleOBS='', styleSTD='--', styleCOR='-', 
                      styleRMS='-', titleCOR='on', titleSTD='on', titleRMS='on', widthOBS=0.7, 
                      widthCOR=0.4, widthRMS=0.3, widthSTD=0.4, tickSTD=[0,0.5,1,1.5,2],
                      tickRMS=[0.5,1,1.5,2], tickCOR=[ 0.99, 0.95, 0.9, 0.8, 0.6, 0.3 ])   
    plt.savefig('Taylor_'+title+'_'+var+'.png', transparent=False, bbox_inches ='tight', pad_inches=0, dpi=250)

if choose_plot == '3':
    plt.title(title +': ' + var +'\n'+ latlon, pad=40)
    sm.target_diagram(bias,crmsd,rmsd,circles=[1,2,3,4,5], markersize=9,ticks=np.arange(-5,6,1),equalAxes='on',circleLineWidth=0.7,markerLabelColor = 'b',markerColor = 'r', markerLegend = 'off')
    
    plt.savefig('Target_'+title+'_'+var+'.png', transparent=False, bbox_inches ='tight', pad_inches=0, dpi=250)
