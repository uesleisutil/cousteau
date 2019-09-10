#!/home/uesleisutil/anaconda3/bin/python python
# -*- coding: utf-8 -*-

"""
File name:      xy_sst_plot.py
Author:         Ueslei Adriano Sutil
Email:          uesleisutil1@gmail.com
Created:        10 July 2019
Last modified:  11 July 2019
Version:        1.0
Python:         3.7.1

Creates XY graph from TRMM (Huffman et al., 2007), in situ data and WRF output.

TRMM:
Huffman, G.J. (1997), Estimates of Root-Mean-Square Random Error for Finite Samples of Estimated Precipitation, J. Appl. Meteor.
Huffman, G.J., R.F. Adler, D.T. Bolvin, G. Gu, E.J. Nelkin, K.P. Bowman, Y. Hong, E.F. Stocker, D.B. Wolff (2007), The TRMM Multi-satellite Precipitation Analysis: Quasi- Global, Multi-Year, Combined-Sensor Precipitation Estimates at Fine Scale., J. Hydrometeor
Huffman, G.J., R.F. Adler, D.T. Bolvin, E.J. Nelkin (2010), The TRMM Multi-satellite Precipitation Analysis (TMPA). Chapter 1 in Satellite Rainfall Applications for Surface Hydrology, doi:10.1007/978-90-481-2915-7
"""

import netCDF4
import numpy as np
lonbounds = [-48.56,-48.56] 
latbounds = [-27.58,-27.58]

trmm_dir  = '/media/ueslei/Ueslei/INPE/PCI/Projetos/SC_2008/Dados/TRMM/trmm.nc'
trmm_file = netCDF4.Dataset(trmm_dir)
trmm_prec = trmm_file.variables['precipitation']
trmm_lat  = trmm_file.variables['lat']
trmm_lon  = trmm_file.variables['lon']


# Florianópolis 83897 Lat: -27.58 Lon:-48.56
# Indaial 83872 Lat: -26.9 Lon: -48.21
# Paranaguá 83844 Lat: -25.53 Lon: -48.51
