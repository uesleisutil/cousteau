#!/home/uesleisutil/anaconda3/bin/python python
# -*- coding: utf-8 -*-

"""
File name:      wrf_vimfc.py
Author:         Ueslei Adriano Sutil
Email:          uesleisutil1@gmail.com
Created:        26 August 2019
Last modified:  26 August 2019
Version:        1.0
Python:         3.7.1

Calculate Vertical Integrated Moisture Flux Convergence from WRF output then
plot as contour.

RAO, V. B.; CHAPA, S.R.; FRAQUITO, S. H.; Decadal variation of atmospheric 
ocean interaction in the Tropical Atlantic Ocean and its relationship to 
the northwest-brazil Rainfall. Journal of Meteorological Society of Japan,
v. 77, n.1, p. 63-77, 1999.
"""

from wrf import to_np, getvar, latlon_coords, smooth2d, extract_times

wrf_file = '/media/ueslei/Ueslei/INPE/PCI/Projetos/SC_2008/Outputs/cold_sst_100/wrf.nc'
bbox      = [-52.5, -45.5, -30.5, -25.5]
uvmet10 = getvar(nc_wrf, "uvmet10", i, units="m s-1")