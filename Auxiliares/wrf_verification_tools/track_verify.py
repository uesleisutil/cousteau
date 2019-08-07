#####################################################################
#
# Tropical Cyclone Tracks (Python)
# Created by Jason W. Godwin (jgodwin@rsmas.miami.edu)
#
# Version: 1.0
# 
# Description: This program uses a CSV containing date and time in
# Unix Epoch format (column 1), latitude (column 2), longitude
# (column 3), min. central pressure (column 3), and max. sustained
# winds (column 4) to plot the tropical cyclone path and colour-coding
# each point according to the storm's Saffir-Simpson Scale category.
#
# Dependencies:
# Python: basemap, matplotlib, and numpy
# Files: cities.csv (contains the city locations that will plot on
# the map)
#
#####################################################################

from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.dates as mdates
import datetime as dt
import time
import sys
import math

# specify csv_file path

################ EDIT THIS BLOCK ONLY! ###########################################

best_track_file = 'irene2011.csv'	# CSV containing best track information
wrf_track_file = 'irene_wrf19km.csv'	# CSV containing WRF track information
tc_name = 'Hurricane Irene'		# TC name (include the title "Hurricane" or "Tropical Storm"
tc_year = '2011'			# TC year
outfile_name = 'irenewrf_verify.png'	# what do you want to call the output file
use_map = 'east_coast'			# what map to use (full, gulf, carib, east_coast)
grid ='off'				# lat/lon grid on or off?

###### Function for computing track errors ########################################################
radius_of_earth = 6378.1 # radius of the earth in km
def track_error(tclat,tclon,wrflat,wrflon):
        # convert lat and lon to spherical coordinate in radians
        degrees_to_radians = math.pi/180.0

        # phi = 90 - latitude
        tcphi = (90.0 - tclat) * degrees_to_radians
        wrfphi = (90.0 - wrflat) * degrees_to_radians

        # theta = longitude
        tctheta = tclon * degrees_to_radians
        wrftheta = wrflon * degrees_to_radians

        # compute spherical distance from spherical coordinates
        cos = (math.sin(tcphi)*math.sin(wrfphi)*math.cos(tctheta-wrftheta) +
                math.cos(tcphi)*math.cos(wrfphi))
        arc = math.acos(cos)
        track_error = arc * radius_of_earth
        return track_error

# open new figure
fig = plt.figure(figsize=(11,8))

####################### Draw background map ######################################################
# setup Lambert conformal basemap
if use_map == 'full':
	print('Using Full Atlantic Basin Map')
	m = Basemap(width=10000000,height=7000000,projection='lcc',resolution='c',lat_0=25,lon_0=-50.)
	citycsv = 'cities/cities.csv'
elif use_map == 'gulf':
	print('Using Gulf of Mexico Map')
	m = Basemap(width=2000000,height=1600000,projection='lcc',resolution='c',lat_0=25, lon_0=-90.)
	citycsv = 'cities/gulf_cities.csv'
	lat_max = 32.0
	lat_min = 17.0
	lon_max = -100.0
	lon_min = -79.0
elif use_map == 'carib':
	print('Using Caribbean Map')
	m = Basemap(width=4000000,height=2500000,projection='lcc',resolution='c',lat_0=15, lon_0=-75.)
	citycsv = 'cities/carib_cities.csv'
	lat_max = 24.0 
	lat_min = 8.0
	lon_max = -89.0
	lon_min = -59.0
elif use_map == 'east_coast':
	print('Using East Coast Map')
	m = Basemap(width=4000000,height=3500000,projection='lcc',resolution='c',lat_0=32, lon_0=-65.)
	citycsv = 'cities/eastus_cities.csv'
	lat_max = 45.0 
	lat_min = 20.0
	lon_max = -82.0
	lon_min = -52.0
else:
	sys.exit('Please use either full, gulf, carib, or east_coast for your map!')

# draw the land-sea mask
print('Drawing map...')
m.drawlsmask(land_color='coral',ocean_color='aqua',lakes='True')

# draw various boundaries
m.drawstates(color='grey')
m.drawcountries(color='white')

# draw and label lat and lon grid
if grid == 'on':
	parallels = np.arange(-80.,81,10.)
	meridians = np.arange(10.,351.,20.)
	m.drawparallels(parallels,labels=[False,True,False,False],color='white')
	m.drawmeridians(meridians,labels=[False,False,False,True],color='white')

# add shaded relief background
m.bluemarble()

###################### Add cities ################################################################
print('Adding cities to map...')
cities = np.recfromcsv(citycsv, unpack=True, names=['city', 'clat', 'clon'], dtype=None)

for i in range(len(cities.city)):
	lon, lat = cities.clon[i], cities.clat[i]
	xpt, ypt = m(lon, lat)
	lonpt, latpt = m(xpt, ypt, inverse=True)
	m.plot(xpt, ypt, 'w+')
	if lon > lon_max and lon < lon_min and lat > lat_min and lat < lat_max:
		plt.text(xpt+50000,ypt+50000,cities.city[i],fontsize=6,color='White')

################### Plot best-track ##############################################################
# unpack the storm and WRF CSV files
print('Plotting storm...')
tc = np.recfromcsv(best_track_file, unpack=True, names=['times', 'tclat', 'tclon',
                       'tcpres', 'tcwind'], dtype=None)
wrf = np.recfromcsv(wrf_track_file, unpack=True, names=['times','wrflat','wrflon',
        'wrfpres','wrfwind'],dtype=None)

# plot the best track
for j in range(len(tc.times)):
	for p in range(len(wrf.times)):
		if tc.times[j] == wrf.times[p]:
			lon, lat = tc.tclon[j], tc.tclat[j]
			xpt, ypt = m(lon, lat)
			lonpt, latpt = m(xpt, ypt, inverse=True)

			m.plot(xpt, ypt, 'ro')

			# plot date for 00 UTC positions
			if tc.times[j] % 1000 == 0:
				# convert epoch time to standard GMT time
				day = time.strftime('%m/%d', time.gmtime(tc.times[j]))
				if lon > lon_max and lon < lon_min and lat > lat_min and lat < lat_max:
					plt.text(xpt+50000,ypt,day,fontsize=6,color='Red')

################ Plot WRF track ###########################################################
track_time = []
error = []
for k in range(len(tc.times)):
	for n in range(len(wrf.times)):
		if wrf.times[n] == tc.times[k]:
			lon, lat = wrf.wrflon[n], wrf.wrflat[n]
			xpt, ypt = m(lon, lat)
			lonpt, latpt = m(xpt, ypt, inverse=True)

			m.plot(xpt, ypt, 'yo')

			if wrf.times[n] % 1000 == 0:
                		# convert epoch time to standard GMT time
                		day = time.strftime('%m/%d', time.gmtime(wrf.times[n]))
                		if lon > lon_max and lon < lon_min and lat > lat_min and lat < lat_max:
                        		plt.text(xpt+50000,ypt,day,fontsize=6,color='Yellow')

			track_time.append(tc.times[k])
			error.append(track_error(tc.tclat[k],tc.tclon[k],wrf.wrflat[n],wrf.wrflon[n]))

print('Creating image files...')

title = tc_name + ' (' + tc_year +')'
plt.title(title)
plt.savefig(outfile_name,orientation='landscape',bbox_inches='tight')

gtime = [dt.datetime.utcfromtimestamp(date) for date in track_time]

fig2 = plt.figure(figsize=(11,8))
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%m/%d/%Y %H:%M'))
plt.gca().xaxis.set_major_locator(mdates.DayLocator())
plt.plot(gtime, error)
plt.gcf().autofmt_xdate()
plt.xlabel('Date and Time (GMT)')
plt.ylabel('Track Error (km)')
plt.grid()
plt.title(tc_name + ' WRF Track Error')
plt.savefig('track_error.png')

print('Great success!!!')
