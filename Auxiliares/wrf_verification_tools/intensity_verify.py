# Jason's WRF Verification Python Program
# Version: 1.0
# Last modified:	12 July 2013
#
# Created by:
# Jason W. Godwin
# University of Miami
# Rosenstiel School of Marine and Atmospheric Sciences
# Division of Meteorology and Physical Oceanography
#
# Version history
#
# 1.0 - Original build
#
# INSTRUCTIONS TO USER
# 
# To use this program, create a comma-separated values (CSV) file with date and time
# in Unix Epoch format (GMT) in the first column, best track pressure in the second
# column, WRF pressure in the third column, best track winds in the fourth, and WRF
# winds in the fifth. You should only have to modify the "edit this block only" block.
#
#################################################################################

# Edit this block only ##########################################################

tc_csv = 'irene2011.csv'	# MODIFY THIS LINE TO CHANGE BEST-TRACK CSV FILENAME
wrf_csv = 'irene_wrf19km.csv'	# MODIFY THIS LINE TO CHANGE WRF CSV FILENAME (INCLUDE .CSV SUFFIX)
storm_name = 'Hurricane Irene'  # MODIFY THIS LINE TO CHANGE STORM NAME (INCLUDE HURRICANE/T.S./T.D.)
storm_year = '2011'             # MODIFY THIS LINE TO CHANGE STORM YEAR
case_name = 'irene19km'		# MODIFY THIS LINE TO CHANGE THE CASE NAME (THIS WILL BE PART OF FILENAME)

#################################################################################

# Create storm title for plots
storm_title = storm_name + ' (' + storm_year + ')' + ' WRF Verification'

# Import the necessary modules
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import datetime as dt

# Import data from CSV: Edit line 70, first argument in np.recfromcsv to change
# the target CSV file to import

print('Reading CSV')
tc = np.recfromcsv(tc_csv, unpack=True, names=['bt_times','bt_lat','bt_lon','bt_pres',
	'bt_wind'],dtype=None)
wrf = np.recfromcsv(wrf_csv, unpack=True, names=['wrf_times','wrf_lat','wrf_lon',
                       'wrf_pres', 'wrf_wind'], dtype=None)

# Calculate WRF error

pres_error = []
wind_error = []
wrf_time = []
bt_pres = []
wrf_pres = []
bt_wind = []
wrf_wind = []

print('Computing WRF Errors')
for j in range(len(tc.bt_times)):
	for i in range(len(wrf.wrf_times)):
		if tc.bt_times[j] == wrf.wrf_times[i]:
			pres_error.append(wrf.wrf_pres[i] - tc.bt_pres[j])
			wind_error.append(wrf.wrf_wind[i] - tc.bt_wind[j])
			wrf_time.append(tc.bt_times[j])
			bt_pres.append(tc.bt_pres[j])
			wrf_pres.append(wrf.wrf_pres[i])
			bt_wind.append(tc.bt_wind[j])
			wrf_wind.append(wrf.wrf_wind[i])
print(wind_error)

gtime = [dt.datetime.utcfromtimestamp(date) for date in wrf_time]

print('Performing Pressure Verification')
# Pressure plot
fig = plt.figure(1)
plt.subplot(211)
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%m/%d/%Y %H:%M'))
plt.gca().xaxis.set_major_locator(mdates.DayLocator())

# Line plots for best track and WRF pressure
plt.plot(gtime,bt_pres,label='Best Track Pressure')
plt.plot(gtime,wrf_pres,label='WRF Forecast Pressure')
plt.gcf().autofmt_xdate()

plt.xlabel('Date and Time (GMT)',fontsize=8)
plt.ylabel('Minimum Central Pressure (mb)',fontsize=8)
plt.grid()
plt.title(storm_title + ' - Min. Pressure')
plt.legend(loc=1,borderaxespad=0.25,prop={'size':8})

ax = plt.subplot(212)
# Bar plot for WRF forecast error
ax.bar(gtime,pres_error,width=0.20,label='WRF Forecast Error')
plt.ylim([-1*max(pres_error)+1,max(pres_error)+1])
plt.legend(loc=1,borderaxespad=0.25,prop={'size':8})

plt.xlabel('Date and Time (GMT)',fontsize=8)
plt.ylabel('WRF Forecast Pressure Error',fontsize=8)
ax.xaxis_date()
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%m/%d/%Y %H:%M'))
plt.gca().xaxis.set_major_locator(mdates.DayLocator())
plt.gcf().autofmt_xdate()
plt.grid()

filename1 = case_name + '_pressure.png'
plt.savefig(filename1)

#---------------------------------------------------------------------------------
print('Performing Wind Verification')
# Wind plot
fig = plt.figure(2)
plt.subplot(211)
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%m/%d/%Y %H:%M'))
plt.gca().xaxis.set_major_locator(mdates.DayLocator())

# Line plots for best track and WRF pressure
plt.plot(gtime,bt_wind,label='Best Track Wind')
plt.plot(gtime,wrf_wind,label='WRF Forecast Wind')
plt.gcf().autofmt_xdate()

plt.xlabel('Date and Time (GMT)',fontsize=8)
plt.ylabel('Maximum Sustained Winds (knots)',fontsize=8)
plt.grid()
plt.title(storm_title + ' - Max Winds')
plt.legend(loc=1,borderaxespad=0.25,prop={'size':8})

ax = plt.subplot(212)
# Bar plot for WRF forecast error
ax.bar(gtime,wind_error,width=0.20,label='WRF Forecast Error')

# Get y limits
derp1 = min(wind_error)
derp2 = max(wind_error)

derp1 = abs(derp1)
derp2 = abs(derp2)

if derp1 > derp2:
	maxy = derp1 - 1
	miny = -1 * derp1 + 1
else:
	maxy = derp2 - 1
	miny = -1 * derp2 + 1

plt.ylim([miny,maxy])
plt.legend(loc=1,borderaxespad=0.25,prop={'size':8})

plt.xlabel('Date and Time (GMT)',fontsize=8)
plt.ylabel('WRF Forecast Wind Error',fontsize=8)
ax.xaxis_date()
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%m/%d/%Y %H:%M'))
plt.gca().xaxis.set_major_locator(mdates.DayLocator())
plt.gcf().autofmt_xdate()
plt.grid()

filename2 = case_name + '_wind.png'

plt.savefig(filename2)
print('Great Success!!!')
