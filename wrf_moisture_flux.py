"""
File name:      wrf_moisture_flux.py
Author:         Ueslei Adriano Sutil
Email:          uesleisutil1@gmail.com
Created:        26 August 2019
Last modified:  16 October 2019
Version:        1.0
Python:         3.7.1

Calculate Vertical Integrated Moisture Flux Convergence from WRF output then
plot as contour.

Banacos, P. C. and Schultz, D. M. The Use of Moisture Flux Convergence in Forecasting 
Convective Initiation: Historical and Poeration Perspectives. Weather and Forecasting. 
v. 20, p. 351-366, 2005.

van Zomeren, J.; van Delden, A. Vertically integrated moisture flux convergence
as a predictor of thunderstorms. Atmospheric Research, v. 83, p. 435-445, 2007.
"""


# Python preable.
import numpy                as     np
from   wrf                  import getvar, to_np, interplevel, latlon_coords, extract_times
import netCDF4
import pandas               as     pd
import matplotlib.pyplot    as     plt
import matplotlib
from   mpl_toolkits.basemap import Basemap
import cmocean
from   progress.bar          import IncrementalBar
import os
matplotlib.use('Agg')
import warnings
warnings.filterwarnings("ignore")

# Set file directory and lon/lat.
wrf_file      = '/media/ueslei/Ueslei/INPE/PCI/Projetos/SC_2008/Outputs/normal/wrf.nc'
nc_wrf        = netCDF4.Dataset(wrf_file)
bbox          = [-53, -40, -32, -23]
clevs         = np.arange(-50,50.2,0.1)
ticks         = np.arange(min(clevs),max(clevs),20)  
cmap          = cmocean.cm.balance   
create_video  = False
ppt_fig       = False # If True, it will generate a figure with transparent background.

# 3. Start looping through time
atemp      = nc_wrf.variables['LH']
ntimes     = len(atemp)
range_loop = [i for i in range(168,ntimes,1)]
bar        = IncrementalBar('', max=len(range_loop))

for i in range_loop:
    timestr1   = extract_times(nc_wrf,timeidx=None,meta=False,do_xtime=False)
    timestr11  = timestr1[i]
    timestr    = pd.to_datetime(timestr11, format="%b %d %Y %H:%M")
    
    ua         = getvar(nc_wrf, "ua",timeidx=i,units="m s-1")
    va         = getvar(nc_wrf, "va",timeidx=i,units="m s-1")
    p          = getvar(nc_wrf,"pres",timeidx=i,units="hPa")
    qvapor     = nc_wrf.variables['QVAPOR'][i,:,:,:]
    qvapor     = qvapor/(1+qvapor)

    u1000 = interplevel(ua, p, 1000.0)
    lat_wrf, lon_wrf = latlon_coords(u1000, as_np=True)
    u1000 = u1000.values
    v1000 = interplevel(va, p, 1000.0)
    v1000 = v1000.values    
    q1000 = interplevel(qvapor, p, 1000.0)
    q1000 = q1000.values

    u925  = interplevel(ua, p, 925.0)
    u925  = u925.values
    v925  = interplevel(va, p, 925.0)
    v925  = v925.values
    q925  = interplevel(qvapor, p, 925.0)
    q925  = q925.values

    u850  = interplevel(ua, p, 850.0)
    u850  = u850.values
    v850  = interplevel(va, p, 850.0)
    v850  = v850.values
    q850  = interplevel(qvapor, p, 850.0)
    q850  = q850.values
    
    u700  = interplevel(ua, p, 700.0)
    u700  = u700.values
    v700  = interplevel(va, p, 700.0)
    v700  = v700.values
    q700  = interplevel(qvapor, p, 700.0)
    q700  = q700.values
      
    u500  = interplevel(ua, p, 500.0)
    u500  = u500.values
    v500  = interplevel(va, p, 500.0)
    v500  = v500.values
    q500  = interplevel(qvapor, p, 500.0)
    q500  = q500.values  

    q1    = (q1000+q925)/2
    u1    = (u1000+u925)/2
    v1    = (v1000+v925)/2
    dp1   = 1000-925   

    q2    = (q925+q850)/2
    u2    = (u925+u850)/2
    v2    = (v925+v850)/2
    dp2   = 925-850 

    q3    = (q850+q700)/2
    u3    = (u850+u700)/2
    v3    = (v850+v700)/2
    dp3   = 850-700

    q4    = (q700+q500)/2
    u4    = (u700+u500)/2
    v4    = (v700+v500)/2
    dp4   = 700-500

    qu1   = ((q1*u1)*dp1)
    qu2   = ((q2*u2)*dp2)
    qu3   = ((q3*u3)*dp3)
    qu4   = ((q4*u4)*dp4)

    qv1   = ((q1*v1)*dp1)
    qv2   = ((q2*v2)*dp2)
    qv3   = ((q3*v3)*dp3)
    qv4   = ((q4*v4)*dp4)

    qu1   = (1/9.81)*qu1
    qu2   = (1/9.81)*qu2
    qu3   = (1/9.81)*qu3
    qu4   = (1/9.81)*qu4

    qv1   = (1/9.81)*qv1
    qv2   = (1/9.81)*qv2
    qv3   = (1/9.81)*qv3
    qv4   = (1/9.81)*qv4

    nan1 = np.isnan(qu1)
    qu1[nan1] = 0
    nan2 = np.isnan(qu2)
    qu2[nan2] = 0
    nan3 = np.isnan(qu3)
    qu3[nan3] = 0
    nan4 = np.isnan(qu4)
    qu4[nan4] = 0
   
    nan1 = np.isnan(qv1)
    qv1[nan1] = 0
    nan2 = np.isnan(qv2)
    qv2[nan2] = 0
    nan3 = np.isnan(qv3)
    qv3[nan3] = 0
    nan4 = np.isnan(qv4)
    qv4[nan4] = 0

    qu = qu1+qu1+qu3+qu4
    qv = qv1+qv1+qv3+qv4    

    vimfc1 = (q1000*u1000)+(q1000+u1000)*dp1
    vimfc2 = (q1000*u1000)+(q1000+u1000)*dp1

    qu = (((u925*q925+u1000*q1000)/2)*(1000 - 925))+(((u850*q850+u925*q925)/2)*(925-850))+(((u700*q700+u850*q850)/2)*(850-700))+(((u500*q500+u700*q700)/2)*(700-500))

    qv = (((v925*q925+v1000*q1000)/2)*(1000 - 925))+(((v850*q850+v925*q925)/2)*(925-850))+(((v700*q700+v850*q850)/2)*(850-700))+(((v500*q500+v700*q700)/2)*(700-500))

    m   = Basemap(projection='merc',llcrnrlat=bbox[2],urcrnrlat=bbox[3],llcrnrlon=bbox[0],urcrnrlon=bbox[1], lat_ts=30,resolution='i')
    fig = plt.figure(1,figsize=(10,8)) 
    plt.xlabel('Longitude'u' [\N{DEGREE SIGN}]',labelpad=18,size=10)
    plt.ylabel('Latitude'u' [\N{DEGREE SIGN}]',labelpad=33,size=10)
    ax  = fig.add_subplot(111)
    plt.title(timestr, fontsize=11)
    m.drawparallels(np.arange(-90.,120.,1.), linewidth=0.00, color='black', labels=[1,0,0,1],labelstyle="N/S",fontsize=10)
    m.drawmeridians(np.arange(-180.,180.,2.), linewidth=0.00,color='black', labels=[1,0,0,1],labelstyle="N/S",fontsize=10)
    m.fillcontinents(color = '#ffffff')
    m.drawcountries(color = '#000000',linewidth=0.5)
    m.drawcoastlines(color = '#000000',linewidth=0.5)
    m.drawstates(color = '#000000',linewidth=0.5)       
 
    nsub  = 7
    scale = 0.02
    lon_wrf, lat_wrf = m(to_np(lon_wrf), to_np(lat_wrf))
    C = ax.quiver(lon_wrf[::nsub,::nsub],lat_wrf[::nsub,::nsub],qu[::nsub,::nsub],qv[::nsub,::nsub],alpha=1,scale=1/scale, zorder=1e35, width=0.0025,color='black',pivot='middle')
    wk = ax.quiverkey(C, .75, -.1, 2, ' Moisture Integrated Flux Vector\n 2 Kg.m⁻¹s⁻¹', coordinates='axes',color='#444444',labelsep=0.05, labelcolor='black',fontproperties={'size': '9'})
       
    try:
        os.makedirs("wrf_vimfc")
    except FileExistsError:
        pass
    if ppt_fig==True:
        plt.savefig('./wrf_vimfc/vimfc_{0:03d}.png'.format(i), transparent=True, bbox_inches = 'tight', pad_inches=0, dpi=250)
    if ppt_fig==False:
        plt.savefig('./wrf_vimfc/vimfc_{0:03d}.png'.format(i), transparent=False, bbox_inches = 'tight', pad_inches=0, dpi=250)        
    plt.clf()
    bar.next()
bar.finish()
# 4. Create mp4 file from figures.
if create_video==True:
    cwd = os.getcwd()
    exists = os.path.isfile('./wrf_vimfc/wrf_vimfc.mp4')
    if exists==True:
        os.system("rm -rf ./wrf_vimfc/wrf_vimfc.mp4")
        os.system("ffmpeg -r 10 -pattern_type glob -i '"+cwd+"/wrf_vimfc/*.png' -c:v libx264 -pix_fmt yuv420p -crf 15 ./wrf_vimfc/wrf_vimfc.mp4")
    else:
        os.system("ffmpeg -r 10 -pattern_type glob -i '"+cwd+"/wrf_vimfc/*.png' -c:v libx264 -pix_fmt yuv420p -crf 15 ./wrf_vimfc/wrf_vimfc.mp4")
        print(bg.da_cyan+'Wish to delete the .png files? (1) Yes or (2) No.'+bg.rs)
        removefiles = input()
        if removefiles=='1':
            os.system("rm -rf ./wrf_vimfc/*.png")
        if removefiles=='2':
            pass


