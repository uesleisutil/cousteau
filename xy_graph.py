#%%
from   wrf                   import ll_to_xy,vinterp,getvar
from   netCDF4               import Dataset
import numpy                 as     np
from   matplotlib            import path 
from   tqdm                  import tqdm
import csv

def bbox2ij(lon,lat,bbox=[-160., -155., 18., 23.]):
    """Return indices for i,j that will completely cover the specified bounding box.

    i0,i1,j0,j1 = bbox2ij(lon,lat,bbox)

    lon,lat = 2D arrays that are the target of the subset
    bbox = list containing the bounding box: [lon_min, lon_max, lat_min, lat_max]

    Example
    -------  
    i0,i1,j0,j1 = bbox2ij(lon_rho,[-71, -63., 39., 46])
    h_subset = nc.variables['h'][j0:j1,i0:i1]       
    """
    bbox=np.array(bbox)
    mypath=np.array([bbox[[0,1,1,0]],bbox[[2,2,3,3]]]).T
    p = path.Path(mypath)
    points = np.vstack((lon.flatten(),lat.flatten())).T
    n,m = np.shape(lon)
    inside = p.contains_points(points).reshape((n,m))
    ii,jj = np.meshgrid(range(m),range(n))
    return min(ii[inside]),max(ii[inside]),min(jj[inside]),max(jj[inside])

def rot2d(x, y, ang):
    """Rotate vectors by geometric angle
        
    This routine is part of Rob Hetland's OCTANT package:
    https://github.com/hetland/octant
    """
    xr = x*np.cos(ang) - y*np.sin(ang)
    yr = x*np.sin(ang) + y*np.cos(ang)
    return xr, yr

def shrink(a,b):
    """Return array shrunk to fit a specified shape by triming or averaging.
    
    a = shrink(array, shape)
    
    array is an numpy ndarray, and shape is a tuple (e.g., from
    array.shape). a is the input array shrunk such that its maximum
    dimensions are given by shape. If shape has more dimensions than
    array, the last dimensions of shape are fit.
    
    as, bs = shrink(a, b)
    
    If the second argument is also an array, both a and b are shrunk to
    the dimensions of each other. The input arrays must have the same
    number of dimensions, and the resulting arrays will have the same
    shape.
    
    This routine is part of Rob Hetland's OCTANT package:
    https://github.com/hetland/octant
    
    Example
    -------
    
    >>> shrink(rand(10, 10), (5, 9, 18)).shape
    (9, 10)
    >>> map(shape, shrink(rand(10, 10, 10), rand(5, 9, 18)))        
    [(5, 9, 10), (5, 9, 10)]   
    """

    if isinstance(b, np.ndarray):
        a = shrink(a,b.shape)
        b = shrink(b,a.shape)
        return (a, b)

    if isinstance(b, int):
        b = (b,)

    if len(a.shape) == 1:                # 1D array is a special case
        dim = b[-1]
        while a.shape[0] > dim:          # only shrink a
            if (dim - a.shape[0]) >= 2:  # trim off edges evenly
                a = a[1:-1]
            else:                        # or average adjacent cells
                a = 0.5*(a[1:] + a[:-1])
    else:
        for dim_idx in range(-(len(a.shape)),0):
            dim = b[dim_idx]
            a = a.swapaxes(0,dim_idx)        # put working dim first
            while a.shape[0] > dim:          # only shrink a
                if (a.shape[0] - dim) >= 2:  # trim off edges evenly
                    a = a[1:-1,:]
                if (a.shape[0] - dim) == 1:  # or average adjacent cells
                    a = 0.5*(a[1:,:] + a[:-1,:])
            a = a.swapaxes(0,dim_idx)        # swap working dim back

    return a


for greatloop in range(0,10,1):
    if greatloop == 0:
        print('Dataset: 100_cold and 100_cold')
        display   = '100_cold and 100_cold'
        roms_file = '/media/ueslei/Ueslei_HD/INPE/PCI/Projetos/SC_2008/Outputs/100_cold/roms.nc'
        wrf_file  = '/media/ueslei/Ueslei_HD/INPE/PCI/Projetos/SC_2008/Outputs/100_cold/wrf.nc'
    if greatloop == 1:
        print('Dataset: 80_cold and 80_cold')
        display   = '80_cold and 80_cold'       
        roms_file = '/media/ueslei/Ueslei_HD/INPE/PCI/Projetos/SC_2008/Outputs/80_cold/roms.nc'
        wrf_file  = '/media/ueslei/Ueslei_HD/INPE/PCI/Projetos/SC_2008/Outputs/80_cold/wrf.nc'
    if greatloop == 2:
        print('Dataset: 60_cold and 60_cold')
        display   = '60_cold and 60_cold'
        roms_file = '/media/ueslei/Ueslei_HD/INPE/PCI/Projetos/SC_2008/Outputs/60_cold/roms.nc'
        wrf_file  = '/media/ueslei/Ueslei_HD/INPE/PCI/Projetos/SC_2008/Outputs/60_cold/wrf.nc'
    if greatloop == 4:
        print('Dataset: 40_cold and 40_col')
        display   = '40_cold and 40_cold'
        roms_file = '/media/ueslei/Ueslei_HD/INPE/PCI/Projetos/SC_2008/Outputs/40_cold/roms.nc'
        wrf_file  = '/media/ueslei/Ueslei_HD/INPE/PCI/Projetos/SC_2008/Outputs/40_cold/wrf.nc'
    if greatloop == 5:
        print('Dataset: 20_cold and 20_cold')
        display   = '20_cold and 20_cold'
        roms_file = '/media/ueslei/Ueslei_HD/INPE/PCI/Projetos/SC_2008/Outputs/40_cold/roms.nc' # ARRUMAR
        wrf_file  = '/media/ueslei/Ueslei_HD/INPE/PCI/Projetos/SC_2008/Outputs/40_cold/wrf.nc' # ARRUMAR
    if greatloop == 6:
        print('Dataset: 100_warm and 100_cold')
        display = '100_warm and 100_cold'
        roms_file = '/media/ueslei/Ueslei_HD/INPE/PCI/Projetos/SC_2008/Outputs/100_warm/roms.nc'
        wrf_file  = '/media/ueslei/Ueslei_HD/INPE/PCI/Projetos/SC_2008/Outputs/100_warm/wrf.nc'
    if greatloop == 7:
        print('Dataset: 80_warm and 50_cold')
        display = '80_warm and 80_cold'
        roms_file = '/media/ueslei/Ueslei_HD/INPE/PCI/Projetos/SC_2008/Outputs/80_warm/roms.nc'
        wrf_file  = '/media/ueslei/Ueslei_HD/INPE/PCI/Projetos/SC_2008/Outputs/80_warm/wrf.nc'
    if greatloop == 8:
        print('Dataset: 60_warm and 60_cold')
        display = '60_warm and 60_cold'
        roms_file = '/media/ueslei/Ueslei_HD/INPE/PCI/Projetos/SC_2008/Outputs/60_warm/roms.nc'
        wrf_file  = '/media/ueslei/Ueslei_HD/INPE/PCI/Projetos/SC_2008/Outputs/60_warm/wrf.nc'
    if greatloop == 9:
        print('Dataset: 40_warm and 40_cold')
        display = '40_warm and 40_cold'
        roms_file = '/media/ueslei/Ueslei_HD/INPE/PCI/Projetos/SC_2008/Outputs/40_warm/roms.nc'
        wrf_file  = '/media/ueslei/Ueslei_HD/INPE/PCI/Projetos/SC_2008/Outputs/40_warm/wrf.nc'
    if greatloop == 10:
        print('Dataset: 20_warm and 20_cold')
        display = '20_warm and 20_cold'
        roms_file = '/media/ueslei/Ueslei_HD/INPE/PCI/Projetos/SC_2008/Outputs/40_warm/roms.nc'
        wrf_file  = '/media/ueslei/Ueslei_HD/INPE/PCI/Projetos/SC_2008/Outputs/40_warm/wrf.nc'
     
    nc_roms   = Dataset(roms_file)
    nc_wrf    = Dataset(wrf_file)
    
    # Set map limits and ROMS level.
    bbox = [-53,-44,-31,-23]
    zlev = -1 # Last sigma layer corresponds to surface
    

    lon_rho = nc_roms.variables['lon_rho'][:]
    lat_rho = nc_roms.variables['lat_rho'][:]
    i0,i1,j0,j1 = bbox2ij(lon_rho,lat_rho,bbox)
    lon_var = lon_rho[j0:j1, i0:i1]
    lat_var = lat_rho[j0:j1, i0:i1]
    lon = lon_rho[(j0+1):(j1-1), (i0+1):(i1-1)]
    lat = lat_rho[(j0+1):(j1-1), (i0+1):(i1-1)]

    sst = nc_roms.variables['temp'][:, zlev,  j0:j1, i0:i1]
    sst[sst > 30000] = np.nan
    sst=np.nanmean(sst)
    
    blumenau = ll_to_xy(nc_wrf,-26.555,-49.355,timeidx=0, squeeze=False, meta=False, stagger=None, as_int=True)
    major = ll_to_xy(nc_wrf,-27.2457,-48.574,timeidx=0, squeeze=False, meta=False, stagger=None, as_int=True)
    
    # Count how many time-steps.
    ntimes = len(nc_wrf.variables['Times'][:])
    loop = [i for i in range(97,ntimes,1)]
    len_loop = len(loop)
    
    # Store variables.
    qa_blumenau = np.zeros([len_loop,1,1])
    wa_blumenau = np.zeros([len_loop,1,1])
    qa_major = np.zeros([len_loop,1,1])
    wa_major = np.zeros([len_loop,1,1])
    
    qa_blumenau_total = np.zeros([10])
    wa_blumenau_total = np.zeros([10])       
    qa_major_total = np.zeros([10])
    wa_major_total = np.zeros([10])      
                      
    for i in tqdm(range(len_loop), desc='2m Spec. Hum. for '+display):
        qa = getvar(nc_wrf, "QVAPOR", loop[i])
        qa = vinterp(nc_wrf,field=qa, vert_coord='ght_msl',interp_levels=[0.002],meta=False,extrapolate=True)
        qa = qa[0,blumenau[0],blumenau[1]]
        qa_blumenau[i,:,:] = qa/(1+qa)*1000 #Calculates Specific Humidity then convert kg/kg to g/kg
       
        qa = getvar(nc_wrf, "QVAPOR", loop[i])
        qa = vinterp(nc_wrf,field=qa, vert_coord='ght_msl',interp_levels=[0.002],meta=False,extrapolate=True)
        qa = qa[0,major[0],major[1]]
        qa_major[i,:,:] = qa/(1+qa)*1000 #Calculates Specific Humidity then convert kg/kg to g/kg
    
        qa_blumenau_2 = np.nanmean(qa_blumenau)
        qa_major_2 = np.nanmean(qa_major)
        qa_blumenau_total[greatloop] = qa_blumenau_2
        qa_major_total[greatloop] = qa_major_2

    for i in tqdm(range(len_loop), desc='10m Vert. Vel. for'+ display):
        wa = getvar(nc_wrf, "wa", timeidx=loop[i],units="m s-1")
        wa = vinterp(nc_wrf,field=wa, vert_coord='ght_msl',interp_levels=[0.01],meta=False,extrapolate=True)
        wa_blumenau = wa[0,blumenau[0],blumenau[1]]
       
        wa = getvar(nc_wrf, "wa", timeidx=loop[i],units="m s-1")
        wa = vinterp(nc_wrf,field=wa, vert_coord='ght_msl',interp_levels=[0.01],meta=False,extrapolate=True)
        wa_major = wa[0,major[0],major[1]]
    
        wa_blumenau_2 = np.nanmean(wa_blumenau)
        wa_major_2 = np.nanmean(wa_major)
        wa_blumenau_total[greatloop] = wa_blumenau_2
        wa_major_total[greatloop] = wa_major_2

try:
    fobj = open('blumenau_wa.csv', 'w')
    fields = ['cold', 'warm']

    blumenau = [
        {'cold_100': wa_blumenau_total[0], 'warm_100': wa_blumenau_total[5]},
        {'cold_80':  wa_blumenau_total[1], 'warm_80': wa_blumenau_total[6]},
        {'cold_60':  wa_blumenau_total[2], 'warm_60': wa_blumenau_total[7]},
        {'cold_40':  wa_blumenau_total[3], 'warm_40': wa_blumenau_total[8]},
        {'cold_20':  wa_blumenau_total[4], 'warm_20': wa_blumenau_total[9]}           
    ]
    writer = csv.DictWriter(fobj, blumenau)
    writer.writeheader()
    for row in blumenau:
        writer.writerow(row)
except:
    print("An error occurred while writing the file.")
finally:
    fobj.close()

try:
    fobj = open('blumenau_qa.csv', 'w')
    fields = ['cold', 'warm']

    blumenau = [
        {'cold_100': qa_blumenau_total[0], 'warm_100':qa_blumenau_total[5]},
        {'cold_80':  qa_blumenau_total[1], 'warm_80': qa_blumenau_total[6]},
        {'cold_60':  qa_blumenau_total[2], 'warm_60': qa_blumenau_total[7]},
        {'cold_40':  qa_blumenau_total[3], 'warm_40': qa_blumenau_total[8]},
        {'cold_20':  qa_blumenau_total[4], 'warm_20': qa_blumenau_total[9]}           
    ]
    writer = csv.DictWriter(fobj, blumenau)
    writer.writeheader()
    for row in blumenau:
        writer.writerow(row)
except:
    print("An error occurred while writing the file.")
finally:
    fobj.close()

try:
    fobj = open('major_qa.csv', 'w')
    fields = ['cold', 'warm']

    major = [
        {'cold_100': qa_major_total[0], 'warm_100':qa_major_total[5]},
        {'cold_80':  qa_major_total[1], 'warm_80': qa_major_total[6]},
        {'cold_60':  qa_major_total[2], 'warm_60': qa_major_total[7]},
        {'cold_40':  qa_major_total[3], 'warm_40': qa_major_total[8]},
        {'cold_20':  qa_major_total[4], 'warm_20': qa_major_total[9]}           
    ]
    writer = csv.DictWriter(fobj, major)
    writer.writeheader()
    for row in major:
        writer.writerow(row)
except:
    print("An error occurred while writing the file.")
finally:
    fobj.close()

try:
    fobj = open('major_wa.csv', 'w')
    fields = ['cold', 'warm']

    major = [
        {'cold_100': wa_major_total[0], 'warm_100':wa_major_total[5]},
        {'cold_80':  wa_major_total[1], 'warm_80': wa_major_total[6]},
        {'cold_60':  wa_major_total[2], 'warm_60': wa_major_total[7]},
        {'cold_40':  wa_major_total[3], 'warm_40': wa_major_total[8]},
        {'cold_20':  wa_major_total[4], 'warm_20': wa_major_total[9]}           
    ]
    writer = csv.DictWriter(fobj, major)
    writer.writeheader()
    for row in major:
        writer.writerow(row)
except:
    print("An error occurred while writing the file.")
finally:
    fobj.close()
