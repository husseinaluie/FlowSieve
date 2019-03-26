from netCDF4 import Dataset
import numpy as np
import scipy.interpolate as spi

# The purpose of this script is to interpolate
#  salinity (which is montly) onto the vel / temp
#  points (which is hourly)

source  = 'salinity.nc'
outfile = 'source.nc'

D2R = np.pi / 180.

rho0 = 1e3
g    = 9.81

# Get the velocity and SSH data
with Dataset(source, 'r') as source:

    # Get the dimensions and sizes
    times_so     = source.variables['time'][:]
    depths_so    = source.variables['depth'][:]
    latitude_so  = source.variables['latitude'][:]
    longitude_so = source.variables['longitude'][:]
    
    # Get up the initial data
    so = source.variables['so'][:]

#T_so, D_so, LAT_so, LON_so = np.meshgrid(times_so, depths_so, latitude_so, longitude_so)
#num_pts_sal = len(times_so) * len(depths_so) * len(latitude_so) * len(longitude_so)
#sal_pts = np.zeros((num_pts_sal, 4))
#sal_pts[:,0] = T_so.ravel()
#sal_pts[:,1] = D_so.ravel()
#sal_pts[:,2] = LAT_so.ravel()
#sal_pts[:,3] = LON_so.ravel()

# Use a linear spline that only interpolates in time
interpolator = spi.interp1d(times_so, so, axis=0, kind='slinear')

print(so.shape)

# Now create a file to run through coarse-graining
with Dataset(outfile, 'a', format='NETCDF4') as fp:
     
    time      = fp.variables['time'     ][:] 
    depth     = fp.variables['depth'    ][:]
    latitude  = fp.variables['latitude' ][:]
    longitude = fp.variables['longitude'][:]
    mask      = fp.variables['uo'][:].mask

    print(times_so, time)
    print(depths_so, depth)
    print(latitude_so.min(), latitude_so.max(), latitude.min(), latitude.max())
    print(longitude_so.min(), longitude_so.max(), longitude.min(), longitude.max())

    #T, D, LAT, LON = np.meshgrid(time, depth, latitude, longitude)
    #
    #num_pts = len(time) * len(depth) * len(latitude) * len(longitude)
    #pts = np.zeros((num_pts, 4))
    #pts[:,0] = T.ravel()
    #pts[:,1] = D.ravel()
    #pts[:,2] = LAT.ravel()
    #pts[:,3] = LON.ravel()
    #
    #s_interp = spi.griddata(sal_pts, so.ravel(), pts, method='linear')

    s_interp = interpolator(time)
    print(s_interp.shape)

    s_interp = np.ma.masked_where(mask, s_interp)

    out_s = fp.createVariable('so', np.float64, 
                        ('time','depth','latitude','longitude'),
                            contiguous = True, fill_value = -32767)

    out_s.scale_factor   = 1.

    out_s[:] = s_interp
    
