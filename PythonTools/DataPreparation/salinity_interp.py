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

    dtype = source.variables['so'].datatype
    
    # Get up the initial data
    so = source.variables['so'][:]

# Use a spline that only interpolates in time
if len(times_so) > 4:
    kind = 'cubic'  # cubic spline
else:
    kind = 'slinear' # linear spline

times_so = np.array(times_so)
interpolator = spi.interp1d(times_so, so, axis=0, kind=kind)
print(times_so[0], times_so[-1])

# Now create a file to run through coarse-graining
with Dataset(outfile, 'a', format='NETCDF4') as fp:
     
    time      = fp.variables['time'     ][:] 
    depth     = fp.variables['depth'    ][:]
    latitude  = fp.variables['latitude' ][:]
    longitude = fp.variables['longitude'][:]
    mask      = fp.variables['uo'][:].mask

    time = np.array(time)
    print(time[0], time[-1])
    s_interp = interpolator(time)
    s_interp = np.ma.masked_where(mask, s_interp)

    out_s = fp.createVariable('so', dtype, 
                        ('time','depth','latitude','longitude'),
                            contiguous = True, fill_value = -32767)

    out_s.scale_factor   = 1.

    out_s[:] = s_interp
    
