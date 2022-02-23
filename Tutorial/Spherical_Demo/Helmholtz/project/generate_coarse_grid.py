import numpy as np
from netCDF4 import Dataset

##
#
# The generation script builds a full-sphere data set
#   at 0.5 degree resolution. It consists of a large
#   collection of eddies of differing sizes, randomly
#   distributed around the globe.
#
##

# Make grid
Nlon, Nlat = int(360*2/3), int(180*2/3)

Llon, Llat = 2 * np.pi, np.pi

dlon = Llon / Nlon
dlat = Llat / Nlat

lon = np.arange( dlon/2, Llon, dlon )       - Llon/2
lat = np.arange( dlat/2, Llat, dlat )[1:-1] - Llat/2  # remove the actual poles
Nlat -= 2

# Save flow to a file
dtype_dim = np.float64

dims = ('time','depth','latitude','longitude')

fill_value = -1e10

with Dataset('coarse_grid.nc', 'w', format='NETCDF4') as fp:

    # time
    dim = 'time'
    t_dim = fp.createDimension(dim, 1)
    t_var = fp.createVariable(dim, dtype_dim, (dim,))
    t_var[:] = 0

    # depth
    dim = 'depth'
    d_dim = fp.createDimension(dim, 1)
    d_var = fp.createVariable(dim, dtype_dim, (dim,))
    d_var[:] = 0

    # lat
    dim = 'latitude'
    lat_dim = fp.createDimension(dim, Nlat)
    lat_var = fp.createVariable(dim, dtype_dim, (dim,))
    lat_var[:] = lat * 180. / np.pi

    # lon
    dim = 'longitude'
    lon_dim = fp.createDimension(dim, Nlon)
    lon_var = fp.createVariable(dim, dtype_dim, (dim,))
    lon_var[:] = lon * 180 / np.pi
