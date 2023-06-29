import numpy as np
from netCDF4 import Dataset

##
#
#   This script defines geographic regions of interest.
#       They will be used in the post-processing. That is,
#       area-averaging will be performed over each of the
#       provided geographic regions.
#
##

vel_file = 'velocity_sample.nc'

with Dataset(vel_file, 'r') as dset:
    latitude  = dset['latitude'][:]
    longitude = dset['longitude'][:]

Nlat = len(latitude)
Nlon = len(longitude)

LON, LAT = np.meshgrid( longitude, latitude )


# Set up the storage information
dtype = np.float32
dtype_dim  = np.float64
fill_value = -1e10

dims = ('region', 'latitude', 'longitude')

# Give our regions meaningful names. This will be stored in the output files as well
#   You can call them whatever you'd like, but try to avoid special characters, for the
#   sake of paranoia.
regions = ['Global', 'Tropics', 'North of Tropics', 'South of Tropics']

# Now create a file to run through coarse-graining
with Dataset('region_definitions.nc', 'w', format='NETCDF4') as fp:

    # region
    dim = 'region'
    r_dim = fp.createDimension(dim, len(regions))
    r_var = fp.createVariable(dim, np.str, (dim,))
    for II, reg in enumerate(regions):
        r_var[II] = reg

    # lat
    dim = 'latitude'
    lat_dim = fp.createDimension(dim, Nlat)
    lat_var = fp.createVariable(dim, dtype_dim, (dim,))
    lat_var[:] = latitude

    # lon
    dim = 'longitude'
    lon_dim = fp.createDimension(dim, Nlon)
    lon_var = fp.createVariable(dim, dtype_dim, (dim,))
    lon_var[:] = longitude

    #
    reg_var = fp.createVariable('region_definition', dtype, dims, contiguous=True, fill_value = fill_value)
    reg_var.scale_factor = 1.

    # These regions are simple geographic definitions
    reg_var[0,:,:] = 1.
    reg_var[1,:,:] = np.abs(LAT) <= 15
    reg_var[2,:,:] =        LAT  >  15
    reg_var[3,:,:] =        LAT  < -15

