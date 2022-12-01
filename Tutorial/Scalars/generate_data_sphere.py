import numpy as np
from netCDF4 import Dataset

##
#
# The generation script builds a full-sphere data set
#   at 0.5 degree resolution. It consists of a large
#   collection of 'balls' of differing sizes and densities,
#   randomly distributed around the globe.
#
##

# Make grid
Nlon, Nlat = int(360*1), int(180*1)

Llon, Llat = 2 * np.pi, np.pi

dlon = Llon / Nlon
dlat = Llat / Nlat

lon = np.arange( dlon/2, Llon, dlon )       - Llon/2
lat = np.arange( dlat/2, Llat, dlat )[1:-1] - Llat/2  # remove the actual poles
Nlat -= 2

LON, LAT = np.meshgrid(lon, lat)

R_earth = 6371e3
dAreas = R_earth**2 * np.cos(LAT) * dlat * dlon

# Define some eddy parameters
ball_scales = [ 250e3, 750e3, 3500e3 ]
ball_counts = [ 1e3,   3e2,   5e1 ]
ball_magnis = [ 0.2,   0.4,   0.1 ]

# Function to compute distances on sphere
def dist( lon0, lat0, LON = LON, LAT = LAT, R = 6371e3 ):

    del_LON = np.abs( LON - lon0 )
    del_LAT = np.abs( LAT - lat0 )

    numer  = ( np.cos( LAT ) * np.sin( del_LON ) )**2
    numer += ( np.cos( lat0 ) * np.sin( LAT ) - np.sin( lat0 ) * np.cos( LAT ) * np.cos( del_LON) )**2
    numer  = np.sqrt(numer)

    denom = np.sin( lat0 ) * np.sin( LAT ) + np.cos( lat0 ) * np.cos( LAT ) * np.cos( del_LON )

    return R * np.arctan2( numer, denom )


# Build streamfunction
rho   = np.zeros( (Nlat,Nlon) )

for scale, count, magnitude in zip( ball_scales, ball_counts, ball_magnis ):

    for ii in range(int(count)):
        W = ( scale * (1 + 0.05 * np.random.randn() ) ) / 2

        lon0 = np.random.choice( lon )  # choose longitude at random
        lat0 = np.random.choice( lat, p = np.cos(lat) / np.sum( np.cos(lat) ) )   # choose latitude, weighting based on area

        D = dist( lon0, lat0 )

        sign = np.random.choice( (-1,1) )

        rho += sign * magnitude * np.exp( -(D/W)**2 )

# Save flow to a file
dtype_dim = np.float64
dtype = np.float32

dims = ('time','depth','latitude','longitude')

fill_value = -1e10

with Dataset('density_sample.nc', 'w', format='NETCDF4') as fp:

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

    #
    rho_var = fp.createVariable('rho', dtype, dims, contiguous=True, fill_value = fill_value)
    rho_var.scale_factor = 1.

    rho_var[0,0,:,:] = rho
