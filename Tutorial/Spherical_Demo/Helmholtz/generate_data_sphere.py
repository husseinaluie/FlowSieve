import numpy as np
from netCDF4 import Dataset
import FiniteDiff

##
#
# The generation script builds a full-sphere data set
#   at 0.5 degree resolution. It consists of a large
#   collection of eddies of differing sizes, randomly
#   distributed around the globe.
#
##

# Make grid
Nlon, Nlat = int(360*2), int(180*2)

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
eddy_scales = [ 250e3, 750e3, 3500e3 ]
eddy_counts = [ 1e3,   3e2,   5e1 ]
eddy_velocs = [ 0.7,   0.4,   0.1 ]
eddy_kinds  = ['Tor',  'Pot', 'Tor']


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
Psi   = np.zeros( (Nlat,Nlon) )
Phi   = np.zeros( (Nlat,Nlon) )
u_lon = np.zeros( (Nlat,Nlon) )
u_lat = np.zeros( (Nlat,Nlon) )

for scale, count, veloc, kind in zip( eddy_scales, eddy_counts, eddy_velocs, eddy_kinds ):

    for ii in range(int(count)):
        W = ( scale * (1 + 0.05 * np.random.randn() ) ) / 2

        lon0 = np.random.choice( lon )  # choose longitude at random
        lat0 = np.random.choice( lat, p = np.cos(lat) / np.sum( np.cos(lat) ) )   # choose latitude, weighting based on area

        D = dist( lon0, lat0 )

        sign = np.random.choice( (-1,1) )

        if kind == 'Tor':
            Psi += sign * veloc * W * np.exp( -(D/W)**2 )
        elif kind == 'Pot':
            Phi += sign * veloc * W * np.exp( -(D/W)**2 )

# Add in continent info
D2R = np.pi / 180.
for coast_lat, coast_wid in zip( [-65*D2R, 65*D2R], [4*D2R, 4*D2R] ):
    ENV = 0.5 * ( 1 + np.tanh( -1 * np.sign(coast_lat) * (LAT - coast_lat) / coast_wid ) )
    Psi *= ENV
    Phi *= ENV

# Get velocity from streamfunction
ddlon = FiniteDiff.FiniteDiff( lon, 8, Uniform = True, Periodic = True,  Sparse = True )
ddlat = FiniteDiff.FiniteDiff( lat, 8, Uniform = True, Periodic = False, Sparse = True )

u_lon = - ddlat.dot( Psi   )   /   6371e3
u_lat =   ddlon.dot( Psi.T ).T / ( 6371e3 * np.cos(LAT) )

u_lon +=   ddlon.dot( Phi.T ).T / ( 6371e3 * np.cos(LAT) )
u_lat +=   ddlat.dot( Phi   )   /   6371e3

# Save flow to a file
dtype_dim = np.float64
dtype = np.float32

dims = ('time','depth','latitude','longitude')

fill_value = -1e10

with Dataset('velocity_sample.nc', 'w', format='NETCDF4') as fp:

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
    uo_var = fp.createVariable('uo', dtype, dims, contiguous=True, fill_value = fill_value)
    uo_var.scale_factor = 1.
    
    vo_var = fp.createVariable('vo', dtype, dims, contiguous=True, fill_value = fill_value)
    vo_var.scale_factor = 1.
    
    psi_var = fp.createVariable('psi', dtype, dims, contiguous=True, fill_value = fill_value)
    psi_var.scale_factor = 1.
    
    phi_var = fp.createVariable('phi', dtype, dims, contiguous=True, fill_value = fill_value)
    phi_var.scale_factor = 1.

    uo_var[ 0,0,:,:] = u_lon
    vo_var[ 0,0,:,:] = u_lat
    psi_var[0,0,:,:] = Psi
    phi_var[0,0,:,:] = Phi
