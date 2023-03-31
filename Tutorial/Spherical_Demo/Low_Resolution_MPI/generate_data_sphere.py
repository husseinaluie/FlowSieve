import numpy as np
from netCDF4 import Dataset
import FiniteDiff

##
#
# This generation script builds a full-sphere data set
#   at 2-degree resolution. 
#
##

# Make grid
Nlon, Nlat, Ndepth = int(360//2), int(180//2), 48

Llon, Llat, Ldepth = 2 * np.pi, np.pi, 500.

dlon   = Llon / Nlon
dlat   = Llat / Nlat
ddepth = Ldepth / Ndepth

lon = np.arange( dlon/2, Llon, dlon )       - Llon/2
lat = np.arange( dlat/2, Llat, dlat )[1:-1] - Llat/2  # remove the actual poles
Nlat -= 2

depth = np.linspace( 0, Ldepth, Ndepth )

LON, LAT = np.meshgrid(lon, lat)

R_earth = 6371e3
dAreas = R_earth**2 * np.cos(LAT) * dlat * dlon

D2R = np.pi / 180.
R2D = 180. / np.pi

# Set up derivative operators
ddlon = FiniteDiff.FiniteDiff( lon, 8, Uniform = True, Periodic = True,  Sparse = True )
ddlat = FiniteDiff.FiniteDiff( lat, 8, Uniform = True, Periodic = False, Sparse = True )

# Build streamfunction and potential function
Psi_full = np.zeros( (Ndepth, Nlat, Nlon) )
Phi_full = np.zeros( (Ndepth, Nlat, Nlon) )

u_lon_full = np.zeros( (Ndepth, Nlat, Nlon) )
u_lat_full = np.zeros( (Ndepth, Nlat, Nlon) )

for Idepth, dep in enumerate( depth ):
    Psi = np.exp( -( (LAT * R2D - 10 * np.cos( 5 * LON )) / 15)**2 )
    Phi = np.exp( -( (np.abs(LAT * R2D) - 40.) / 10 )**2 ) * np.sin( 2 * LON )

    Psi = Psi - np.mean(Psi)
    Phi = Phi - np.mean(Phi)

    Psi = Psi * ( 1 - 0.5 * np.tanh( (dep - 200.)/25 ) )
    Phi = Phi * np.exp( -(dep/75)**2 )

    Psi_full[ Idepth,:,: ] = Psi
    Phi_full[ Idepth,:,: ] = Phi

    # streamfunction contribution
    u_lon = - ddlat.dot( Psi   )   /   6371e3
    u_lat =   ddlon.dot( Psi.T ).T / ( 6371e3 * np.cos(LAT) )

    # potential function contribution
    u_lon += ddlon.dot( Phi.T ).T / ( 6371e3 * np.cos(LAT) )
    u_lat += ddlat.dot( Phi   )   /   6371e3

    u_lon_full[ Idepth,:,: ] = u_lon
    u_lat_full[ Idepth,:,: ] = u_lat

# re-scale to have max velocity of 0.1 at the surface
rms_vel = np.sqrt( ( u_lon_full[0,:,:]**2 + u_lat_full[0,:,:]**2 ) * dAreas / np.sum(dAreas) )
u_lon *= 0.1 / rms_vel.max()
u_lat *= 0.1 / rms_vel.max()
Psi   *= 0.1 / rms_vel.max()
Phi   *= 0.1 / rms_vel.max()


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
    d_dim = fp.createDimension(dim, Ndepth)
    d_var = fp.createVariable(dim, dtype_dim, (dim,))
    d_var[:] = depth

    # lat
    dim = 'latitude'
    lat_dim = fp.createDimension(dim, Nlat)
    lat_var = fp.createVariable(dim, dtype_dim, (dim,))
    lat_var[:] = lat * R2D

    # lon
    dim = 'longitude'
    lon_dim = fp.createDimension(dim, Nlon)
    lon_var = fp.createVariable(dim, dtype_dim, (dim,))
    lon_var[:] = lon * R2D

    #
    uo_var = fp.createVariable('uo', dtype, dims, contiguous=True, fill_value = fill_value)
    uo_var.scale_factor = 1.
    
    vo_var = fp.createVariable('vo', dtype, dims, contiguous=True, fill_value = fill_value)
    vo_var.scale_factor = 1.
    
    psi_var = fp.createVariable('psi', dtype, dims, contiguous=True, fill_value = fill_value)
    psi_var.scale_factor = 1.

    phi_var = fp.createVariable('phi', dtype, dims, contiguous=True, fill_value = fill_value)
    phi_var.scale_factor = 1.

    uo_var[ 0,:,:,:] = u_lon_full
    vo_var[ 0,:,:,:] = u_lat_full
    psi_var[0,:,:,:] = Psi_full
    phi_var[0,:,:,:] = Phi_full
