import numpy as np
from netCDF4 import Dataset
from scipy.fftpack import fftn, ifftn, fftfreq, fftshift


# Make grid
Nx, Ny = 128,128

Lx, Ly = 400e3, 150e3

dx = Lx / Nx
dy = Ly / Ny

x = np.arange( dx/2, Lx, dx )
y = np.arange( dy/2, Ly, dy ) - Ly/2

X, Y = np.meshgrid(x, y)


# Set some parameters
N_osc = 2
L_jet = 20e3
L_meander = 100e3
W_meander = 4e3

U0 = 1.0
U_pert = 0.2


# Make 'jet'
u_vel = np.zeros( (Ny,Nx) )
v_vel = np.zeros( (Ny,Nx) )

for Ix in range(Nx):
    phase = 2 * np.pi * x[Ix] / L_meander
    
    y0_jet = W_meander * np.sin( phase )
    angle = np.arctan( np.cos( phase ) )
    
    vel = U0 * np.exp( -( (y - y0_jet) / L_jet )**2 )
    
    u_vel[:,Ix] = vel * np.cos(angle)
    v_vel[:,Ix] = vel * np.sin(angle)


def get_random_vel( slope = -1 ):
    
    kx = fftfreq( Nx, d = dx )
    ky = fftfreq( Ny, d = dy )

    KX, KY = np.meshgrid( kx, ky )
    K = np.sqrt( KX**2 + KY**2 )

    phases = 2 * np.pi * np.random.rand( Ny,Nx )

    powers = np.zeros( (Ny,Nx) )
    powers[:,:Nx//2] = K[:,:Nx//2]**slope
    powers[K==0] = 0

    coeffs = np.exp( 1j * phases ) * powers

    # Adjust coefficients to force to be real
    coeffs[Ny//2+1:,Nx//2+1:] = np.conj( coeffs[1:Ny//2,1:Nx//2][::-1,::-1] )
    coeffs[1:Ny//2, Nx//2+1:] = np.conj( coeffs[Ny//2+1:,1:Nx//2][::-1,::-1] )
    coeffs[0,:] = 0.
    coeffs[:,0] = 0.
    
    # Get the corresponding real field, and normalize
    u_pert = ifftn(coeffs).real
    u_pert *= U_pert / np.max(np.abs(u_pert))
    
    return u_pert

# Get perturbation flow
u_pert = get_random_vel()
v_pert = get_random_vel()


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
    lat_dim = fp.createDimension(dim, Ny)
    lat_var = fp.createVariable(dim, dtype_dim, (dim,))
    lat_var[:] = y

    # lon
    dim = 'longitude'
    lon_dim = fp.createDimension(dim, Nx)
    lon_var = fp.createVariable(dim, dtype_dim, (dim,))
    lon_var[:] = x

    #
    uo_var = fp.createVariable('uo', dtype, dims, contiguous=True, fill_value = fill_value)
    uo_var.scale_factor = 1.
    
    vo_var = fp.createVariable('vo', dtype, dims, contiguous=True, fill_value = fill_value)
    vo_var.scale_factor = 1.

    uo_var[0,0,:,:] = u_vel + u_pert
    vo_var[0,0,:,:] = v_vel + v_pert
