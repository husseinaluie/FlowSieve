from netCDF4 import Dataset
import numpy as np

infile  = 'interp.nc'
outfile = 'input.nc'

# Get the surface density
interp = Dataset('interp.nc', 'r')
source = Dataset('source.nc', 'r')
        
time      = interp.variables['time'][:]
depth     = interp.variables['depth'][:]
longitude = interp.variables['longitude'][:]
latitude  = interp.variables['latitude' ][:]

rho = interp.variables['density_interp'][:]
p   = interp.variables['pressure_interp'][:]

uo  = source.variables['uo'][:]
vo  = source.variables['vo'][:]

dims = ('time', 'depth', 'latitude', 'longitude')
    
# Now create a file to run through coarse-graining
with Dataset(outfile, 'w', format='NETCDF4') as fp:
     
    # Create dimension objects
    for dim in dims:
        shape = len(interp.variables[dim])
        dim_var = fp.createDimension(dim, shape)
        var = fp.createVariable(dim, np.float64, (dim,))
        var[:] = interp.variables[dim][:]

    var = fp.createVariable('rho', np.float64, dims, 
            contiguous=True, fill_value = -32767)
    var.scale_factor = 1.
    var[:] = rho

    var = fp.createVariable('p', np.float64, dims,
            contiguous=True, fill_value = -32767)
    var.scale_factor = 1.
    var[:] = p

    var = fp.createVariable('uo', np.float64, dims,
            contiguous=True, fill_value = -32767)
    var.scale_factor = 1.
    var[:] = uo

    var = fp.createVariable('vo', np.float64, dims,
            contiguous=True, fill_value = -32767)
    var.scale_factor = 1.
    var[:] = vo
