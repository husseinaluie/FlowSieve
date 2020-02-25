from netCDF4 import Dataset
import numpy as np

infile  = 'input.nc'
outfile = 'subset.nc'

# Get the surface density
source = Dataset(infile, 'r')
        
time      = source.variables['time'][:]
depth     = source.variables['depth'][:]
longitude = source.variables['longitude'][:]
latitude  = source.variables['latitude' ][:]

ref_var = 'uo'
dims  = source.variables[ref_var].dimensions
dtype = source.variables[ref_var].datatype

time_lb   = 0
time_ub   = 1#len(time)
time_step = 1

#lon_lb = np.argmin(np.abs(longitude - (-90)))
#lon_ub = np.argmin(np.abs(longitude - ( 60)))
lon_lb = 0
lon_ub = len(longitude)

#lat_lb = np.argmin(np.abs(latitude - (-20)))
#lat_ub = np.argmin(np.abs(latitude - ( 50)))
lat_lb = 0
lat_ub = len(latitude)
    
# Now create a file to run through coarse-graining
with Dataset(outfile, 'w', format='NETCDF4') as fp:
     
    # Create dimension objects
    for dim in source.variables[ref_var].dimensions:
        if dim == 'longitude':
            shape = len(longitude[lon_lb:lon_ub])
        elif dim == 'latitude':
            shape = len(latitude[lat_lb:lat_ub])
        elif dim == 'time':
            shape = len(time[time_lb:time_ub:time_step])
        else:
            shape = len(source.variables[dim])
        dim_var = fp.createDimension(dim, shape)

    for var in source.variables:
        if var in dims:
            var_var = fp.createVariable(var, dtype, source.variables[var].dimensions)
        else:
            var_var = fp.createVariable(var, dtype, source.variables[var].dimensions,
                    fill_value = -32767, contiguous = True)
            var_var.scale_factor = 1.

        if var == 'longitude':
            var_var[:] = longitude[lon_lb:lon_ub]
        elif var == 'latitude':
            var_var[:] = latitude[lat_lb:lat_ub]
        elif var == 'time':
            var_var[:] = time[time_lb:time_ub:time_step]
        elif var == 'depth':
            var_var[:] = depth
        else:
            var_var[:] = source.variables[var][
                    time_lb:time_ub:time_step,
                    :,
                    lat_lb:lat_ub,
                    lon_lb:lon_ub]

