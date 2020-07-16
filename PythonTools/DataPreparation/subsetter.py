from netCDF4 import Dataset
import numpy as np

#infile  = 'filter_50km_red.nc'
#outfile = 'hourly_50km.nc'

infile  = '../coarse_input.nc'
outfile = 'coarse_in_subs.nc'

# Get the surface density
source = Dataset(infile, 'r')
        
time      = source.variables['time'][:]
depth     = source.variables['depth'][:]
longitude = source.variables['longitude'][:]
latitude  = source.variables['latitude' ][:]

ref_var = 'uo'
dims  = source.variables[ref_var].dimensions
#dtype = source.variables[ref_var].datatype
dtype = np.float32

time_lb   = 0
time_ub   = 4#len(time)
time_step = 1

lon_lb   = 0
lon_ub   = len(longitude)
lon_step = 1

lat_lb   = 0
lat_ub   = len(latitude)
lat_step = 1
    
# Now create a file to run through coarse-graining
with Dataset(outfile, 'w', format='NETCDF4') as fp:
     
    # Create dimension objects
    for dim in source.variables[ref_var].dimensions:
        if dim == 'longitude':
            shape = len(longitude[lon_lb:lon_ub:lon_step])
        elif dim == 'latitude':
            shape = len(latitude[ lat_lb:lat_ub:lat_step])
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
            var_var[:] = longitude[lon_lb:lon_ub:lon_step]
        elif var == 'latitude':
            var_var[:] = latitude[ lat_lb:lat_ub:lat_step]
        elif var == 'time':
            var_var[:] = time[time_lb:time_ub:time_step]
        elif var == 'depth':
            var_var[:] = depth
        else:
            out_ind = 0
            for in_ind in range(time_lb, time_ub, time_step):
                var_var[out_ind,:,:,:] = source.variables[var][
                        in_ind,
                        :,
                        lat_lb:lat_ub:lat_step,
                        lon_lb:lon_ub:lon_step]
                out_ind += 1

