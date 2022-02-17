import numpy as np
import argparse
import sys
from scipy.interpolate import interp1d
from netCDF4 import Dataset

parser = argparse.ArgumentParser(description='Merge downsampling results.')

parser.add_argument('--input_files', metavar='file1.nc file2.nc', type=str, nargs='+', required = True,
        help='Filenames to downsampling results to be merged. e.g. "results_1.nc results_4.nc"')

parser.add_argument('--output_filename', metavar='output', type=str, nargs=1, required = True,
        help='Filename pattern for merged set. e.g. "postprocess.nc"')

parser.add_argument('--print_level', metavar='debug', type=int, nargs=1, default = [0,],
        help='String to indicate how much printing to do. Options are 0, 1, 2 [higher value means more printed].')

args = parser.parse_args()

if type(args.print_level) == type(0):
    print_level = args.print_level
elif type(args.print_level) == type([0,]):
    print_level = args.print_level[0]

print("Attempting to merge {0} into output file {1}.".format( args.input_files, args.output_filename[0] ))

## We need to determine the resolution order for the downsampling files
Nlats = [ len(Dataset(fp, 'r')['latitude'][:]) for fp in args.input_files ]
res_order = np.argsort( Nlats )[::-1]

sorted_inputs = [ args.input_files[ind] for ind in res_order ]

print("Determined the resolution ordering (highest to lowest) to be: {0}".format( sorted_inputs ))
with Dataset( sorted_inputs[0], 'r' ) as dset:
    finest_latitude_grid = dset['latitude'][:]
    finest_longitude_grid = dset['longitude'][:]


## Get a list of the dimensions and variables
with Dataset( args.input_files[0], 'r' ) as dset:
    Nvars = len( dset.variables.copy() )
    Ndims = len( dset.dimensions.copy() )

print("  Identified {0:g} dimensions and {1:g} variables to copy into merged file.".format( Ndims, Nvars ), flush = True)


## Get the filter scale from each file, and sort them in ascending order
ell_sets = []
for fp in sorted_inputs:
    with Dataset( fp, 'r' ) as dset:
        ell_sets = ell_sets + [ dset['ell'][:], ]
ells = np.unique( np.concatenate( ell_sets ) )

print("  Identified {0:g} unique filter scales.".format( len(ells) ), flush = True)

## Create the output file
with Dataset( args.output_filename[0], 'w', format='NETCDF4') as out_fp:

    dtype     = np.float32
    dtype_dim = np.float64

    # Create ell dimension
    dim = 'ell'
    ell_dim = out_fp.createDimension(dim, len(ells))
    ell_var = out_fp.createVariable(dim, dtype_dim, (dim,))
    ell_var[:] = ells

    # Reproduce all previously-existing dimensions, except for ell (since we're merging on that dimension)
    with Dataset( sorted_inputs[0], 'r' ) as dset:
        all_dims = dset.dimensions.copy()
        for dim in all_dims:
            if not( dim == 'ell' ):
                if print_level >= 0:
                    print("  .. copying dimension " + dim, flush = True)
                dim_dim = out_fp.createDimension( dim, all_dims[dim].size )

    # Now loop through all previously-existing variables, and re-create them will an ell-dimension prepended.
    #   The iterate through all files and copy in the data to the new file.
    print("  Preparing to copy variables...")
    nc_var_objs = dict()
    with Dataset( sorted_inputs[0], 'r' ) as dset:
        all_vars = dset.variables.copy()
        for varname in all_vars:
            if not( varname == 'ell' ):

                # Extract the dimensions (in order) for varname
                var_dims = all_vars[varname].dimensions

                # Get the variable type
                var_dtype = all_vars[varname].dtype

                if print_level >= 1:
                    print("  .. initializing variable " + varname + " with dimensions {0}".format(var_dims), flush = True)

                # Create netcdf object for variable
                nc_var_objs[varname] = out_fp.createVariable( varname, var_dtype, var_dims )

                # Copy over any attributes that the variable had, except for factor/offset info
                var_attrs = all_vars[varname].ncattrs()
                for attr in var_attrs:
                    if not( attr in ['scale_factor', 'add_offset'] ):
                        nc_var_objs[varname].setncattr(attr, all_vars[varname].getncattr( attr ) )

        
    # Copy of the non-ell dimensions from the highest resolution grid (since those aren't changing)
    with Dataset( sorted_inputs[0], 'r' ) as in_dset:
        for Ivar,varname in enumerate(all_dims):
            if not( varname == 'ell' ):
                nc_var_objs[varname][:] = in_dset[varname][:]

    # Now, just iterate through the ells and copy in the data from the highest resolution option
    if print_level >= 1:
        print("  Copying variable data", flush = True)

    prog = 5
    for Iell, ell in enumerate( ells ):

        if print_level == 0:
            while 100. * Iell / float(len(ells)) >= prog:
                if prog == 5:
                    print("  .", end = '', flush = True)
                elif prog % 25 == 0:
                    print(" | ", end = '', flush = True)
                else:
                    print(".", end = '', flush = True)
                prog += 5

        # Find the first (highest resolution) dataset that has ell
        has_ell = [ (ell in ell_set) for ell_set in ell_sets ]
        set_index = np.where(has_ell)[0][0]

        Iell_in = np.where( ell == ell_sets[set_index] )[0][0]

        if print_level >= 1:
            print("    Pulling ell = {0:g} from {1}".format(ell, sorted_inputs[set_index]), end = '' )
            print("  .. out_index is {0:d} or {1:d} and input_index is {2:d}".format( Iell+1, len(ells), Iell_in ) )
        
        with Dataset( sorted_inputs[set_index], 'r' ) as in_dset:
            for Ivar,varname in enumerate(all_vars):
                if not( varname in all_dims):

                    var_dims = in_dset.variables[varname].dimensions

                    data_to_copy = in_dset[varname][Iell_in,:]

                    if ( 'latitude' in var_dims ):
                        if print_level >= 2:
                            if (Iell == 0):
                                print("    .. {0} has latitude dependence - will interpolate onto finest grid.".format(varname) )
                        # If there is latitude dependence, we jeed to interpolated onto the finest grid
                        lat_axis = var_dims.index('latitude')
                        interp_func = interp1d( in_dset['latitude'][:], data_to_copy, axis = lat_axis-1, 
                                                    kind = 'nearest', fill_value = 'extrapolate' )
                        data_to_copy = interp_func( finest_latitude_grid.data )

                    if ( 'longitude' in var_dims ):
                        if print_level >= 2:
                            if (Iell == 0):
                                print("    .. {0} has longitude dependence - will interpolate onto finest grid.".format(varname) )
                        # If there is longitude dependence, we jeed to interpolated onto the finest grid
                        lon_axis = var_dims.index('longitude')
                        interp_func = interp1d( in_dset['longitude'][:], data_to_copy, axis = lon_axis-1, 
                                                    kind = 'nearest', fill_value = 'extrapolate' )
                        data_to_copy = interp_func( finest_longitude_grid.data )

                    # Once all interpolations are done, write to the output
                    nc_var_objs[varname][Iell,:] = data_to_copy

if print_level == 0:
    print('\n')
print("Done.")
