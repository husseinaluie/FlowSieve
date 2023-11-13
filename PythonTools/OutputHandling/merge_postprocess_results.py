import numpy as np
import argparse
import glob
from netCDF4 import Dataset

parser = argparse.ArgumentParser(description='Merge postprocessing.')

parser.add_argument('--file_pattern', metavar='search_pattern', type=str, nargs=1, required = True,
        help='Filename pattern for merging. e.g. "postprocess_*.nc"')

parser.add_argument('--output_filename', metavar='output', type=str, nargs=1, required = True,
        help='Filename pattern for merged set. e.g. "postprocess.nc"')

parser.add_argument('--print_level', metavar='debug', type=int, nargs=1, default = 0,
        help='String to indicate how much printing to do. Options are 0, 1, 2 [higher value means more printed].')

# Pass --exclude_time to exclude time means
parser.add_argument('--exclude_time_means', dest='copy_time_means', action='store_false',
        help="Pass '--exclude_time_means' to have time means ignored when merging postprocessing files.")
parser.set_defaults(copy_time_means=True)

# Pass --exclude_OkuboWeiss to exclude OkuboWeiss histograms
parser.add_argument('--exclude_OkuboWeiss', dest='copy_OkuboWeiss', action='store_false',
        help="Pass '--exclude_OkuboWeiss' to have OkuboWeiss histograms ignored when merging postprocessing files.")
parser.set_defaults(copy_OkuboWeiss=True)

# Pass --exclude_zonal_mean to exclude zonal means
parser.add_argument('--exclude_zonal_means', dest='copy_zonal_means', action='store_false',
        help="Pass '--exclude_zonal_means' to have zonal means ignored when merging postprocessing files.")
parser.set_defaults(copy_zonal_means=True)

# Pass --exclude_coarse_maps to exclude coarsened maps
parser.add_argument('--exclude_coarse_maps', dest='copy_coarse_maps', action='store_false',
        help="Pass '--exclude_coarse_maps' to have coarsened_maps ignored when merging postprocessing files.")
parser.set_defaults(copy_coarse_maps=True)

# Now actually parse in command-line flags
args = parser.parse_args()

if type(args.print_level) == type(0):
    print_level = args.print_level
elif type(args.print_level) == type([0,]):
    print_level = args.print_level[0]

print("Attempting to merge all files matching pattern {0} into output file {1}.".format( args.file_pattern[0], args.output_filename[0] ))

DO_NOT_COPY_TIME_MEANS = not( args.copy_time_means )
if DO_NOT_COPY_TIME_MEANS:
    print("  Will not merge time means maps.")

DO_NOT_COPY_OKUBOWEISS = not( args.copy_OkuboWeiss )
if DO_NOT_COPY_OKUBOWEISS:
    print("  Will not merge OkuboWeiss histograms.")

DO_NOT_COPY_ZONAL_MEANS = not( args.copy_zonal_means )
if DO_NOT_COPY_ZONAL_MEANS:
    print("  Will not merge zonal means.")

DO_NOT_COPY_COARSE_MAPS = not( args.copy_coarse_maps )
if DO_NOT_COPY_COARSE_MAPS:
    print("  Will not merge coarsened maps.")



## First, find all files matching the requested pattern
all_fps = glob.glob( args.file_pattern[0] )
Nfps = len(all_fps)
print("  Identified {0:g} files for merging.".format(Nfps), flush = True)


## Get a list of the dimensions and variables
with Dataset( all_fps[0], 'r' ) as dset:
    Nvars = len( dset.variables.copy() )
    Ndims = len( dset.dimensions.copy() )

print("  Identified {0:g} dimensions and {1:g} variables to copy into merged file.".format( Ndims, Nvars ), flush = True)


## Get the filter scale from each file, and sort them in ascending order
ells = np.zeros( Nfps )
for Ifp, fp in enumerate(all_fps):
    ells[Ifp] = Dataset( fp, 'r' ).filter_scale

ell_sort_inds = np.argsort( ells )
ells = ells[ ell_sort_inds ]
sorted_fps = [ all_fps[ind] for ind in ell_sort_inds ]


## Create the output file
with Dataset( args.output_filename[0], 'w', format='NETCDF4') as out_fp:

    dtype     = np.float32
    dtype_dim = np.float64

    # Create ell dimension
    dim = 'ell'
    ell_dim = out_fp.createDimension(dim, Nfps)
    ell_var = out_fp.createVariable(dim, dtype_dim, (dim,))
    ell_var[:] = ells

    # Reproduce all previously-existing dimensions
    with Dataset( all_fps[0], 'r' ) as dset:
        all_dims = dset.dimensions.copy()
        for dim in all_dims:
            if ( (DO_NOT_COPY_OKUBOWEISS ) and ('OkuboWeiss' == dim) ):
                continue
            if print_level >= 0:
                print("  .. copying dimension " + dim, flush = True)
            dim_dim = out_fp.createDimension( dim, all_dims[dim].size )

    # Now loop through all previously-existing variables, and re-create them will an ell-dimension prepended.
    #   The iterate through all files and copy in the data to the new file.
    print("  Preparing to copy variables...")
    nc_var_objs = dict()
    with Dataset( all_fps[0], 'r' ) as dset:
        all_vars = dset.variables.copy()
        for varname in all_vars:

            # Don't include time averages
            if ( (DO_NOT_COPY_TIME_MEANS ) and ('time_average' in varname) ):
                continue

            # Don't include OkuboWeiss histograms
            if ( (DO_NOT_COPY_OKUBOWEISS ) and ('OkuboWeiss' in varname) ):
                continue

            # Don't include time averages
            if ( (DO_NOT_COPY_ZONAL_MEANS ) and ('zonal_average' in varname) ):
                continue

            # Don't include time averages
            if ( (DO_NOT_COPY_COARSE_MAPS ) and ('coarsened_map' in varname) ):
                continue

            # Extract the dimensions (in order) for varname
            var_dims = all_vars[varname].dimensions

            # Prepend ell dimension to others, if not a dimension variable
            if not( varname in all_dims):
                dims = ( 'ell', *var_dims )
            else:
                dims = var_dims

            # Get the variable type
            var_dtype = all_vars[varname].dtype

            if print_level >= 1:
                print("  .. initializing variable " + varname + " with dimensions {0}".format(dims), flush = True)

            # Create netcdf object for variable
            params = dict()
            if '_FillValue' in all_vars[varname].ncattrs():
                params['fill_value'] = all_vars[varname].getncattr('_FillValue')
            nc_var_objs[varname] = out_fp.createVariable( varname, var_dtype, dims, **params )

            # Copy over any attributes that the variable had, except for factor/offset info
            var_attrs = all_vars[varname].ncattrs()
            for attr in var_attrs:
                if not( attr in ['scale_factor', 'add_offset', '_FillValue'] ):
                    nc_var_objs[varname].setncattr(attr, all_vars[varname].getncattr( attr ) )


    # Now all that remains is to iterate through the files and copy over the variable information
    #   again, handling dimension variables differently
    if print_level >= 1:
        print("  Copying variable data", flush = True)
    prog = 5
    for Ivar,varname in enumerate(all_vars):

        # Don't include time averages
        if ( (DO_NOT_COPY_TIME_MEANS ) and ('time_average' in varname) ):
            continue

        # Don't include OkuboWeiss histograms
        if ( (DO_NOT_COPY_OKUBOWEISS ) and ('OkuboWeiss' in varname) ):
            continue

        # Don't include time averages
        if ( (DO_NOT_COPY_ZONAL_MEANS ) and ('zonal_average' in varname) ):
            continue

        # Don't include time averages
        if ( (DO_NOT_COPY_COARSE_MAPS ) and ('coarsened_map' in varname) ):
            continue

        if print_level >= 2:
            print("  .. .. copying data for " + varname, flush = True)
        else:
            while 100 * Ivar / Nvars >= prog:
                if prog == 5:
                    print("  .", end = '', flush = True)
                elif prog % 25 == 0:
                    print(" | ", end = '', flush = True)
                else:
                    print(".", end = '', flush = True)
                prog += 5
            
        if not( varname in all_dims):
            for Ifp, fp in enumerate(sorted_fps):
                with Dataset( fp, 'r' ) as in_dset:
                    nc_var_objs[varname][Ifp,:] = in_dset[varname][:]
        else:
            dims = var_dims
            with Dataset( all_fps[0], 'r' ) as in_dset:
                nc_var_objs[varname][:] = in_dset[varname][:]

if print_level >= 2:
    print("Done.")
else:
    print("\nDone.")
