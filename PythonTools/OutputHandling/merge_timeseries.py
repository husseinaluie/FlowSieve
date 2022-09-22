import numpy as np
import argparse
import glob
from netCDF4 import Dataset

parser = argparse.ArgumentParser(description='Merge timeseries.')

parser.add_argument('--file_pattern', metavar='search_pattern', type=str, nargs=1, required = True,
        help='Filename pattern for merging. e.g. "postprocess_*.nc"')

parser.add_argument('--output_filename', metavar='output', type=str, nargs=1, required = True,
        help='Filename pattern for merged set. e.g. "postprocess.nc"')

parser.add_argument('--print_level', metavar='debug', type=int, nargs=1, default = 0,
        help='String to indicate how much printing to do. Options are 0, 1, 2 [higher value means more printed].')

# Pass --exclude_time to exclude time means
parser.add_argument('--exclude_time', dest='copy_time_means', action='store_false',
        help="Pass '--exclude_time' to have time means ignored when merging postprocessing files.")
parser.set_defaults(copy_time_means=True)

# Pass --exclude_OkuboWeiss to exclude OkuboWeiss histograms
parser.add_argument('--exclude_OkuboWeiss', dest='copy_OkuboWeiss', action='store_false',
        help="Pass '--exclude_OkuboWeiss' to have OkuboWeiss histograms ignored when merging postprocessing files.")
parser.set_defaults(copy_OkuboWeiss=True)

# Pass --exclude_zonal_mean to exclude zonal means
parser.add_argument('--exclude_zonal_means', dest='copy_zonal_means', action='store_false',
        help="Pass '--exclude_zonal_means' to have zonal means ignored when merging postprocessing files.")
parser.set_defaults(copy_zonal_means=True)

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



## First, find all files matching the requested pattern
all_fps = glob.glob( args.file_pattern[0] )
Nfps = len(all_fps)
print("  Identified {0:g} files for merging.".format(Nfps), flush = True)


## Get a list of the dimensions and variables
with Dataset( all_fps[0], 'r' ) as dset:
    Nvars = len( dset.variables.copy() )
    Ndims = len( dset.dimensions.copy() )

print("  Identified {0:g} dimensions and {1:g} variables to copy into merged file.".format( Ndims, Nvars ), flush = True)


## Get the time from each file, and sort them in ascending order
#   for the sorting, we will assume that the it's enough to sort
#   the files, and not that we need to sort within or between files.
##
t0s = np.zeros( Nfps )
Nts = np.zeros( Nfps, dtype = np.int )
for Ifp, fp in enumerate(all_fps):
    with Dataset( fp, 'r' ) as dset:
        t0s[Ifp] = dset['time'][0]
        Nts[Ifp] = len( dset['time'] )

t0s_sort_inds = np.argsort( t0s )
sorted_fps = [ all_fps[ind] for ind in t0s_sort_inds ]
Ntime = np.sum( Nts )


## Create the output file
with Dataset( args.output_filename[0], 'w', format='NETCDF4') as out_fp:

    dtype     = np.float32
    dtype_dim = np.float64

    # Reproduce all previously-existing dimensions EXCEPT TIME
    with Dataset( all_fps[0], 'r' ) as dset:
        all_dims = dset.dimensions.copy()
        for dim in all_dims:
            if ( 'time' == dim ):
                continue
            if ( (DO_NOT_COPY_OKUBOWEISS ) and ('OkuboWeiss' == dim) ):
                continue
            if print_level >= 0:
                print("  .. copying dimension " + dim, flush = True)
            dim_dim = out_fp.createDimension( dim, all_dims[dim].size )

    # Create new time dimension
    dim = 'time'
    time_dim = out_fp.createDimension(dim, Ntime)
    time_var = out_fp.createVariable(dim, dtype_dim, (dim,))

    # Now loop through all previously-existing variables, and re-create them 
    #   Then iterate through all files and copy in the data to the new file.
    print("  Preparing to copy variables...")
    nc_var_objs = dict()
    with Dataset( all_fps[0], 'r' ) as dset:
        all_vars = dset.variables.copy()
        for varname in all_vars:

            # Skip time
            if (varname == 'time'):
                continue

            # Don't include time averages
            if ( (DO_NOT_COPY_TIME_MEANS ) and ('time_average' in varname) ):
                continue

            # Don't include OkuboWeiss histograms
            if ( (DO_NOT_COPY_OKUBOWEISS ) and ('OkuboWeiss' in varname) ):
                continue

            # Don't include zonal means
            if ( (DO_NOT_COPY_ZONAL_MEANS ) and ('zonal_average' in varname) ):
                continue

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

        # Don't include zonal means
        if ( (DO_NOT_COPY_ZONAL_MEANS ) and ('zonal_average' in varname) ):
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
            # If it's a variable
            
            # We need to do some unpleasant stuff to handle time being
            # anywhere in the dimension ordering *ugh*
            # the 'indices' variable contains that magic

            with Dataset( all_fps[0], 'r' ) as dset:
                all_vars = dset.variables.copy()
                vardims = all_vars[varname].dimensions
                time_dim_index = vardims.index('time')

            time_ind_0 = 0
            for Ifp, fp in enumerate(sorted_fps):
                time_ind_1 = time_ind_0 + Nts[Ifp]
                indices = tuple([ slice(None) for null in np.arange(time_dim_index) ]) + ( slice(time_ind_0,time_ind_1), Ellipsis )

                # Now do the actual copy
                with Dataset( fp, 'r' ) as in_dset:
                    nc_var_objs[varname][indices] = in_dset[varname][Ellipsis]

                # And increment time
                time_ind_0 += Nts[Ifp]
        elif (varname != 'time'):
            # If it's a dimension
            dims = var_dims
            with Dataset( all_fps[0], 'r' ) as in_dset:
                nc_var_objs[varname][:] = in_dset[varname][:]
        else:
            # if it's time dimension
            time_ind_0 = 0
            for Ifp, fp in enumerate(sorted_fps):
                time_ind_1 = time_ind_0 + Nts[Ifp]
                with Dataset( fp, 'r' ) as in_dset:
                    time_var[time_ind_0:time_ind_1] = in_dset['time'][:]
                time_ind_0 += Nts[Ifp]

if print_level >= 2:
    print("Done.")
else:
    print("\nDone.")
