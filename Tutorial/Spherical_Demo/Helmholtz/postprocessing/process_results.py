import numpy as np

import matplotlib
import matplotlib.pyplot as plt

from netCDF4 import Dataset
import glob, sys


# This is a helper function that extracts filter scale from output files 
#   and sort the files accordingly
def sort_files_by_scale(list_of_files):
        
    scales = np.array([Dataset(file, 'r').filter_scale for file in list_of_files])
    inds = np.argsort(scales)
    scales = scales[inds]
    sorted_files = [list_of_files[ind] for ind in inds]
                                    
    return sorted_files, scales

# Get all of the postprocess files and their associated filtering scales
files, space_filters = sort_files_by_scale( glob.glob('postprocess_*km.nc') )


# Also extract the region labels
with Dataset(files[0], 'r') as dset:
    regions = dset['region'][:]

# Get dimensions
Nscales = len(space_filters)
Nregions = len(regions)

# Create storage arrays to hold the variables that we want to plot
coarse_KE = np.zeros(( Nscales, Nregions ))
fine_KE   = np.zeros(( Nscales, Nregions ))
Pi        = np.zeros(( Nscales, Nregions ))


# Now loop through each of the files and pull out the data that we want
for Iscale, fp in enumerate(files):

    with Dataset(fp, 'r') as dset:
        # Recall that dimension order is time - depth - region, so here
        #   we grab the values for the first time, first depth, and all regions
        coarse_KE[Iscale,:] = dset['coarse_KE_area_average'][0,0,:]
        fine_KE[  Iscale,:] = dset['fine_KE_area_average'  ][0,0,:]
        Pi[       Iscale,:] = dset['Pi_area_average'       ][0,0,:]


# Now, simply loop through all of our regions and create a figure!
for Iregion, region in enumerate(regions):

    fig, axes = plt.subplots( 2, 1, figsize = (5,3), sharex = True,
            gridspec_kw = dict( left = 0.1, right = 0.975, bottom = 0.15, top = 0.975 ) )

    axes[0].plot( space_filters / 1e3, coarse_KE[:,Iregion], '-o', label = 'Large-scale' )
    axes[0].plot( space_filters / 1e3, fine_KE[  :,Iregion], '-o', label = 'Small-scale' )
    axes[0].plot( space_filters / 1e3, coarse_KE[:,Iregion] + fine_KE[  :,Iregion], '-o', label = 'Combined' )
    axes[1].plot( space_filters / 1e3, Pi[       :,Iregion], '-o' )

    axes[1].set_xscale('log')
    axes[1].set_xlabel('Filter Scale (km)')

    axes[0].set_ylabel('KE')
    axes[1].set_ylabel('Energy Cascade (Pi)')

    for ax in axes:
        ax.grid(True, which = 'major', linestyle = '--', alpha = 0.7)
        ax.grid(True, which = 'minor', linestyle = ':',  alpha = 0.7)

    # A little bit of string-handling to deal with spaces in our file names
    plt.savefig( region.strip().replace(' ', '_') + '.pdf' )
    plt.close()
