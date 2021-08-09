import numpy as np
from netCDF4 import Dataset
import glob

####
##
####

# This is a helper function that extracts filter scale from output files 
#   and sort the files accordingly
def sort_files_by_scale(list_of_files):
        
    scales = np.zeros( len(list_of_files) )
    for ind, fp in enumerate( list_of_files ):
        scales[ind] = Dataset(fp, 'r').filter_scale
    inds = np.argsort(scales)
    scales = scales[inds]
    sorted_files = [list_of_files[ind] for ind in inds]
                                    
    return sorted_files, scales

###
##
###

fields_to_plot = ['Pi_area_average', 'Pi_Helm_area_average']

# Read in from one file to get region and scale information
files_tor, ells = sort_files_by_scale( glob.glob('postprocess_toroidal_*km.nc') )
files_pot, ells = sort_files_by_scale( glob.glob('postprocess_potential_*km.nc') )
files_tot, ells = sort_files_by_scale( glob.glob('postprocess_full_*km.nc') )

with Dataset(files_tor[0], 'r') as dset:
    time = dset['time'][:]

Ntime = len( time )
Nell  = len( ells )

dataset = {}

Pi_tor  = np.zeros( Nell )
Pi_pot  = np.zeros( Nell )
Pi_tot  = np.zeros( Nell )
Pi_Helm = np.zeros( Nell )

with Dataset( files_tor[0], 'r' ) as dset:
    water_areas = dset['region_areas_water_only'][0,0,0]
    total_areas = dset['region_areas'           ][0,0,0]

    area_adjustment = total_areas / water_areas
print(area_adjustment)

for Iell, fp in enumerate(files_tot):
    with Dataset(fp, 'r') as dset:
        Pi_tot[Iell,:] = dset['Pi_area_average'][0,0,0] * area_adjustment

for Iell, fp in enumerate(files_pot):
    with Dataset(fp, 'r') as dset:
        Pi_pot[Iell,:] = dset['Pi_area_average'][0,0,0] * area_adjustment

for Iell, fp in enumerate(files_tor):
    with Dataset(fp, 'r') as dset:
        Pi_tor[ Iell,:] = dset['Pi_area_average'     ][0,0,0] * area_adjustment
        Pi_Helm[Iell,:] = dset['Pi_Helm_area_average'][0,0,0] * area_adjustment

###
##      Create plots
###

print('Done, ready to plot.')
import matplotlib.pyplot as plt

fig, ax = plt.subplots( 1, 1, figsize = (5,2.5), gridspec_kw = dict( left = 0.075, right = 0.95, bottom = 0.075, top = 0.95 ) )

ax.plot( 1e3 / ells, Pi_tor,  label = 'Toroidal' )
ax.plot( 1e3 / ells, Pi_pot,  label = 'Potential' )
ax.plot( 1e3 / ells, Pi_tot,  label = 'Old bar(uiuj)' )
ax.plot( 1e3 / ells, Pi_Helm, label = 'Helm. bar(uiuj)' )

ax.legend( loc = 'best', ncol = 2 )

ax.set_xscale('log')
ax.grid( which = 'major', linestyle = '--', alpha = 0.7 )
ax.grid( which = 'minor', linestyle = ':',  alpha = 0.4 )

ax.set_xlabel('$\ell^{-1}$ [km$^{-1}$]')

plt.savefig('Pi_spectra.pdf')
plt.close()
