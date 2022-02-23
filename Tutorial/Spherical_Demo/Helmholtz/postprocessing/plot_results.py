import numpy as np
from FiniteDiff import FiniteDiff

import matplotlib
import matplotlib.pyplot as plt

from netCDF4 import Dataset
import glob, sys


with Dataset('results_toroidal.nc', 'r') as dset:
    regions = dset['region'][:]

    ell = dset['ell'][:]

    lat = dset['latitude'][:]

    coarse_KE_region_avg = dset['coarse_KE_area_average'][:,0,0,:]
    coarse_KE_zonal_avg = dset['coarse_KE_zonal_average'][:,0,0,:]

    fine_KE_region_avg = dset['fine_KE_area_average'][:,0,0,:]
    fine_KE_zonal_avg  = dset['fine_KE_zonal_average'][:,0,0,:]

# region = "              Global", "             Tropics",
#          "    North of Tropics", "    South of Tropics", "  25th KE Percentile",
#          "  50th KE Percentile", "  75th KE Percentile" ;

ddk = FiniteDiff( 1. / ell, 4, Sparse = True, Uniform = False, Periodic = False )


## Plot the area-averaged fields

fig, axes = plt.subplots( 2, 1, figsize = (4,3.5), sharex = True, gridspec_kw = dict( left = 0.15, bottom = 0.1, right = 0.95, top = 0.95, hspace = 0.05 ) )

axes[0].plot( ell / 1e3,          coarse_KE_region_avg[:,4:], marker = '.', markersize = 1  )
axes[1].plot( ell / 1e3, ddk.dot( coarse_KE_region_avg[:,4:] ) / 1e7, marker = '.', markersize = 1 )
#axes[2].plot( ell / 1e3, np.tile( (1./ell).reshape((len(ell),1)), (1,3)) * ddk.dot( coarse_KE_region_avg[:,4:] ), marker = '.', markersize = 1 )

for ax in axes:
    ax.set_ylim( ax.get_ylim() )
    ax.plot( [ 250, 250], ax.get_ylim(), '--k' )
    ax.plot( [ 750, 750], ax.get_ylim(), '--k' )
    ax.plot( [3500,3500], ax.get_ylim(), '--k' )

    ax.set_xlim( 40, ell.max() / 1e3 )
    ax.set_xscale('log')

    ax.grid( which = 'major', alpha = 0.8, linestyle = '-' )

axes[0].legend( ['25th Percentile', '50th Percentile', '75th Percentile'] )

axes[-1].set_xlabel('$\ell$ [km]')

axes[0].set_ylabel('Cumulative Coarse KE\nDensity [J/m$^3$]')
axes[1].set_ylabel('Spectral Coarse KE\nDensity [$10^7$J/m$^2$]')

plt.savefig('coarse_KE_spectra.pdf')
plt.close()




## Plot the zonally-averaged fields

fig, axes = plt.subplots( 2, 1, figsize = (4,4), sharex = True )

qm0 = axes[0].pcolormesh( ell / 1e3, lat,          coarse_KE_zonal_avg.T,   cmap = 'plasma' )
qm1 = axes[1].pcolormesh( ell / 1e3, lat, ddk.dot( coarse_KE_zonal_avg ).T, cmap = 'plasma' )

for ax in axes:
    ax.set_ylim( ax.get_ylim() )
    ax.plot( [ 250, 250], ax.get_ylim(), '--w' )
    ax.plot( [ 750, 750], ax.get_ylim(), '--w' )
    ax.plot( [3500,3500], ax.get_ylim(), '--w' )

    ax.set_xscale('log')

    ax.grid( which = 'major', alpha = 0.8, linestyle = '-' )

plt.colorbar( qm0, ax = axes[0] )
plt.colorbar( qm1, ax = axes[1] )

plt.savefig('coarse_KE_zonal_avg.png', dpi = 350)
plt.close()
