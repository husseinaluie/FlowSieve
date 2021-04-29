import numpy as np

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors

from netCDF4 import Dataset

from PlotTools import ScientificCbar
from pyproj import Proj

# Get the grid and apply a mapping projection
with Dataset('coarsened_sample.nc', 'r') as dset:
    lat = dset['latitude' ][:]
    lon = dset['longitude'][:]
    
LON, LAT = np.meshgrid(lon, lat)

proj = Proj( proj = 'wag7', lon_0 = 0., lat_0 = 0. )
Xp_coarse, Yp_coarse = proj(LON, LAT, inverse=False)


with Dataset('velocity_sample.nc', 'r') as dset:
    lat = dset['latitude' ][:]
    lon = dset['longitude'][:]
    
LON, LAT = np.meshgrid(lon, lat)

proj = Proj( proj = 'wag7', lon_0 = 0., lat_0 = 0. )
Xp_fine, Yp_fine = proj(LON, LAT, inverse=False)


def get_results( res, kind ):

    tor_fname = 'toroidal_projection.nc'  if (res == 'fine') else 'coarse_toroidal_projection.nc'
    pot_fname = 'potential_projection.nc' if (res == 'fine') else 'coarse_potential_projection.nc'
    src_fname = 'velocity_sample.nc'      if (res == 'fine') else 'coarsened_sample.nc'

    with Dataset(tor_fname, 'r') as tor_set:
        u_lon_tor = tor_set['u_lon'][0,0,:,:]
        u_lat_tor = tor_set['u_lat'][0,0,:,:]

    with Dataset(pot_fname, 'r') as pot_set:
        u_lon_pot = pot_set['u_lon'][0,0,:,:]
        u_lat_pot = pot_set['u_lat'][0,0,:,:]

    u_lon_proj = u_lon_tor + u_lon_pot
    u_lat_proj = u_lat_tor + u_lat_pot

    with Dataset(src_fname, 'r') as src:
        u_lon_src = src['uo'][0,0,:,:]
        u_lat_src = src['vo'][0,0,:,:]

    u_lon_diff = u_lon_src - u_lon_proj
    u_lat_diff = u_lat_src - u_lat_proj

    KE_src  = 0.5 * ( u_lon_src**2  + u_lat_src**2  )
    KE_tor  = 0.5 * ( u_lon_tor**2  + u_lat_tor**2  )
    KE_pot  = 0.5 * ( u_lon_pot**2  + u_lat_pot**2  )
    KE_proj = 0.5 * ( u_lon_proj**2 + u_lat_proj**2 )
    KE_diff = 0.5 * ( u_lon_diff**2 + u_lat_diff**2 )

    if kind == 'u_lon':
        return u_lon_tor, u_lon_pot, u_lon_src, u_lon_proj, u_lon_diff
    elif kind == 'u_lat':
        return u_lat_tor, u_lat_pot, u_lat_src, u_lat_proj, u_lat_diff
    elif kind == 'KE':
        return KE_tor, KE_pot, KE_src, KE_proj, KE_diff


for kind in ['u_lon', 'u_lat', 'KE']:

    coarse_tor, coarse_pot, coarse_src, coarse_proj, coarse_diff = get_results( 'coarse', kind )
    fine_tor,   fine_pot,   fine_src,   fine_proj,   fine_diff   = get_results( 'fine', kind )

    cmap = 'viridis' if (kind == 'KE') else 'bwr'
    norm = colors.LogNorm( vmin = fine_src.max() / 1e5, vmax = fine_src.max() ) if (kind == 'KE') else colors.Normalize( vmin = -1, vmax = 1 )

    gridspec_props = dict(wspace = 0.075, hspace = 0.05, left = 0.05, right = 0.95, bottom = 0.05, top = 0.9)

    fig, axes = plt.subplots(5, 2, sharex=True, sharey=True, figsize=(8, 6), gridspec_kw = gridspec_props)

    qms = np.zeros(axes.shape, dtype='object')

    qms[0,0] = axes[0,0].pcolormesh( Xp_coarse, Yp_coarse, coarse_tor,  cmap = cmap, norm = norm )
    qms[1,0] = axes[1,0].pcolormesh( Xp_coarse, Yp_coarse, coarse_pot,  cmap = cmap, norm = norm )
    qms[2,0] = axes[2,0].pcolormesh( Xp_coarse, Yp_coarse, coarse_proj, cmap = cmap, norm = norm )
    qms[3,0] = axes[3,0].pcolormesh( Xp_coarse, Yp_coarse, coarse_src,  cmap = cmap, norm = norm )
    qms[4,0] = axes[4,0].pcolormesh( Xp_coarse, Yp_coarse, coarse_diff, cmap = cmap, norm = norm )

    qms[0,1] = axes[0,1].pcolormesh( Xp_fine, Yp_fine, fine_tor,  cmap = cmap, norm = norm )
    qms[1,1] = axes[1,1].pcolormesh( Xp_fine, Yp_fine, fine_pot,  cmap = cmap, norm = norm )
    qms[2,1] = axes[2,1].pcolormesh( Xp_fine, Yp_fine, fine_proj, cmap = cmap, norm = norm )
    qms[3,1] = axes[3,1].pcolormesh( Xp_fine, Yp_fine, fine_src,  cmap = cmap, norm = norm )
    qms[4,1] = axes[4,1].pcolormesh( Xp_fine, Yp_fine, fine_diff, cmap = cmap, norm = norm )


    cb = plt.colorbar( qms[0,0], ax = axes )
    ScientificCbar(cb)
            
    for ax in axes.ravel():
        ax.set_xticklabels([])
        ax.set_yticklabels([])
            
    axes[0,0].set_title('Coarse Grid')
    axes[0,1].set_title('Full Resolution')

    axes[0,0].set_ylabel('Tor.')
    axes[1,0].set_ylabel('Pot.')
    axes[2,0].set_ylabel('Tor. + Pot.')
    axes[3,0].set_ylabel('Full')
    axes[4,0].set_ylabel('Difference')

    plt.savefig( kind + '_projection_results.png', dpi = 250)
    plt.close()




fig, axes = plt.subplots( 3, 2, sharex = True, sharey = True, figsize = (6,5) )

with Dataset('toroidal_projection.nc', 'r') as tor_set:
    F_tor       = tor_set['F'     ][0,0,:,:]
    F_tor_seed  = tor_set['F_seed'][0,0,:,:]

with Dataset('potential_projection.nc', 'r') as pot_set:
    F_pot       = pot_set['F'     ][0,0,:,:]
    F_pot_seed  = pot_set['F_seed'][0,0,:,:]

with Dataset('coarse_toroidal_projection.nc', 'r') as tor_set:
    F_tor_coarse = tor_set['F'][0,0,:,:]

with Dataset('coarse_potential_projection.nc', 'r') as pot_set:
    F_pot_coarse = pot_set['F'][0,0,:,:]

norm_tor = colors.Normalize( vmin = F_tor.min(), vmax = F_tor.max() )
norm_pot = colors.Normalize( vmin = F_pot.min(), vmax = F_pot.max() )
cmap = 'viridis'

qms = np.zeros(axes.shape, dtype='object')
qms[0,0] = axes[0,0].pcolormesh( Xp_coarse, Yp_coarse, F_tor_coarse,  cmap = cmap, norm = norm_tor )
qms[0,1] = axes[0,1].pcolormesh( Xp_coarse, Yp_coarse, F_pot_coarse,  cmap = cmap, norm = norm_pot )

qms[1,0] = axes[1,0].pcolormesh( Xp_fine, Yp_fine, F_tor_seed,  cmap = cmap, norm = norm_tor )
qms[1,1] = axes[1,1].pcolormesh( Xp_fine, Yp_fine, F_pot_seed,  cmap = cmap, norm = norm_pot )

qms[2,0] = axes[2,0].pcolormesh( Xp_fine, Yp_fine, F_tor,  cmap = cmap, norm = norm_tor )
qms[2,1] = axes[2,1].pcolormesh( Xp_fine, Yp_fine, F_pot,  cmap = cmap, norm = norm_pot )

axes[0,0].set_ylabel('Coarse Result')
axes[1,0].set_ylabel('Fine Seed')
axes[2,0].set_ylabel('Fine Result')

cb0 = plt.colorbar( qms[0,0], ax = axes[:,0] )
cb1 = plt.colorbar( qms[0,1], ax = axes[:,1] )

ScientificCbar(cb0)
ScientificCbar(cb1)

for ax in axes.ravel():
    ax.set_xticklabels([])
    ax.set_yticklabels([])

plt.savefig('Psi_Phi.png', dpi = 350)
plt.close()
