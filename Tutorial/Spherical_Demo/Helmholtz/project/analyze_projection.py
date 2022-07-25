import numpy as np

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors

import sys
sys.path.append('../../../../PythonTools')
import PlotTools.ScientificCbar as ScientificCbar

from netCDF4 import Dataset

from pyproj import Proj

# Get the grid and apply a mapping projection
with Dataset('coarsened_sample.nc', 'r') as dset:
    lat = dset['latitude' ][:]
    lon = dset['longitude'][:]
    
LON, LAT = np.meshgrid(lon, lat)

proj = Proj( proj = 'wag7', lon_0 = 0., lat_0 = 0. )
Xp_coarse, Yp_coarse = proj(LON, LAT, inverse=False)


with Dataset('../velocity_sample.nc', 'r') as dset:
    lat = dset['latitude' ][:]
    lon = dset['longitude'][:]
    
LON, LAT = np.meshgrid(lon, lat)

proj = Proj( proj = 'wag7', lon_0 = 0., lat_0 = 0. )
Xp_fine, Yp_fine = proj(LON, LAT, inverse=False)


def get_results( res, kind ):

    ui_proj_fname   = 'projection_ui.nc'      if (res == 'fine') else 'coarse_projection_ui.nc'
    src_fname       = '../velocity_sample.nc' if (res == 'fine') else 'coarsened_sample.nc'

    with Dataset(ui_proj_fname, 'r') as dset:
        u_lon_tor = dset['u_lon_tor'][0,0,:,:]
        u_lat_tor = dset['u_lat_tor'][0,0,:,:]

        u_lon_pot = dset['u_lon_pot'][0,0,:,:]
        u_lat_pot = dset['u_lat_pot'][0,0,:,:]

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

def get_results_uiuj( res = 'fine' ):

    uiuj_proj_fname = 'projection_uiuj.nc'    if (res == 'fine') else 'coarse_projection_uiuj.nc'
    src_fname       = '../velocity_sample.nc' if (res == 'fine') else 'coarsened_sample.nc'

    with Dataset(uiuj_proj_fname, 'r') as dset:
        uu = dset['uu'][0,0,:,:]
        uv = dset['uv'][0,0,:,:]
        vv = dset['vv'][0,0,:,:]

    with Dataset(src_fname, 'r') as src:
        u = src['uo'][0,0,:,:]
        v = src['vo'][0,0,:,:]

    return u, v, uu, uv, vv


for kind in ['u_lon', 'u_lat', 'KE']:

    coarse_tor, coarse_pot, coarse_src, coarse_proj, coarse_diff = get_results( 'coarse', kind )
    fine_tor,   fine_pot,   fine_src,   fine_proj,   fine_diff   = get_results( 'fine', kind )

    cmap = 'viridis' if (kind == 'KE') else 'bwr'
    cv = np.max(np.abs(fine_src))
    norm = colors.LogNorm( vmin = cv / 1e5, vmax = cv ) if (kind == 'KE') else colors.Normalize( vmin = -cv, vmax = cv )

    cv_diff = max( np.max(np.abs( fine_diff ) ), np.max(np.abs( coarse_diff ) ) )
    norm_diff = colors.LogNorm( vmin = cv_diff / 1e5, vmax = cv_diff ) if (kind == 'KE') else  colors.Normalize( vmin = -cv_diff, vmax = cv_diff )

    gridspec_props = dict(wspace = 0.075, hspace = 0.05, left = 0.05, right = 0.95, bottom = 0.05, top = 0.9)

    fig, axes = plt.subplots(5, 2, sharex=True, sharey=True, figsize=(8, 6), gridspec_kw = gridspec_props)

    qms = np.zeros(axes.shape, dtype='object')

    qms[0,0] = axes[0,0].pcolormesh( Xp_coarse, Yp_coarse, coarse_tor,  cmap = cmap, norm = norm )
    qms[1,0] = axes[1,0].pcolormesh( Xp_coarse, Yp_coarse, coarse_pot,  cmap = cmap, norm = norm )
    qms[2,0] = axes[2,0].pcolormesh( Xp_coarse, Yp_coarse, coarse_proj, cmap = cmap, norm = norm )
    qms[3,0] = axes[3,0].pcolormesh( Xp_coarse, Yp_coarse, coarse_src,  cmap = cmap, norm = norm )
    qms[4,0] = axes[4,0].pcolormesh( Xp_coarse, Yp_coarse, coarse_diff, cmap = cmap, norm = norm_diff )

    qms[0,1] = axes[0,1].pcolormesh( Xp_fine, Yp_fine, fine_tor,  cmap = cmap, norm = norm )
    qms[1,1] = axes[1,1].pcolormesh( Xp_fine, Yp_fine, fine_pot,  cmap = cmap, norm = norm )
    qms[2,1] = axes[2,1].pcolormesh( Xp_fine, Yp_fine, fine_proj, cmap = cmap, norm = norm )
    qms[3,1] = axes[3,1].pcolormesh( Xp_fine, Yp_fine, fine_src,  cmap = cmap, norm = norm )
    qms[4,1] = axes[4,1].pcolormesh( Xp_fine, Yp_fine, fine_diff, cmap = cmap, norm = norm_diff )


    cb0 = plt.colorbar( qms[ 0,0], ax = axes[:-1,:] )
    cb1 = plt.colorbar( qms[-1,0], ax = axes[ -1,:] )
    if kind != 'KE':
        ScientificCbar(cb0)
        ScientificCbar(cb1)
            
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

# we've removed the UiUj decomposition part for now, so don't do it.
sys.exit()

###
##
###

for res in ['coarse', 'fine']:

    u, v, uu, uv, vv = get_results_uiuj( res = res )

    Xp = Xp_fine if (res == 'fine') else Xp_coarse
    Yp = Yp_fine if (res == 'fine') else Yp_coarse

    cv = max( np.max( np.abs(u*u) ), np.max(np.abs( u*v ) ), np.max(np.abs( v*v ) ) )
    cmap = 'bwr'
    norm = colors.SymLogNorm( vmin = -cv, vmax = cv, linthresh = cv / 1e2, linscale = 0.2 )
    #norm = colors.Normalize( vmin = -cv, vmax = cv )

    gridspec_props = dict(wspace = 0.075, hspace = 0.05, left = 0.05, right = 0.95, bottom = 0.05, top = 0.9)

    fig, axes = plt.subplots(3, 3, sharex=True, sharey=True, figsize=(9, 6), gridspec_kw = gridspec_props)

    qms = np.zeros(axes.shape, dtype='object')

    qms[0,0] = axes[0,0].pcolormesh( Xp, Yp, u*u,       cmap = cmap, norm = norm )
    qms[0,1] = axes[0,1].pcolormesh( Xp, Yp, uu,        cmap = cmap, norm = norm )
    qms[0,2] = axes[0,2].pcolormesh( Xp, Yp, u*u - uu,  cmap = cmap, norm = norm )

    qms[1,0] = axes[1,0].pcolormesh( Xp, Yp, u*v,       cmap = cmap, norm = norm )
    qms[1,1] = axes[1,1].pcolormesh( Xp, Yp, uv,        cmap = cmap, norm = norm )
    qms[1,2] = axes[1,2].pcolormesh( Xp, Yp, u*v - uv,  cmap = cmap, norm = norm )

    qms[2,0] = axes[2,0].pcolormesh( Xp, Yp, v*v,       cmap = cmap, norm = norm )
    qms[2,1] = axes[2,1].pcolormesh( Xp, Yp, vv,        cmap = cmap, norm = norm )
    qms[2,2] = axes[2,2].pcolormesh( Xp, Yp, v*v - vv,  cmap = cmap, norm = norm )


    cb = plt.colorbar( qms[0,0], ax = axes )
    #ScientificCbar( cb )
            
    for ax in axes.ravel():
        ax.set_xticklabels([])
        ax.set_yticklabels([])
            
    axes[0,0].set_title('Original')
    axes[0,1].set_title('Helmholtz')
    axes[0,2].set_title('Difference')

    axes[0,0].set_ylabel('u*u')
    axes[1,0].set_ylabel('u*v')
    axes[2,0].set_ylabel('v*v')

    plt.savefig( 'uiuj_{0}_projection_results.png'.format(res), dpi = 250)
    plt.close()
