import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import cmocean, sys
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
import PlotTools

dpi = PlotTools.dpi

fp = 'filter_output.nc'
results = Dataset(fp, 'r')
source  = Dataset('input.nc', 'r')

R_earth = 6371e3

D2R = np.pi / 180
R2D = 180 / np.pi

eps = 1e-10

# Get the grid from the first filter
latitude  = results.variables['latitude'][:] * R2D
longitude = results.variables['longitude'][:] * R2D
scales    = results.variables['scale'][:]
depth     = results.variables['depth'][:]
time      = results.variables['time'][:]
mask      = results.variables['mask'][:]
transfer  = results.variables['energy_transfer'][:, 0, 0, :, :]
if 'baroclinic_transfer' in results.variables:
    bc_transfer = results.variables['baroclinic_transfer'][:, 0, 0, :, :]
    vort_r      = results.variables['vort_r'][:, 0, 0, :, :]
    rho         = source.variables['rho'][0, 0, :, :]

uo = source.variables['uo'][0, 0, :, :]
vo = source.variables['vo'][0, 0, :, :]
Full_KE = 0.5 * (uo**2 + vo**2)

num_scales = len(scales)-1

LON, LAT = np.meshgrid(longitude * D2R, latitude * D2R)

Nlat = len(latitude)
Nlon = len(longitude)

map_settings = PlotTools.MapSettings(longitude, latitude)

meridians = np.round(np.linspace(longitude.min(), longitude.max(), 5))
parallels = np.round(np.linspace(latitude.min(),  latitude.max(),  5))

cbar_props     = dict(pad = 0.1, shrink = 0.85, orientation = 'vertical')
gridspec_props = dict(wspace = 0.05, hspace = 0.05, left = 0.02, right = 0.98, bottom = 0.02, top = 0.98)

##
## Begin Plotting
##

# Pi vs l^2 * Lambda^m
if True:

    # Initialize figure
    fig, axes = plt.subplots(2, 2, 
            gridspec_kw = dict(left = 0.1, right = 0.95, bottom = 0.1, top = 0.95,
                hspace=0.02, wspace=0.02),
            figsize=(7.5, 6) )

    x_data =    transfer[:-1,:,:] * 1e6
    y_data = bc_transfer[:-1,:,:] * 1e6

    tiled_mask = np.tile(mask.reshape(1,Nlat,Nlon), (num_scales,1,1))

    tiled_scales = np.tile(scales[:-1].reshape(num_scales,1,1), (1,Nlat,Nlon))
    y_data = y_data * (tiled_scales**2)

    x_data = x_data[tiled_mask == 1]
    y_data = y_data[tiled_mask == 1]

    PlotTools.SignedLogScatter_hist(x_data.ravel(), y_data.ravel(), axes,
            force_equal = True, nbins_x = 300, nbins_y = 300)

    axes[0,0].set_xticklabels([])
    axes[0,1].set_xticklabels([])
    axes[0,1].set_yticklabels([])
    axes[1,1].set_yticklabels([])

    for ax in axes.ravel():
        ax.set_aspect('equal')
        ax.tick_params(axis='both', left=True, right=True, bottom=True, top=True)

    for ax in axes.ravel():
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        ax.plot(xlim, xlim,'--c', label='$1:1$')
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

    axes[0,1].legend(loc='best')

    # xlabel
    mid_x = 0.5 * ( axes[0,0].get_position().x0 + axes[1,1].get_position().x1 )
    plt.figtext(mid_x, 0.025, '$\Pi$ $(\mathrm{W}\cdot\mathrm{km}^{-2}\cdot\mathrm{m}^{-1})$',
             horizontalalignment='center', verticalalignment='top', rotation='horizontal', fontsize=16)

    # ylabel
    mid_y = 0.5 * ( axes[0,0].get_position().y0 + axes[1,1].get_position().y1 )
    plt.figtext(0.025, mid_y, '$l^2\Lambda^m$ $(\mathrm{W}\cdot\mathrm{km}^{-2}\cdot\mathrm{m}^{-1})$',
           horizontalalignment='right', verticalalignment='center', rotation='vertical', fontsize=16)

    plt.savefig('Figures/transfers_comparison_hist.png', dpi=dpi)
    plt.close()

# Lambda_m vs vorticity
if 'baroclinic_transfer' in results.variables:

    # Initialize figure
    fig, axes = plt.subplots(2, 2, 
            gridspec_kw = dict(left = 0.1, right = 0.95, bottom = 0.1, top = 0.95,
                hspace=0.02, wspace=0.02),
            figsize=(7.5, 6) )

    tiled_scales = np.tile(scales[:-1].reshape(num_scales,1,1), (1,Nlat,Nlon))
    tiled_mask   = np.tile(mask.reshape(1,Nlat,Nlon), (num_scales,1,1))
    tiled_rho    = np.tile(rho.reshape(1,Nlat,Nlon), (num_scales,1,1))

    x_data = bc_transfer[:-1,:,:] / tiled_rho
    #x_data = x_data * (tiled_scales**2)
    y_data = np.tile(np.sum(vort_r, axis=0).reshape(1,Nlat,Nlon), (num_scales,1,1))


    x_data = x_data[tiled_mask == 1]
    y_data = y_data[tiled_mask == 1]

    PlotTools.SignedLogScatter_hist(x_data.ravel(), y_data.ravel(), axes,
            force_equal = True, nbins_x = 300, nbins_y = 300)

    axes[0,0].set_xticklabels([])
    axes[0,1].set_xticklabels([])
    axes[0,1].set_yticklabels([])
    axes[1,1].set_yticklabels([])

    for ax in axes.ravel():
        ax.set_aspect('equal')
        ax.tick_params(axis='both', left=True, right=True, bottom=True, top=True)

    for ax in axes.ravel():
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        #ax.plot(xlim, xlim,'--c', label='$1:1$')
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

    #axes[0,1].legend(loc='best')

    # xlabel
    mid_x = 0.5 * ( axes[0,0].get_position().x0 + axes[1,1].get_position().x1 )
    plt.figtext(mid_x, 0.025, '$\Lambda^m/\\rho$ $(\mathrm{s}^{-3})$',
             horizontalalignment='center', verticalalignment='top', rotation='horizontal', fontsize=16)

    # ylabel
    mid_y = 0.5 * ( axes[0,0].get_position().y0 + axes[1,1].get_position().y1 )
    plt.figtext(0.025, mid_y, '$\omega_r$ $(\mathrm{s}^{-1})$',
           horizontalalignment='right', verticalalignment='center', rotation='vertical', fontsize=16)

    plt.savefig('Figures/transfers_comparison2_hist.png', dpi=dpi)
    plt.close()
