import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import PlotTools, cmocean
import matpy as mp
import sys

fp = 'filter_output.nc'
results = Dataset(fp, 'r')

R_earth = 6371e3
D2R     = np.pi / 180
R2D     = 180 / np.pi
eps     = 1e-10

# Get the grid from the first filter
latitude  = results.variables['latitude'][:] * R2D
longitude = results.variables['longitude'][:] * R2D
scales    = results.variables['scale'][:]
depth     = results.variables['depth'][:]
time      = results.variables['time'][:] / (60 * 60) # time was in hours
mask      = results.variables['mask'][:]
u_r       = results.variables['u_r'  ][:]
u_lon     = results.variables['u_lon'][:]
u_lat     = results.variables['u_lat'][:]
Pi        = results.variables['energy_transfer'][:]

Pi = Pi[:-1,:,:,:,:]

Nscales, Ntime, Ndepth, Nlat, Nlon = u_r.shape

dlat = (latitude[1]  - latitude[0] ) * D2R
dlon = (longitude[1] - longitude[0]) * D2R
LON, LAT = np.meshgrid(longitude * D2R, latitude * D2R)
dAreas = R_earth**2 * np.cos(LAT) * dlat * dlon

dArea = np.tile((mask*dAreas).reshape(1,Nlat,Nlon), (Ndepth,1,1)) 
mask  = np.tile(mask.reshape(1,1,Nlat,Nlon), (Ntime, Ndepth, 1, 1))

Dt = mp.FiniteDiff(time, 1, spb=False)

for iS in range(Nscales-1):
    # First, process data
    KE_from_vel = 0.5 * (  np.sum(u_r[  iS+1:,:,:,:,:], axis=0)**2 
                         + np.sum(u_lat[iS+1:,:,:,:,:], axis=0)**2 
                         + np.sum(u_lon[iS+1:,:,:,:,:], axis=0)**2 )
    KE_from_vel = KE_from_vel.reshape(Ntime, Ndepth*Nlat*Nlon)
    KE_flux = np.matmul(Dt, KE_from_vel)

    KE_flux =  KE_flux.ravel()[mask.ravel() == 1]
    Pi_sel  = -Pi.ravel()[mask.ravel() == 1]

    # Then plot
    # Initialize figure
    fig, axes = plt.subplots(2, 2, 
            gridspec_kw = dict(left = 0.1, right = 0.95, bottom = 0.1, top = 0.95,
                hspace=0.02, wspace=0.02),
            figsize=(7.5, 6) )
    
    PlotTools.SignedLogScatter_hist(Pi_sel, KE_flux, axes,
            force_equal = True, num_ords_x = 12, num_ords_y = 12,
            nbins_x = 200, nbins_y = 200)
    
    for II in range(2):
        axes[II,0].set_ylabel('$\\frac{d}{dt}\left( \\frac{1}{2}\overline{u}\cdot\overline{u} \\right)$ $(\mathrm{W}\cdot\mathrm{m}^{-3})$')
        axes[1,II].set_xlabel('$\Pi$ $(\mathrm{W}\cdot\mathrm{m}^{-3})$')
    
    for ax in axes.ravel():
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        ax.plot(xlim, xlim,'--c', label='$1:1$')
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
    
    axes[0,1].legend(loc='best')
    
    axes[0,0].set_xticklabels([])
    axes[0,1].set_xticklabels([])
    axes[0,1].set_yticklabels([])
    axes[1,1].set_yticklabels([])
    
    plt.savefig('Figures/KE_fluxes_{0:.3g}km.png'.format(scales[iS]/1e3), dpi=500)
    plt.close()

