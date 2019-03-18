import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import PlotTools, cmocean
import matpy as mp

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
v_r       = results.variables['vort_r'  ][:]
v_lon     = results.variables['vort_lon'][:]
v_lat     = results.variables['vort_lat'][:]
Lambda    = results.variables['baroclinic_transfer'][:]

Lambda = Lambda[:-1,:,:,:,:]

Nscales, Ntime, Ndepth, Nlat, Nlon = v_r.shape

dlat = (latitude[1]  - latitude[0] ) * D2R
dlon = (longitude[1] - longitude[0]) * D2R
LON, LAT = np.meshgrid(longitude * D2R, latitude * D2R)
dAreas = R_earth**2 * np.cos(LAT) * dlat * dlon

mask = np.tile(mask.reshape(1,1,Nlat,Nlon), (Ntime, Ndepth, 1, 1))

Dt = mp.FiniteDiff(time, 1, spb=False)

## Scatter of fluxes vs Lambda
for iS in range(Nscales - 1):
    # First, sort out data
    EN_from_vort = 0.5 * (   np.sum(v_r[  iS+1:,:,:,:,:], axis=0)**2 
                           + np.sum(v_lat[iS+1:,:,:,:,:], axis=0)**2 
                           + np.sum(v_lon[iS+1:,:,:,:,:], axis=0)**2 )
    EN_from_vort = EN_from_vort.reshape(Ntime, Ndepth*Nlat*Nlon)
    EN_flux = np.matmul(Dt, EN_from_vort)

    EN_flux    = EN_flux.ravel()[mask.ravel() == 1]
    Lambda_sel = Lambda.ravel()[mask.ravel() == 1]

    # Then plot
    # Initialize figure
    fig, axes = plt.subplots(2, 2, 
            gridspec_kw = dict(left = 0.1, right = 0.95, bottom = 0.1, top = 0.95,
                hspace=0.02, wspace=0.02),
            figsize=(7.5, 6) )
    
    PlotTools.SignedLogScatter_hist(Lambda_sel, EN_flux, axes,
            force_equal = True, num_ords_x = 16, num_ords_y = 16,
            nbins_x = 200, nbins_y = 200)
    
    for II in range(2):
        axes[II,0].set_ylabel('$\\frac{d}{dt}\left(\\frac{1}{2}\overline{\omega}\cdot\overline{\omega}\\right)$ $(\mathrm{W}^{\omega}\cdot\mathrm{m}^{-3})$')
        axes[1,II].set_xlabel('$\Lambda$ $(\mathrm{W}^{\omega}\cdot\mathrm{m}^{-3})$')
    
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
    
    plt.savefig('Figures/EN_fluxes_{0:.3g}km.png'.format(scales[iS]/1e3), dpi=500)
    plt.close()

