import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import cmocean, sys, PlotTools
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import ListedColormap

fp = 'filter_output.nc'
results = Dataset(fp, 'r')

R_earth = 6371e3
D2R = np.pi / 180
R2D = 180 / np.pi
eps = 1e-10

# Create cmap for mask data
ref_cmap = cmocean.cm.gray
mask_cmap = ref_cmap(np.arange(ref_cmap.N))
mask_cmap[:,-1] = np.linspace(1, 0, ref_cmap.N)
mask_cmap = ListedColormap(mask_cmap)

# Get the grid from the first filter
latitude  = results.variables['latitude'][:] * R2D
longitude = results.variables['longitude'][:] * R2D
scales    = results.variables['scale'][:]
depth     = results.variables['depth'][:]
time      = results.variables['time'][:]
vort_r    = results.variables['vort_r'  ][:, 0, 0, :, :]
vort_lat  = results.variables['vort_lat'][:, 0, 0, :, :]
vort_lon  = results.variables['vort_lon'][:, 0, 0, :, :]
mask      = results.variables['mask'][:]

num_scales = len(scales)

LON, LAT = np.meshgrid(longitude * D2R, latitude * D2R)

# Some parameters for plotting
map_settings = PlotTools.MapSettings(longitude, latitude)

meridians = np.round(np.linspace(longitude.min(), longitude.max(), 5))
parallels = np.round(np.linspace(latitude.min(),  latitude.max(),  5))

cbar_props     = dict(pad = 0.1, shrink = 0.85, orientation = 'horizontal')
gridspec_props = dict(wspace = 0.05, hspace = 0.05, left = 0.02, right = 0.98, bottom = 0.02, top = 0.98)


##
## Begin Plotting
##

## Vorticity binning

# Initialize figure
fig, axes = plt.subplots(1, num_scales,
        sharex=True, sharey=True, squeeze=False,
        gridspec_kw = gridspec_props,
        figsize=(4*num_scales, 14/3.))

# Plot each band
for ii in range(num_scales):
    for jj in range(1):  # Right now, just plot vort_r
    
        if jj == 0:
            to_plot = vort_r[ii,:,:]
        if jj == 1:
            to_plot = vort_lon[ii,:,:]
        if jj == 2:
            to_plot = vort_lat[ii,:,:]
        to_plot = np.ma.masked_where(mask==0, to_plot)

        m  = Basemap(ax = axes[jj,ii], **map_settings)

        CV  = np.percentile(np.abs(to_plot), 99.9)
        qm  = m.pcolormesh(LON*R2D, LAT*R2D, to_plot, 
                cmap='cmo.balance', 
                vmin = -CV, vmax = CV, latlon = True)
    
        cbar = plt.colorbar(qm, ax = axes[jj,ii], **cbar_props)
        PlotTools.ScientificCbar(cbar, units='m/s', orientation='horizontal')

        # Add coastlines and lat/lon lines
        m.drawcoastlines(linewidth=0.1)

        if ii == num_scales - 1:
            m.drawparallels(parallels, linewidth=0.5, labels=[0,1,0,0], color='g')
        else:
            m.drawparallels(parallels, linewidth=0.5, labels=[0,0,0,0], color='g')
        if jj == 2:
            m.drawmeridians(meridians, linewidth=0.5, labels=[0,0,0,1], color='g')
        else:
            m.drawmeridians(meridians, linewidth=0.5, labels=[0,0,0,0], color='g')
    
        # Draw the mask back on
        m.pcolormesh(LON*R2D, LAT*R2D, mask, vmin=-1, vmax=1, cmap=mask_cmap, latlon=True)

        if (ii == 0):
            if (jj == 0):
                axes[jj,ii].set_ylabel('$\Omega_r$', fontsize=16)
            if (jj == 1):
                axes[jj,ii].set_ylabel('$\Omega_\lambda$', fontsize=16)
            if (jj == 2):
                axes[jj,ii].set_ylabel('$\Omega_\phi$', fontsize=16)

    
    if (ii == 0):
        axes[0,ii].set_title('Below {0:0.1f} km'.format(scales[0] / 1e3))
    elif (ii == num_scales-1):
        axes[0,ii].set_title('Above {0:0.1f} km'.format(scales[ii-1] / 1e3))
    else:
        axes[0,ii].set_title('{0:.1f} to {1:0.1f} km'.format(scales[ii-1] / 1e3, scales[ii] / 1e3))
        
    
plt.savefig('Figures/vorticity_filter_bands.png', dpi=500)
plt.close()


## Dichotomies

# Initialize figure
fig, axes = plt.subplots(num_scales-1, 2,
        sharex=True, sharey=True, squeeze=False,
        gridspec_kw = gridspec_props,
        figsize=(10, 4*num_scales-1))

# Plot each band
for ii in range(num_scales-1):
    
    to_plot_below = np.sum(vort_r[:ii+1, :,:], axis=0)
    to_plot_above = np.sum(vort_r[ ii+1:,:,:], axis=0)

    to_plot_below = np.ma.masked_where(mask==0, to_plot_below)
    to_plot_above = np.ma.masked_where(mask==0, to_plot_above)

    m_a = Basemap(ax = axes[ii,1], **map_settings)
    m_b = Basemap(ax = axes[ii,0], **map_settings)

    CV_a = np.nanmax(np.abs(to_plot_above))
    CV_b = np.nanmax(np.abs(to_plot_below))

    vmax = max(CV_a, CV_b)
    vmin = 10**(np.log10(vmax) - 3)

    PlotTools.SignedLogPlot_onMap(LON * R2D, LAT * R2D, to_plot_above, 
            axes[ii,1], fig, m_a, num_ords = 2, vmax=vmax)
    PlotTools.SignedLogPlot_onMap(LON * R2D, LAT * R2D, to_plot_below, 
            axes[ii,0], fig, m_b, num_ords = 2, vmax=vmax)

    # Add coastlines and lat/lon lines
    for m in [m_a, m_b]:
        m.drawcoastlines(linewidth=0.1)
        m.drawparallels(parallels, linewidth=0.5, labels=[0,0,0,0], color='g')
        m.drawmeridians(meridians, linewidth=0.5, labels=[0,0,0,0], color='g')
        m.pcolormesh(LON*R2D, LAT*R2D, mask, vmin=-1, vmax=1, cmap=mask_cmap, latlon=True)
    
    axes[ii,0].set_ylabel('{0:.1f} km'.format(scales[ii] / 1e3))

axes[0,0].set_title('Fine $(<l)$')
axes[0,1].set_title('Coarse $(>l)$')
        
    
plt.savefig('Figures/vorticity_dichotomies.png', dpi=500)
plt.close()

