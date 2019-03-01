import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import cmocean, sys, PlotTools
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap

# The purpose of this post-processing script is to read in the results
#   from a series of filterings and produce images of the
#   band-filtered kinetic energy.
# A major underlying assumption is that the grid is unchanged
#   between the filterings so that subtraction etc.
#   is trivial.

fp = 'filter_output.nc'
results = Dataset(fp, 'r')

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
u_r       = results.variables['u_r'  ][:, 0, 0, :, :]
u_lat     = results.variables['u_lat'][:, 0, 0, :, :]
u_lon     = results.variables['u_lon'][:, 0, 0, :, :]
# Ordering: scale, time, depth, lat, lon

num_scales = len(scales)

LON, LAT = np.meshgrid(longitude * D2R, latitude * D2R)

# Some parameters for plotting
map_settings = PlotTools.MapSettings(longitude, latitude)
cbar_props     = dict(pad = 0.1, shrink = 0.85, orientation = 'horizontal')
gridspec_props = dict(wspace = 0.05, hspace = 0.05, left = 0.02, right = 0.98, bottom = 0.02, top = 0.98)

meridians = np.round(np.linspace(longitude.min(), longitude.max(), 5))
parallels = np.round(np.linspace(latitude.min(),  latitude.max(), 5))

# Get the filter profiles


##
## Begin Plotting
##

## First plot: straight KE binning

# Initialize figure
fig, axes = plt.subplots(3, num_scales,
        sharex=True, sharey=True, 
        gridspec_kw = gridspec_props,
        figsize=(4*num_scales, 14))

# Plot each band
for ii in range(num_scales):
    for jj in range(3):
    
        if jj == 0:
            to_plot = u_r[ii,:,:]
        if jj == 1:
            to_plot = u_lon[ii,:,:]
        if jj == 2:
            to_plot = u_lat[ii,:,:]
        to_plot = np.ma.masked_where(mask==0, to_plot)

        m  = Basemap(ax = axes[jj,ii], **map_settings)

        CV  = np.max(np.abs(to_plot))
        qm  = m.pcolormesh(LON*R2D, LAT*R2D, to_plot, 
                cmap='cmo.balance', 
                vmin = -CV, vmax = CV, latlon = True)
    
        cbar = plt.colorbar(qm, ax = axes[jj,ii], **cbar_props)
        PlotTools.ScientificCbar(cbar, units='m/s', orientation='horizontal')

        # Add coastlines and lat/lon lines
        m.drawcoastlines(linewidth=0.5)

        if ii == num_scales - 1:
            m.drawparallels(parallels, linewidth=0.5, labels=[0,1,0,0], color='g')
        else:
            m.drawparallels(parallels, linewidth=0.5, labels=[0,0,0,0], color='g')
        if jj == 2:
            m.drawmeridians(meridians, linewidth=0.5, labels=[0,0,0,1], color='g')
        else:
            m.drawmeridians(meridians, linewidth=0.5, labels=[0,0,0,0], color='g')
    
        # Draw the mask back on
        m.contourf(LON*R2D, LAT*R2D, mask, [-0.5, 0.5], colors='gray', hatches=['','///\\\\\\'], latlon=True)
        
        if (ii == 0):
            if (jj == 0):
                axes[jj,ii].set_ylabel('$u_r$', fontsize=16)
            if (jj == 1):
                axes[jj,ii].set_ylabel('$u_\lambda$', fontsize=16)
            if (jj == 2):
                axes[jj,ii].set_ylabel('$u_\phi$', fontsize=16)

    
    if (ii == 0):
        axes[0,ii].set_title('Below {0:0.1f} km'.format(scales[0] / 1e3))
    elif (ii == num_scales-1):
        axes[0,ii].set_title('Above {0:0.1f} km'.format(scales[ii-1] / 1e3))
    else:
        axes[0,ii].set_title('{0:.1f} to {1:0.1f} km'.format(scales[ii-1] / 1e3, scales[ii] / 1e3))
        
    
plt.savefig('Figures/velocity_filter_bands.png', dpi=500)
plt.close()

