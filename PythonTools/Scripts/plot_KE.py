import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import cmocean, sys
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
import PlotTools

# The purpose of this post-processing script is to read in the results
#   from a series of filterings and produce images of the
#   band-filtered kinetic energy.
# A major underlying assumption is that the grid is unchanged
#   between the filterings so that subtraction etc.
#   is trivial.

source = Dataset('input.nc', 'r')

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
time      = results.variables['time'][:]
mask      = results.variables['mask'][:]

num_scales = len(scales)-1

dlat = latitude[1]  - latitude[0]
dlon = longitude[1] - longitude[0]
LON, LAT = np.meshgrid(longitude * D2R, latitude * D2R)
dAreas = R_earth**2 * np.cos(LAT * np.pi / 180) * dlat * dlon

KE = results.variables['KE'][:, 0, 0, :, :]
uo = source.variables['uo'][0, 0, :, :]
vo = source.variables['vo'][0, 0, :, :]
Full_KE = 0.5 * (uo**2 + vo**2)


# Some parameters for plotting
map_settings = PlotTools.MapSettings(longitude, latitude)

meridians = np.round(np.linspace(longitude.min(), longitude.max(), 5))
parallels = np.round(np.linspace(latitude.min(),  latitude.max(),  5))

cbar_props = dict(pad = 0.1, shrink = 0.85, orientation = 'vertical')
gridspec_props = dict(wspace = 0.05, hspace = 0.05, left = 0.02, right = 0.98, bottom = 0.02, top = 0.98)


##
## Begin Plotting
##

## First plot: straight KE binning

# Initialize figure
fig, axes = plt.subplots(num_scales, 1,
        sharex=True, sharey=True, 
        gridspec_kw = gridspec_props,
        figsize=(6,4*num_scales))

# Plot each band
for ii in range(num_scales):
    
    to_plot = KE[ii,:,:]
    to_plot = np.ma.masked_where(mask==0, to_plot)

    m  = Basemap(ax = axes[ii], **map_settings)

    CV  = np.nanmax(np.abs(to_plot))
    if (CV == 0):
        CV = 1
    KE_min = np.min(to_plot)
    if (KE_min < 0) :
        print("min(KE) = {0:.4g} < 0!".format(KE_min), ii)
        qm = m.pcolormesh(LON*R2D, LAT*R2D, to_plot, cmap='cmo.balance', vmin = -CV, vmax = CV, latlon = True)
    else:
        qm = m.pcolormesh(LON*R2D, LAT*R2D, to_plot, cmap='cmo.amp', vmin = 0, vmax = CV, latlon = True)
    
    cbar = plt.colorbar(qm, ax = axes[ii], **cbar_props)
    PlotTools.ScientificCbar(cbar, units='')

    # Add coastlines and lat/lon lines
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(parallels, linewidth=0.5, labels=[0,0,0,0], color='g')
    m.drawmeridians(meridians, linewidth=0.5, labels=[0,0,0,0], color='g')
    m.contourf(LON*R2D, LAT*R2D, mask, [-0.5, 0.5], colors='gray', hatches=['','///\\\\\\'], latlon=True)

    # Also contour the KE for interests sake
    m.contour(LON*R2D, LAT*R2D, Full_KE, 
            levels = np.array([0, 0.025, 0.1, 0.2]) * np.max(Full_KE * mask),
            cmap='cmo.algae', latlon=True, linewidths=0.2)
        
    if (ii == 0):
        axes[ii].set_title('Below {0:0.1f} km'.format(scales[0] / 1e3))
    elif (ii == num_scales-1):
        axes[ii].set_title('Above {0:0.1f} km'.format(scales[ii-1] / 1e3))
    else:
        axes[ii].set_title('{0:.1f} to {1:0.1f} km'.format(scales[ii-1] / 1e3, scales[ii] / 1e3))
        
    
plt.savefig('Figures/KE_band_filter.png', dpi=500)
plt.close()


## Second plot: normalized KE binning

# Initialize figure
fig, axes = plt.subplots(num_scales, 1,
        sharex=True, sharey=True, 
        gridspec_kw = gridspec_props,
        figsize=(6,4*num_scales))

# Plot each band
band_sum = np.sum(KE, axis=0)
for ii in range(num_scales):
    
    to_plot = KE[ii,:,:] / ( band_sum + eps )
    to_plot = np.ma.masked_where(mask==0, to_plot)

    m  = Basemap(ax = axes[ii], **map_settings)

    CV  = np.nanmax(np.abs(to_plot))
    if (CV > 1) :
        print("CV = {0:.4g} > 1!".format(CV), ii)
    qm = m.pcolormesh(LON*R2D, LAT*R2D, to_plot, cmap='cmo.amp', vmin = 0, vmax = 1, latlon = True)
    
    cbar = plt.colorbar(qm, ax = axes[ii], **cbar_props)
    PlotTools.ScientificCbar(cbar, units='')

    # Add coastlines and lat/lon lines
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(parallels, linewidth=0.5, labels=[0,0,0,0], color='g')
    m.drawmeridians(meridians, linewidth=0.5, labels=[0,0,0,0], color='g')
    m.contourf(LON*R2D, LAT*R2D, mask, [-0.5, 0.5], colors='gray', hatches=['','///\\\\\\'], latlon=True)

    # Also contour the KE for interests sake
    m.contour(LON*R2D, LAT*R2D, Full_KE, 
            levels = np.array([0, 0.025, 0.1, 0.2]) * np.max(Full_KE * mask),
            cmap='cmo.algae', latlon=True, linewidths=0.2)
        
    if (ii == 0):
        axes[ii].set_title('Below {0:0.1f} km'.format(scales[0] / 1e3))
    elif (ii == num_scales-1):
        axes[ii].set_title('Above {0:0.1f} km'.format(scales[ii-1] / 1e3))
    else:
        axes[ii].set_title('{0:.1f} to {1:0.1f} km'.format(scales[ii-1] / 1e3, scales[ii] / 1e3))
        
    
plt.savefig('Figures/KE_band_filter_relative.png', dpi=500)
plt.close()


## Third plot: KE distribution

## Now create a histogram of the energy in each bin
Tot_KEs = np.zeros(num_scales)
labels  = []

for ii in range(num_scales):
    Tot_KEs[ii] = np.sum(KE[ii,:,:] * mask * dAreas)
    
    if (ii == 0):
        lab = 'Below {0:0.1f} km'.format(scales[0] / 1e3)
    elif (ii == num_scales-1):
        lab = 'Above {0:0.1f} km'.format(scales[ii-1] / 1e3)
    else:
        lab = '{0:.1f} to {1:0.1f} km'.format(scales[ii-1] / 1e3, scales[ii] / 1e3)
    labels += [lab]

Tot_KE = np.sum(Tot_KEs)

plt.figure(figsize=(5,5))

ax = plt.gca()

ax.bar(np.arange(num_scales), Tot_KEs / Tot_KE, tick_label = labels)
ax.set_xticklabels(ax.get_xticklabels(), ha='right')

ax.set_ylabel('Relative Proportion of KE\nwithin Filter Bands')
ax.tick_params(axis='x', rotation=30)
plt.tight_layout(True)

plt.savefig('Figures/KE_distribution.pdf')
plt.close()
