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

uo = source.variables['uo'][0, 0, :, :]
vo = source.variables['vo'][0, 0, :, :]
Full_KE = 0.5 * (uo**2 + vo**2)

num_scales = len(scales)-1

dlat = latitude[1] - latitude[0]
dlon = longitude[1] - longitude[0]
LON, LAT = np.meshgrid(longitude * D2R, latitude * D2R)
dAreas = R_earth**2 * np.cos(LAT * np.pi / 180) * dlat * dlon

map_settings = PlotTools.MapSettings(longitude, latitude)

meridians = np.round(np.linspace(longitude.min(), longitude.max(), 5))
parallels = np.round(np.linspace(latitude.min(),  latitude.max(),  5))

cbar_props     = dict(pad = 0.1, shrink = 0.85, orientation = 'vertical')
gridspec_props = dict(wspace = 0.05, hspace = 0.05, left = 0.02, right = 0.98, bottom = 0.02, top = 0.98)

##
## Begin Plotting
##

## First plot: straight KE binning

# Initialize figure
fig, axes = plt.subplots(num_scales, 1,
        sharex=True, sharey=True, 
        gridspec_kw = gridspec_props,
        figsize=(6, 4*num_scales))

# Plot each band
for ii in range(num_scales):
    
    to_plot = transfer[ii,:,:] * 1e6
    to_plot = np.ma.masked_where(mask==0, to_plot)

    m  = Basemap(ax = axes[ii], **map_settings)

    PlotTools.SignedLogPlot_onMap(LON * R2D, LAT * R2D, to_plot, axes[ii], fig, m, num_ords = 3)

    # Add coastlines, lat/lon lines, and draw the map
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(parallels, linewidth=0.5, labels=[0,0,0,0], color='k')
    m.drawmeridians(meridians, linewidth=0.5, labels=[0,0,0,0], color='k')
    m.contourf(LON*R2D, LAT*R2D, mask, [-0.5, 0.5], colors='gray', hatches=['','///\\\\\\'], latlon=True)

    # Also contour the KE for interests sake
    m.contour(LON*R2D, LAT*R2D, Full_KE, 
            levels = np.array([0, 0.025, 0.1, 0.2]) * np.max(Full_KE * mask),
            cmap='cmo.algae', latlon=True, linewidths=0.2)
        
    axes[ii].set_title('Through {0:0.1f} km'.format(scales[ii] / 1e3))
        
    
plt.savefig('Figures/energy_transfers.png', dpi=500)
plt.close()

