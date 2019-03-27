import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import cmocean, sys, PlotTools
from netCDF4 import Dataset
from matplotlib.colors import ListedColormap

dpi = PlotTools.dpi

fp = 'filter_output.nc'
results = Dataset(fp, 'r')
source = Dataset('input.nc', 'r')

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
mask      = results.variables['mask'][:]
u_r       = results.variables['u_r'  ][:, 0, 0, :, :]
u_lat     = results.variables['u_lat'][:, 0, 0, :, :]
u_lon     = results.variables['u_lon'][:, 0, 0, :, :]
# Ordering: scale, time, depth, lat, lon

uo = source.variables['uo'][0,0,:,:]
vo = source.variables['vo'][0,0,:,:]

num_scales = len(scales)

LON, LAT = np.meshgrid(longitude * D2R, latitude * D2R)

# Some parameters for plotting
proj = PlotTools.MapProjection(longitude, latitude)
Xp, Yp = proj(LON*R2D, LAT*R2D, inverse=False)

cbar_props     = dict(pad = 0.1, shrink = 0.85, orientation = 'horizontal')
gridspec_props = dict(wspace = 0.05, hspace = 0.05, left = 0.05, right = 0.95, bottom = 0.05, top = 0.95)

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

        CV  = np.max(np.abs(to_plot))
        qm  = axes[jj,ii].pcolormesh(Xp, Yp, to_plot, 
                cmap='cmo.balance', vmin = -CV, vmax = CV)
    
        cbar = plt.colorbar(qm, ax = axes[jj,ii], **cbar_props)
        PlotTools.ScientificCbar(cbar, units='m/s', orientation='horizontal')

        # Add coastlines and lat/lon lines
        axes[jj,ii].pcolormesh(Xp, Yp, mask, vmin=-1, vmax=1, cmap=mask_cmap)
        PlotTools.AddParallels_and_Meridians(axes[jj,ii], proj, 
            parallels, meridians, latitude, longitude)
        
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
        
    
plt.savefig('Videos/velocity_filter_bands.png', dpi=dpi)
plt.close()

