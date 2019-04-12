import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import cmocean, sys, PlotTools, os, shutil, datetime
from netCDF4 import Dataset
from matplotlib.colors import ListedColormap

dpi = PlotTools.dpi

try: # Try using mpi
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    num_procs = comm.Get_size()
except:
    rank = 0
    num_procs = 1
print("Proc {0:d} of {1:d}".format(rank+1,num_procs))

# If the Figures directory doesn't exist, create it.
# Same with the Figures/tmp
out_direct = os.getcwd() + '/Videos'
tmp_direct = out_direct + '/tmp'

if (rank == 0):
    print("Saving outputs to " + out_direct)
    print("  will use temporary directory " + tmp_direct)

    if not(os.path.exists(out_direct)):
        os.makedirs(out_direct)

    if not(os.path.exists(tmp_direct)):
        os.makedirs(tmp_direct)

fp = 'filter_output.nc'
results = Dataset(fp, 'r')
source  = Dataset('input.nc', 'r')

try:
    units = source.variables['latitude'].units
except:
    units = ''

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
if units == 'm':
    latitude  = results.variables['latitude'][:] / 1e3
    longitude = results.variables['longitude'][:] / 1e3
else:
    latitude  = results.variables['latitude'][:] * R2D
    longitude = results.variables['longitude'][:] * R2D
LON, LAT = np.meshgrid(longitude * D2R, latitude * D2R)
scales    = results.variables['scale'][:]
depth     = results.variables['depth'][:]
time      = results.variables['time'][:] * (60*60) # convert hours to second
div       = results.variables['vel_div'][:, :, 0, :, :]
mask      = results.variables['mask'][:]

# Do some time handling tp adjust the epochs
# appropriately
epoch = datetime.datetime(1950,1,1)   # the epoch of the time dimension
dt_epoch = datetime.datetime.fromtimestamp(0)  # the epoch used by datetime
epoch_delta = dt_epoch - epoch  # difference
time = time - epoch_delta.total_seconds()  # shift

num_scales = len(scales)
Ntime = len(time)

# Some parameters for plotting
proj = PlotTools.MapProjection(longitude, latitude)
Xp, Yp = proj(LON * R2D, LAT * R2D, inverse=False)

meridians = np.round(np.linspace(longitude.min(), longitude.max(), 5))
parallels = np.round(np.linspace(latitude.min(),  latitude.max(),  5))

cbar_props     = dict(pad = 0.1, shrink = 0.85)
gridspec_props = dict(wspace = 0.05, hspace = 0.05, left = 0.1, right = 0.9, bottom = 0.1, top = 0.9)


##
## Begin Plotting
##

## Dichotomies
for Itime in range(rank, Ntime, num_procs):    

    timestamp = datetime.datetime.fromtimestamp(time[Itime])
    sup_title = "{0:02d} - {1:02d} - {2:04d} ( {3:02d}:{4:02d} )".format(
            timestamp.day, timestamp.month, timestamp.year, 
            timestamp.hour, timestamp.minute)
    
    # Plot each band
    for ii in range(num_scales-1):

        # Initialize figure
        fig, axes = plt.subplots(2, 1,
            sharex=True, sharey=True, squeeze=False,
            gridspec_kw = gridspec_props,
            figsize=(10, 10))

        fig.suptitle(sup_title)
        
        to_plot_below = np.sum(div[:ii+1 , Itime, :, :], axis=0)
        to_plot_above = np.sum(div[ ii+1:, Itime, :, :], axis=0)
    
        to_plot_below = np.ma.masked_where(mask==0, to_plot_below)
        to_plot_above = np.ma.masked_where(mask==0, to_plot_above)
    
        CV_a = np.nanpercentile(np.abs(to_plot_above), 99)
        CV_b = np.nanpercentile(np.abs(to_plot_below), 99)
    
        qm_a  = axes[0,0].pcolormesh(Xp, Yp, to_plot_above, 
                    cmap='cmo.balance', vmin = -CV_a, vmax = CV_a)
        qm_b  = axes[1,0].pcolormesh(Xp, Yp, to_plot_below, 
                    cmap='cmo.balance', vmin = -CV_b, vmax = CV_b)

        cbar_a = plt.colorbar(qm_a, ax = axes[0,0], **cbar_props)
        cbar_b = plt.colorbar(qm_b, ax = axes[1,0], **cbar_props)
        PlotTools.ScientificCbar(cbar_a, units='')
        PlotTools.ScientificCbar(cbar_b, units='')
    
        # Add coastlines and lat/lon lines
        for ax in axes[:,0]:
            ax.pcolormesh(Xp, Yp, mask, vmin=-1, vmax=1, cmap=mask_cmap)
            PlotTools.AddParallels_and_Meridians(ax, proj, 
                parallels, meridians, latitude, longitude)
        
        axes[0,0].set_ylabel('Coarse $(>l)$')
        axes[1,0].set_ylabel('Fine $(<l)$')
        
        plt.savefig(tmp_direct + '/{0:.4g}_div_vel_dichotomies_{1:04d}.png'.format(scales[ii]/1e3, Itime), dpi=dpi)
        plt.close()

    
if (rank > 0):
    sys.exit()

# If more than one time point, create mp4s
if Ntime > 1:
    for ii in range(num_scales-1):
        PlotTools.merge_to_mp4(tmp_direct + '/{0:.04g}_div_vel_dichotomies_%04d.png'.format(scales[ii]/1e3),    
                out_direct + '/{0:.04g}km/div_vel_dichotomies.mp4'.format(scales[ii]/1e3), fps=12)
    
else:
    shutil.move(tmp_direct + '/vorticity_bands_0000.png',
            out_direct + '/vorticity_bands.png')
    for ii in range(num_scales-1):
        shutil.move(tmp_direct + '/{0:.04g}_vorticity_dichotomies_0000.png'.format(scales[ii]/1e3),
                out_direct + '/{0:.04g}km/vorticity_dichotomies.png'.format(scales[ii]/1e3))

# Now delete the frames
shutil.rmtree(tmp_direct)
