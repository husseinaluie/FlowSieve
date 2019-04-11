import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import cmocean, sys, PlotTools, os, shutil, datetime, glob
from netCDF4 import Dataset
from matplotlib.colors import LogNorm
from matplotlib.colors import ListedColormap

dpi = PlotTools.dpi

# List of variables to (try to) plot
variables = ['energy_transfer', 'Lambda_m', 'PEtoKE', 'div_Jtransport']

try: # Try using mpi
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    num_procs = comm.Get_size()
except:
    rank = 0
    num_procs = 1
print("Proc {0:d} of {1:d}".format(rank+1,num_procs))

# Get the available filter files
files = glob.glob('filter_*.nc')

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

source = Dataset('input.nc', 'r')
try:
    units = source.variables['latitude'].units
except:
    units = ''

# Create cmap for mask data
ref_cmap = cmocean.cm.gray
mask_cmap = ref_cmap(np.arange(ref_cmap.N))
mask_cmap[:,-1] = np.linspace(1, 0, ref_cmap.N)
mask_cmap = ListedColormap(mask_cmap)

# Some parameters for plotting
cbar_props     = dict(pad = 0.02, shrink = 0.85)
gridspec_props = dict(wspace = 0.15, hspace = 0.15, left = 0.1, right = 0.95, bottom = 0.1, top = 0.9)

# Loop through filters
for fp in files:
    with Dataset(fp, 'r') as results:
        scale = results.filter_scale

        # Get the grid from the first filter
        if units == 'm':
            latitude  = results.variables['latitude'][:] / 1e3
            longitude = results.variables['longitude'][:] / 1e3
        else:
            latitude  = results.variables['latitude'][:]
            longitude = results.variables['longitude'][:]
        LON, LAT = np.meshgrid(longitude, latitude)

        depth = results.variables['depth'][:]
        time  = results.variables['time'][:] * (60*60) # convert hours to seconds
        mask  = results.variables['mask'][:]

        Ntime = len(time)

        # Do some time handling tp adjust the epochs
        # appropriately
        epoch       = datetime.datetime(1950,1,1)   # the epoch of the time dimension
        dt_epoch    = datetime.datetime.fromtimestamp(0)  # the epoch used by datetime
        epoch_delta = dt_epoch - epoch  # difference
        time        = time - epoch_delta.total_seconds()  # shift

        # lat/lon lines to draw
        meridians = np.round(np.linspace(longitude.min(), longitude.max(), 5))
        parallels = np.round(np.linspace(latitude.min(),  latitude.max(), 5))

        # Map projection
        proj = PlotTools.MapProjection(longitude, latitude)
        Xp, Yp = proj(LON, LAT, inverse=False)

        rat   = (Xp.max() - Xp.min()) / (Yp.max() - Yp.min())
        rat  *= 0.5 * (1.2)
        fig_h = 8.

        # Plot stuff
        for Itime in range(rank, Ntime, num_procs):    

            timestamp = datetime.datetime.fromtimestamp(time[Itime])
            sup_title = "{0:02d} - {1:02d} - {2:04d} ( {3:02d}:{4:02d} )".format(
                timestamp.day, timestamp.month, timestamp.year, 
                timestamp.hour, timestamp.minute)

            u_r   = results.variables['fine_u_r'  ][Itime,0,:,:]
            u_lon = results.variables['fine_u_lon'][Itime,0,:,:]
            u_lat = results.variables['fine_u_lat'][Itime,0,:,:]
            to_plot_below = 0.5 * ( u_r**2 + u_lon**2 + u_lat**2)

            u_r   = results.variables['coarse_u_r'  ][Itime,0,:,:]
            u_lon = results.variables['coarse_u_lon'][Itime,0,:,:]
            u_lat = results.variables['coarse_u_lat'][Itime,0,:,:]
            to_plot_above = 0.5 * ( u_r**2 + u_lon**2 + u_lat**2)

            # Initialize figure
            fig, axes = plt.subplots(2, 1,
                sharex=True, sharey=True, squeeze=False, 
                gridspec_kw = gridspec_props,
                figsize=(fig_h*rat, fig_h))

            fig.suptitle(sup_title)

            CV_a = np.nanmax(np.abs(to_plot_above))
            CV_b = np.nanmax(np.abs(to_plot_below))
    
            vmax = max(CV_a, CV_b)
            vmin = 10**(np.log10(vmax) - 3)

            qm_a = axes[0,0].pcolormesh(Xp, Yp, to_plot_above, cmap='cmo.amp',
                norm=LogNorm(vmin=vmin, vmax=vmax))
    
            qm_b = axes[1,0].pcolormesh(Xp, Yp, to_plot_below, cmap='cmo.amp',
                norm=LogNorm(vmin=vmin, vmax=vmax))

            cbar = plt.colorbar(qm_b, ax = axes[:,0], **cbar_props)
    
            # Add land and lat/lon lines
            for ax in axes[:,0]:
                ax.pcolormesh(Xp, Yp, mask, vmin=-1, vmax=1, cmap=mask_cmap)
                PlotTools.AddParallels_and_Meridians(ax, proj, 
                    parallels, meridians, latitude, longitude)

            for ax in axes.ravel():
                ax.set_aspect('equal')


            plt.savefig(tmp_direct + '/{0:.4g}_KE_dichotomies_from_vels_{1:04d}.png'.format(
                scale/1e3, Itime), dpi=dpi)
            plt.close()

    
if (rank > 0):
    sys.exit()

# If more than one time point, create mp4s
if Ntime > 1:
    for fp in files:
        with Dataset(fp, 'r') as results:
            scale = results.filter_scale
            PlotTools.merge_to_mp4(tmp_direct + '/{0:.04g}_KE_dichotomies_from_vels_%04d.png'.format(scale/1e3),
                    out_direct + '/{0:.04g}km/KE_dichotomies_from_vels.mp4'.format(scale/1e3), fps=12)
else:
    for fp in files:
        with Dataset(fp, 'r') as results:
            scale = results.filter_scale
            shutil.move(tmp_direct + '/{0:.04g}_KE_dichotomies_from_vels_0000.png'.format(scale/1e3),
                    out_direct + '/{0:.04g}km/KE_dichotomies_from_vels.png'.format(scale/1e3))

# Now delete the frames
shutil.rmtree(tmp_direct)

# Plot time-mean
for fp in files:
    with Dataset(fp, 'r') as results:

        scale = results.filter_scale

        # Get the grid from the first filter
        if units == 'm':
            latitude  = results.variables['latitude'][:] / 1e3
            longitude = results.variables['longitude'][:] / 1e3
        else:
            latitude  = results.variables['latitude'][:]
            longitude = results.variables['longitude'][:]
        LON, LAT = np.meshgrid(longitude, latitude)

        depth = results.variables['depth'][:]
        time  = results.variables['time'][:] * (60*60) # convert hours to seconds
        mask  = results.variables['mask'][:]

        Ntime = len(time)

        # Do some time handling tp adjust the epochs
        # appropriately
        epoch       = datetime.datetime(1950,1,1)   # the epoch of the time dimension
        dt_epoch    = datetime.datetime.fromtimestamp(0)  # the epoch used by datetime
        epoch_delta = dt_epoch - epoch  # difference
        time        = time - epoch_delta.total_seconds()  # shift

        # lat/lon lines to draw
        meridians = np.round(np.linspace(longitude.min(), longitude.max(), 5))
        parallels = np.round(np.linspace(latitude.min(),  latitude.max(), 5))

        # Map projection
        proj = PlotTools.MapProjection(longitude, latitude)
        Xp, Yp = proj(LON, LAT, inverse=False)

        rat   = (Xp.max() - Xp.min()) / (Yp.max() - Yp.min())
        rat  *= 0.5 * (1.2)
        fig_h = 8.

        ## Vorticity dichotomies
        if (Ntime > 1):

            # Initialize figure
            fig, axes = plt.subplots(2, 1,
                sharex=True, sharey=True, squeeze=False, 
                gridspec_kw = gridspec_props,
                figsize=(fig_h*rat, fig_h))

            fig.suptitle('Time average')

            u_r   = results.variables['fine_u_r'  ][:,0,:,:]
            u_lon = results.variables['fine_u_lon'][:,0,:,:]
            u_lat = results.variables['fine_u_lat'][:,0,:,:]
            to_plot_below = 0.5 * ( u_r**2 + u_lon**2 + u_lat**2)
            to_plot_below = np.mean(to_plot_below, axis=0)

            u_r   = results.variables['coarse_u_r'  ][:,0,:,:]
            u_lon = results.variables['coarse_u_lon'][:,0,:,:]
            u_lat = results.variables['coarse_u_lat'][:,0,:,:]
            to_plot_above = 0.5 * ( u_r**2 + u_lon**2 + u_lat**2)
            to_plot_above = np.mean(to_plot_above, axis=0)

            CV_a = np.nanmax(np.abs(to_plot_above))
            CV_b = np.nanmax(np.abs(to_plot_below))
    
            vmax = max(CV_a, CV_b)
            vmin = 10**(np.log10(vmax) - 3)

            qm_a = axes[0,0].pcolormesh(Xp, Yp, to_plot_above, cmap='cmo.amp',
                norm=LogNorm(vmin=vmin, vmax=vmax))
    
            qm_b = axes[1,0].pcolormesh(Xp, Yp, to_plot_below, cmap='cmo.amp',
                norm=LogNorm(vmin=vmin, vmax=vmax))

            cbar = plt.colorbar(qm_b, ax = axes[:,0], **cbar_props)

            # Add land and lat/lon lines
            for ax in axes[:,0]:
                ax.pcolormesh(Xp, Yp, mask, vmin=-1, vmax=1, cmap=mask_cmap)
                PlotTools.AddParallels_and_Meridians(ax, proj, 
                    parallels, meridians, latitude, longitude)

            for ax in axes.ravel():
                ax.set_aspect('equal')
        
            plt.savefig(out_direct + '/{0:.4g}km/AVE_KE_dichotomies_from_vels.png'.format(
                scale/1e3), dpi=dpi)
            plt.close()
