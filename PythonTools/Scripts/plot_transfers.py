import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import cmocean, sys
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
import PlotTools, shutil, os, datetime
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
time      = results.variables['time'][:] * (60*60) # hours to seconds
mask      = results.variables['mask'][:]
transfer  = results.variables['energy_transfer'][:, :, 0, :, :]
if 'baroclinic_transfer' in results.variables:
    bc_transfer = results.variables['baroclinic_transfer'][:, :, 0, :, :]

# Do some time handling tp adjust the epochs
# appropriately
epoch = datetime.datetime(1950,1,1)   # the epoch of the time dimension
dt_epoch = datetime.datetime.fromtimestamp(0)  # the epoch used by datetime
epoch_delta = dt_epoch - epoch  # difference
time = time - epoch_delta.total_seconds()  # shift

Ntime = len(time)
num_scales = len(scales)-1

LON, LAT = np.meshgrid(longitude * D2R, latitude * D2R)

Nlat = len(latitude)
Nlon = len(longitude)

map_settings = PlotTools.MapSettings(longitude, latitude)
plt.figure()
ax = plt.gca()
proj = Basemap(ax = ax, **map_settings)
Xp, Yp = proj(LON*R2D, LAT*R2D)

meridians = np.round(np.linspace(longitude.min(), longitude.max(), 5))
parallels = np.round(np.linspace(latitude.min(),  latitude.max(),  5))

cbar_props     = dict(pad = 0.1, shrink = 0.85, orientation = 'vertical')
gridspec_props = dict(wspace = 0.05, hspace = 0.05, left = 0.1, right = 0.9, bottom = 0.1, top = 0.9)

##
## Begin Plotting
##

for Itime in range(rank, Ntime, num_procs):    

    timestamp = datetime.datetime.fromtimestamp(time[Itime])
    sup_title = "{0:02d} - {1:02d} - {2:04d} ( {3:02d}:{4:02d} )".format(
            timestamp.day, timestamp.month, timestamp.year, 
            timestamp.hour, timestamp.minute)
    
    # Plot each band
    for ii in range(num_scales):
    
        # Initialize figure
        fig, axes = plt.subplots(1, 1,
            sharex=True, sharey=True, squeeze=False,
            gridspec_kw = gridspec_props,
            figsize=(6, 4))

        fig.suptitle(sup_title)
        
        to_plot = transfer[ii,Itime,:,:] * 1e6
        to_plot = np.ma.masked_where(mask==0, to_plot)
    
        if np.max(np.abs(to_plot)) > 0:
            PlotTools.SignedLogPlot_onMap(LON * R2D, LAT * R2D, to_plot, axes[0,0], fig, proj, num_ords = 5)
    
        # Add coastlines, lat/lon lines, and draw the map
        axes[0,0].pcolormesh(Xp, Yp, mask, vmin=-1, vmax=1, cmap=mask_cmap)
        PlotTools.AddParallels_and_Meridians(axes[0,0], proj, 
            parallels, meridians, latitude, longitude)
            
        plt.savefig(tmp_direct + '/{0:.4g}_Pi_{1:04d}.png'.format(scales[ii]/1e3,Itime), dpi=dpi)
        plt.close()
    
    
# If baroclinic transfers are there, use them.
if 'baroclinic_transfer' in results.variables:
    for Itime in range(rank, Ntime, num_procs):    

        timestamp = datetime.datetime.fromtimestamp(time[Itime])
        sup_title = "{0:02d} - {1:02d} - {2:04d} ( {3:02d}:{4:02d} )".format(
            timestamp.day, timestamp.month, timestamp.year, 
            timestamp.hour, timestamp.minute)
    
        # Plot each band
        for ii in range(num_scales):
    
            # Initialize figure
            fig, axes = plt.subplots(1, 1,
                sharex=True, sharey=True, squeeze=False,
                gridspec_kw = gridspec_props,
                figsize=(6, 4))

            fig.suptitle(sup_title)
        
            to_plot = bc_transfer[ii,Itime,:,:] * 1e6
            to_plot = np.ma.masked_where(mask==0, to_plot)
    
            if np.max(np.abs(to_plot)) > 0:
                PlotTools.SignedLogPlot_onMap(LON * R2D, LAT * R2D, to_plot, axes[0,0], fig, proj, num_ords = 4, percentile=99.99)
    
            # Add coastlines, lat/lon lines, and draw the map
            Xp, Yp = proj(LON*R2D, LAT*R2D)
            axes[0,0].pcolormesh(Xp, Yp, mask, vmin=-1, vmax=1, cmap=mask_cmap)
            PlotTools.AddParallels_and_Meridians(axes[0,0], proj, 
                parallels, meridians, latitude, longitude)
            
            plt.savefig(tmp_direct + '/{0:.4g}_Lambda_m_{1:04d}.png'.format(scales[ii]/1e3,Itime), dpi=dpi)
            plt.close()


# If more than one time point, create mp4s
for ii in range(num_scales):
    if Ntime > 1:
        PlotTools.merge_to_mp4(tmp_direct + '/{0:.4g}_Pi_%04d.png'.format(scales[ii]/1e3),    
                out_direct + '/{0:.4g}km/Pi.mp4'.format(scales[ii]/1e3),    fps=12)
        if 'baroclinic_transfer' in results.variables:
            PlotTools.merge_to_mp4(tmp_direct + '/{0:.4g}_Lambda_m_%04d.png'.format(scales[ii]/1e3), 
                    out_direct + '/{0:.4g}km/Lambda_m.mp4'.format(scales[ii]/1e3), fps=12)
    else:
        shutil.move(tmp_direct + '/{0:.4g}_Pi_0000.png'.format(scales[ii]/1e3), 
                out_direct + '/{0:.4g}km/Pi.png'.format(scales[ii]/1e3))
        if 'p' in source.variables:
            shutil.move(tmp_direct + '/{0:.4g}_Lambda_m_0000.png'.format(scales[ii]/1e3), 
                    out_direct + '/{0:.4g}km/Lambda_m.png'.format(scales[ii]/1e3))

# Now delete the frames
shutil.rmtree(tmp_direct)

## Plot time averages if Ntime > 1

if (Ntime > 1):
    
    # Plot each band
    for ii in range(num_scales):

        # Initialize figure
        fig, axes = plt.subplots(1, 1,
            sharex=True, sharey=True, squeeze=False,
            gridspec_kw = gridspec_props,
            figsize=(6, 4))

        fig.suptitle('Time average')
        
        to_plot = np.mean(transfer[ii,:,:,:] * 1e6, axis=0)
        to_plot = np.ma.masked_where(mask==0, to_plot)
    
        if np.max(np.abs(to_plot)) > 0:
            PlotTools.SignedLogPlot_onMap(LON * R2D, LAT * R2D, to_plot, axes[0,0], fig, proj, num_ords = 5)
    
        # Add coastlines, lat/lon lines, and draw the map
        axes[0,0].pcolormesh(Xp, Yp, mask, vmin=-1, vmax=1, cmap=mask_cmap)
        PlotTools.AddParallels_and_Meridians(axes[0,0], proj, 
            parallels, meridians, latitude, longitude)
            
        plt.savefig(out_direct + '/{0:.4g}km/AVE_Pi.png'.format(scales[ii]/1e3), dpi=dpi)
        plt.close()
    
    
if 'baroclinic_transfer' in results.variables:
    if (Ntime > 1):
    
        # Plot each band
        for ii in range(num_scales):

            # Initialize figure
            fig, axes = plt.subplots(1, 1,
                sharex=True, sharey=True, squeeze=False,
                gridspec_kw = gridspec_props,
                figsize=(6, 4))

            fig.suptitle('Time average')
        
            to_plot = np.mean(bc_transfer[ii,:,:,:] * 1e6, axis=0)
            to_plot = np.ma.masked_where(mask==0, to_plot)
    
            if np.max(np.abs(to_plot)) > 0:
                PlotTools.SignedLogPlot_onMap(LON * R2D, LAT * R2D, to_plot, axes[0,0], fig, proj, num_ords = 4, percentile=99.99)
    
            # Add coastlines, lat/lon lines, and draw the map
            axes[0,0].pcolormesh(Xp, Yp, mask, vmin=-1, vmax=1, cmap=mask_cmap)
            PlotTools.AddParallels_and_Meridians(axes[0,0], proj, 
                parallels, meridians, latitude, longitude)
        
            plt.savefig(out_direct + '/{0:.4g}km/AVE_Lambda_m.png'.format(scales[ii]/1e3), dpi=dpi)
            plt.close()
