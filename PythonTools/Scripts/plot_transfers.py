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

meridians = np.round(np.linspace(longitude.min(), longitude.max(), 5))
parallels = np.round(np.linspace(latitude.min(),  latitude.max(),  5))

cbar_props     = dict(pad = 0.1, shrink = 0.85, orientation = 'vertical')
gridspec_props = dict(wspace = 0.05, hspace = 0.05, left = 0.05, right = 0.95, bottom = 0.05, top = 0.95)

##
## Begin Plotting
##

for Itime in range(rank, Ntime, num_procs):    

    timestamp = datetime.datetime.fromtimestamp(time[Itime])
    sup_title = "{0:02d} - {1:02d} - {2:04d} ( {3:02d}:{4:02d} )".format(
            timestamp.day, timestamp.month, timestamp.year, 
            timestamp.hour, timestamp.minute)
    
    # Initialize figure
    fig, axes = plt.subplots(num_scales, 1,
            sharex=True, sharey=True, squeeze=False,
            gridspec_kw = gridspec_props,
            figsize=(6, 4*num_scales))

    fig.suptitle(sup_title)
    
    # Plot each band
    for ii in range(num_scales):
        
        to_plot = transfer[ii,Itime,:,:] * 1e6
        to_plot = np.ma.masked_where(mask==0, to_plot)
    
        m  = Basemap(ax = axes[ii,0], **map_settings)
    
        if np.max(np.abs(to_plot)) > 0:
            PlotTools.SignedLogPlot_onMap(LON * R2D, LAT * R2D, to_plot, axes[ii,0], fig, m, num_ords = 5)
    
        # Add coastlines, lat/lon lines, and draw the map
        m.drawcoastlines(linewidth=0.1)
        m.drawparallels(parallels, linewidth=0.5, labels=[0,0,0,0], color='k')
        m.drawmeridians(meridians, linewidth=0.5, labels=[0,0,0,0], color='k')
        m.pcolormesh(LON*R2D, LAT*R2D, mask, vmin=-1, vmax=1, cmap=mask_cmap, latlon=True)
            
        axes[ii,0].set_title('Across {0:0.1f} km'.format(scales[ii] / 1e3))
            
        
    plt.savefig(tmp_direct + '/Pi_{0:04d}.png'.format(Itime), dpi=dpi)
    plt.close()
    
    
# If baroclinic transfers are there, use them.
if 'baroclinic_transfer' in results.variables:
    for Itime in range(rank, Ntime, num_procs):    

        timestamp = datetime.datetime.fromtimestamp(time[Itime])
        sup_title = "{0:02d} - {1:02d} - {2:04d} ( {3:02d}:{4:02d} )".format(
            timestamp.day, timestamp.month, timestamp.year, 
            timestamp.hour, timestamp.minute)
    
        # Initialize figure
        fig, axes = plt.subplots(num_scales, 1,
                sharex=True, sharey=True, squeeze=False,
                gridspec_kw = gridspec_props,
                figsize=(6, 4*num_scales))

        fig.suptitle(sup_title)
    
        # Plot each band
        for ii in range(num_scales):
        
            to_plot = bc_transfer[ii,Itime,:,:] * 1e6
            to_plot = np.ma.masked_where(mask==0, to_plot)
    
            m  = Basemap(ax = axes[ii,0], **map_settings)
    
            if np.max(np.abs(to_plot)) > 0:
                PlotTools.SignedLogPlot_onMap(LON * R2D, LAT * R2D, to_plot, axes[ii,0], fig, m, num_ords = 4, percentile=99.99)
    
            # Add coastlines, lat/lon lines, and draw the map
            m.drawcoastlines(linewidth=0.1)
            m.drawparallels(parallels, linewidth=0.5, labels=[0,0,0,0], color='k')
            m.drawmeridians(meridians, linewidth=0.5, labels=[0,0,0,0], color='k')
            m.pcolormesh(LON*R2D, LAT*R2D, mask, vmin=-1, vmax=1, cmap=mask_cmap, latlon=True)
            
            axes[ii,0].set_title('Across {0:0.1f} km'.format(scales[ii] / 1e3))
            
        
        plt.savefig(tmp_direct + '/Lambda_m_{0:04d}.png'.format(Itime), dpi=dpi)
        plt.close()


# If more than one time point, create mp4s
if Ntime > 1:
    PlotTools.merge_to_mp4(tmp_direct + '/Pi_%04d.png',    
            out_direct + '/Pi.mp4',    fps=12)
    if 'baroclinic_transfer' in source.variables:
        PlotTools.merge_to_mp4(tmp_direct + '/Lambda_m_%04d.png', 
                out_direct + '/Lambda_m.mp4', fps=12)
else:
    shutil.move(tmp_direct + '/Pi_0000.png', 
            out_direct + '/Pi.png')
    if 'p' in source.variables:
        shutil.move(tmp_direct + '/Lambda_m_0000.png', 
                out_direct + '/Lambda_m.png')

# Now delete the frames
shutil.rmtree(tmp_direct)

## Plot time averages if Ntime > 1

if (Ntime > 1):

    # Initialize figure
    fig, axes = plt.subplots(num_scales, 1,
            sharex=True, sharey=True, squeeze=False,
            gridspec_kw = gridspec_props,
            figsize=(6, 4*num_scales))

    fig.suptitle('Time average')
    
    # Plot each band
    for ii in range(num_scales):
        
        to_plot = np.mean(transfer[ii,:,:,:] * 1e6, axis=0)
        to_plot = np.ma.masked_where(mask==0, to_plot)
    
        m  = Basemap(ax = axes[ii,0], **map_settings)
    
        if np.max(np.abs(to_plot)) > 0:
            PlotTools.SignedLogPlot_onMap(LON * R2D, LAT * R2D, to_plot, axes[ii,0], fig, m, num_ords = 5)
    
        # Add coastlines, lat/lon lines, and draw the map
        m.drawcoastlines(linewidth=0.1)
        m.drawparallels(parallels, linewidth=0.5, labels=[0,0,0,0], color='k')
        m.drawmeridians(meridians, linewidth=0.5, labels=[0,0,0,0], color='k')
        m.pcolormesh(LON*R2D, LAT*R2D, mask, vmin=-1, vmax=1, cmap=mask_cmap, latlon=True)
            
        axes[ii,0].set_title('Across {0:0.1f} km'.format(scales[ii] / 1e3))
            
        
    plt.savefig(out_direct + '/AVE_Pi.png', dpi=dpi)
    plt.close()
    
    
if 'baroclinic_transfer' in results.variables:
    if (Ntime > 1):

        # Initialize figure
        fig, axes = plt.subplots(num_scales, 1,
                sharex=True, sharey=True, squeeze=False,
                gridspec_kw = gridspec_props,
                figsize=(6, 4*num_scales))

        fig.suptitle('Time average')
    
        # Plot each band
        for ii in range(num_scales):
        
            to_plot = np.mean(bc_transfer[ii,:,:,:] * 1e6, axis=0)
            to_plot = np.ma.masked_where(mask==0, to_plot)
    
            m  = Basemap(ax = axes[ii,0], **map_settings)
    
            if np.max(np.abs(to_plot)) > 0:
                PlotTools.SignedLogPlot_onMap(LON * R2D, LAT * R2D, to_plot, axes[ii,0], fig, m, num_ords = 4, percentile=99.99)
    
            # Add coastlines, lat/lon lines, and draw the map
            m.drawcoastlines(linewidth=0.1)
            m.drawparallels(parallels, linewidth=0.5, labels=[0,0,0,0], color='k')
            m.drawmeridians(meridians, linewidth=0.5, labels=[0,0,0,0], color='k')
            m.pcolormesh(LON*R2D, LAT*R2D, mask, vmin=-1, vmax=1, cmap=mask_cmap, latlon=True)
            
            axes[ii,0].set_title('Across {0:0.1f} km'.format(scales[ii] / 1e3))
        
        plt.savefig(out_direct + '/AVE_Lambda_m.png', dpi=dpi)
        plt.close()
