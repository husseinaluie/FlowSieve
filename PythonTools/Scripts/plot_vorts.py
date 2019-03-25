import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import cmocean, sys, PlotTools, os, shutil, datetime
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
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
time      = results.variables['time'][:] * (60*60) # convert hours to second
vort_r    = results.variables['vort_r'  ][:, :, 0, :, :]
vort_lat  = results.variables['vort_lat'][:, :, 0, :, :]
vort_lon  = results.variables['vort_lon'][:, :, 0, :, :]
mask      = results.variables['mask'][:]

# Do some time handling tp adjust the epochs
# appropriately
epoch = datetime.datetime(1950,1,1)   # the epoch of the time dimension
dt_epoch = datetime.datetime.fromtimestamp(0)  # the epoch used by datetime
epoch_delta = dt_epoch - epoch  # difference
time = time - epoch_delta.total_seconds()  # shift

num_scales = len(scales)
Ntime = len(time)

LON, LAT = np.meshgrid(longitude * D2R, latitude * D2R)

# Some parameters for plotting
map_settings = PlotTools.MapSettings(longitude, latitude)

meridians = np.round(np.linspace(longitude.min(), longitude.max(), 5))
parallels = np.round(np.linspace(latitude.min(),  latitude.max(),  5))

cbar_props     = dict(pad = 0.1, shrink = 0.85, orientation = 'horizontal')
gridspec_props = dict(wspace = 0.05, hspace = 0.05, left = 0.05, right = 0.95, bottom = 0.05, top = 0.95)


##
## Begin Plotting
##

## Vorticity binning
for Itime in range(rank, Ntime, num_procs):    

    timestamp = datetime.datetime.fromtimestamp(time[Itime])
    sup_title = "{0:02d} - {1:02d} - {2:04d} ( {3:02d}:{4:02d} )".format(
            timestamp.day, timestamp.month, timestamp.year, 
            timestamp.hour, timestamp.minute)
    
    # Initialize figure
    fig, axes = plt.subplots(1, num_scales,
            sharex=True, sharey=True, squeeze=False,
            gridspec_kw = gridspec_props,
            figsize=(4*num_scales, 14/3.))

    fig.suptitle(sup_title)
    
    # Plot each band
    for ii in range(num_scales):

        CV_r   = np.nanpercentile(np.abs(vort_r[  ii, :, :]), 99)
        CV_lon = np.nanpercentile(np.abs(vort_lon[ii, :, :]), 99)
        CV_lat = np.nanpercentile(np.abs(vort_lat[ii, :, :]), 99)

        for jj in range(1):  # Right now, just plot vort_r
        
            if jj == 0:
                to_plot = vort_r[ii,  Itime,:,:]
                CV = CV_r
            if jj == 1:
                to_plot = vort_lon[ii,Itime,:,:]
                CV = CV_lon
            if jj == 2:
                to_plot = vort_lat[ii,Itime,:,:]
                CV = CV_lat

            to_plot = np.ma.masked_where(mask==0, to_plot)
        
            m  = Basemap(ax = axes[jj,ii], **map_settings)
    
            qm  = m.pcolormesh(LON*R2D, LAT*R2D, to_plot, 
                    cmap='cmo.balance', 
                    vmin = -CV, vmax = CV, latlon = True)
        
            cbar = plt.colorbar(qm, ax = axes[jj,ii], **cbar_props)
            PlotTools.ScientificCbar(cbar, units='$\mathrm{s}^{-1}$', orientation='horizontal')
    
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
        
    
    plt.savefig(tmp_direct + '/vorticity_bands_{0:04d}.png'.format(Itime), dpi=dpi)
    plt.close()


## Dichotomies
for Itime in range(rank, Ntime, num_procs):    

    timestamp = datetime.datetime.fromtimestamp(time[Itime])
    sup_title = "{0:02d} - {1:02d} - {2:04d} ( {3:02d}:{4:02d} )".format(
            timestamp.day, timestamp.month, timestamp.year, 
            timestamp.hour, timestamp.minute)

    # Initialize figure
    fig, axes = plt.subplots(num_scales-1, 2,
            sharex=True, sharey=True, squeeze=False,
            gridspec_kw = gridspec_props,
            figsize=(10, 4*num_scales-1))

    fig.suptitle(sup_title)
    
    # Plot each band
    for ii in range(num_scales-1):
        
        to_plot_below = np.sum(vort_r[:ii+1, Itime,:,:], axis=0)
        to_plot_above = np.sum(vort_r[ ii+1:,Itime,:,:], axis=0)
    
        to_plot_below = np.ma.masked_where(mask==0, to_plot_below)
        to_plot_above = np.ma.masked_where(mask==0, to_plot_above)
    
        m_a = Basemap(ax = axes[ii,0], **map_settings)
        m_b = Basemap(ax = axes[ii,1], **map_settings)
    
        CV_a = np.nanpercentile(np.abs(to_plot_above), 99)
        CV_b = np.nanpercentile(np.abs(to_plot_below), 99)
    
        vmax = max(CV_a, CV_b)
    
        qm_a  = m_a.pcolormesh(LON*R2D, LAT*R2D, to_plot_above, 
                    cmap='cmo.balance', 
                    vmin = -CV_a, vmax = CV_a, latlon = True)
        qm_b  = m_b.pcolormesh(LON*R2D, LAT*R2D, to_plot_below, 
                    cmap='cmo.balance', 
                    vmin = -CV_b, vmax = CV_b, latlon = True)

        cbar_a = plt.colorbar(qm_a, ax = axes[ii,0], **cbar_props)
        cbar_b = plt.colorbar(qm_b, ax = axes[ii,1], **cbar_props)
        PlotTools.ScientificCbar(cbar_a, units='$\mathrm{s}^{-1}$', orientation='horizontal')
        PlotTools.ScientificCbar(cbar_b, units='$\mathrm{s}^{-1}$', orientation='horizontal')
    
        # Add coastlines and lat/lon lines
        for m in [m_a, m_b]:
            m.drawcoastlines(linewidth=0.1)
            m.drawparallels(parallels, linewidth=0.5, labels=[0,0,0,0], color='g')
            m.drawmeridians(meridians, linewidth=0.5, labels=[0,0,0,0], color='g')
            m.pcolormesh(LON*R2D, LAT*R2D, mask, vmin=-1, vmax=1, cmap=mask_cmap, latlon=True)
        
        axes[ii,0].set_ylabel('{0:.1f} km'.format(scales[ii] / 1e3))
    
    axes[0,0].set_title('Coarse $(>l)$')
    axes[0,1].set_title('Fine $(<l)$')
        
    plt.savefig(tmp_direct + '/vorticity_dichotomies_{0:04d}.png'.format(Itime), dpi=dpi)
    plt.close()



# If more than one time point, create mp4s
if Ntime > 1:
    PlotTools.merge_to_mp4(tmp_direct + '/vorticity_bands_%04d.png',    
            out_direct + '/vorticity_bands.mp4', fps=12)
    PlotTools.merge_to_mp4(tmp_direct + '/vorticity_dichotomies_%04d.png',    
            out_direct + '/vorticity_dichotomies.mp4', fps=12)

    # Now delete the frames
    shutil.rmtree(tmp_direct)
    
else:
    shutil.move(tmp_direct + '/vorticity_bands_0000.png',
            out_direct + '/vorticity_bands.png')
    shutil.move(tmp_direct + '/vorticity_dichotomies_0000.png',
            out_direct + '/vorticity_dichotomies.png')
