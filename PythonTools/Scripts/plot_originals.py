import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import cmocean, sys, datetime
from netCDF4 import Dataset
import PlotTools, subprocess, shutil, os
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

source = Dataset('input.nc', 'r')

fp = 'filter_output.nc'
results = Dataset(fp, 'r')

R_earth = 6371e3
D2R     = np.pi / 180
R2D     = 180 / np.pi
eps     = 1e-10

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
mask      = results.variables['mask'][:]

# Do some time handling tp adjust the epochs
# appropriately
epoch = datetime.datetime(1950,1,1)   # the epoch of the time dimension
dt_epoch = datetime.datetime.fromtimestamp(0)  # the epoch used by datetime
epoch_delta = dt_epoch - epoch  # difference
time = time - epoch_delta.total_seconds()  # shift

num_scales = len(scales)-1
Ntime = len(time)

LON, LAT = np.meshgrid(longitude * D2R, latitude * D2R)

uo = source.variables['uo'][:, 0, :, :]
vo = source.variables['vo'][:, 0, :, :]
Full_KE = 0.5 * (uo**2 + vo**2)


# Some parameters for plotting
proj = PlotTools.MapProjection(longitude, latitude)
Xp, Yp = proj(LON * R2D, LAT * R2D, inverse=False)

meridians = np.round(np.linspace(longitude.min(), longitude.max(), 5))
parallels = np.round(np.linspace(latitude.min(),  latitude.max(),  5))

cbar_props = dict(pad = 0.1, shrink = 0.85, orientation = 'vertical')
gridspec_props = dict(wspace = 0.05, hspace = 0.05, left = 0.1, right = 0.9, bottom = 0.1, top = 0.9)

def plot(LON, LAT, to_plot, filename, 
        one_sided=False, vmin=None, vmax=None, units='',
        cmap=None, title=None):

    plt.figure()
    ax = plt.subplot(1,1,1)
    
    if cmap == None:
        if one_sided:
            cmap = 'cmo.amp'
        else:
            cmap = 'cmo.balance'
    
    if (vmin == None) or (vmax == None):
        if one_sided:
            vmax = np.nanmax(to_plot)
            vmin = np.nanmin(to_plot)
        else:
            vmax = np.nanmax(np.abs(to_plot))
            vmin = -vmax

    qm = ax.pcolormesh(Xp, Yp, to_plot, cmap=cmap, vmin = vmin, vmax = vmax)
        
    cbar = plt.colorbar(qm, ax = ax, **cbar_props)
    PlotTools.ScientificCbar(cbar, units=units)
    
    # Add coastlines and lat/lon lines
    ax.pcolormesh(Xp, Yp, mask, vmin=-1, vmax=1, cmap=mask_cmap)
    PlotTools.AddParallels_and_Meridians(ax, proj, 
            parallels, meridians, latitude, longitude)

    if not(title == None):
        ax.set_title(title)
    
    plt.savefig(filename, dpi=dpi)
    plt.close()

for Itime in range(rank, Ntime, num_procs):    

    timestamp = datetime.datetime.fromtimestamp(time[Itime])
    sup_title = "{0:02d} - {1:02d} - {2:04d} ( {3:02d}:{4:02d} )".format(
            timestamp.day, timestamp.month, timestamp.year, 
            timestamp.hour, timestamp.minute)

    # Plot KE
    to_plot = Full_KE[Itime,:,:]
    to_plot = np.ma.masked_where(mask==0, to_plot)

    CV = np.max(np.abs(Full_KE))

    plot(LON*R2D, LAT*R2D, to_plot, tmp_direct + '/KE_{0:04d}.png'.format(Itime),
            one_sided=True, vmin = 0, vmax = CV, cmap='cmo.dense', title=sup_title)
    
    ## Full uo
    to_plot = uo[Itime,:,:]
    to_plot = np.ma.masked_where(mask==0, to_plot)

    CV = np.max(np.abs(uo))

    plot(LON*R2D, LAT*R2D, to_plot, tmp_direct + '/u_lon_{0:04d}.png'.format(Itime),
            one_sided=False, vmin = -CV, vmax = CV, title=sup_title)
        
    ## Full vo
    to_plot = vo[Itime,:,:]
    to_plot = np.ma.masked_where(mask==0, to_plot)

    CV = np.max(np.abs(vo))

    plot(LON*R2D, LAT*R2D, to_plot, tmp_direct + '/u_lat_{0:04d}.png'.format(Itime),
            one_sided=False, vmin = -CV, vmax = CV, title=sup_title)
    
    ## Rho, if available
    if 'rho' in source.variables:
        to_plot = source.variables['rho'][Itime,0,:,:]
        to_plot = np.ma.masked_where(mask==0, to_plot)

        mu = np.mean(source.variables['rho'][:,0,:,:])
        CV = np.max(np.abs(source.variables['rho'][:,0,:,:] - mu))

        plot(LON*R2D, LAT*R2D, to_plot-mu, tmp_direct + '/rho_{0:04d}.png'.format(Itime),
                one_sided=False, vmin=-CV, vmax=CV, title=sup_title, units=' + {0:.4g}'.format(mu))
    
    ## Pressure, if available
    if 'p' in source.variables:
        to_plot = source.variables['p'][Itime,0,:,:]
        to_plot = np.ma.masked_where(mask==0, to_plot)

        mu = np.mean(source.variables['p'][:,0,:,:])
        CV = np.max(np.abs(source.variables['p'][:,0,:,:] - mu))

        plot(LON*R2D, LAT*R2D, to_plot, tmp_direct + '/pressure_{0:04d}.png'.format(Itime), 
                one_sided=True, vmin=mu-CV, vmax=mu+CV, cmap='cmo.thermal', title=sup_title)
    

# If more than one time point, create mp4s
if Ntime > 1:
    PlotTools.merge_to_mp4(tmp_direct + '/KE_%04d.png',    
            out_direct + '/KE.mp4',    fps=12)
    PlotTools.merge_to_mp4(tmp_direct + '/u_lon_%04d.png', 
            out_direct + '/u_lon.mp4', fps=12)
    PlotTools.merge_to_mp4(tmp_direct + '/u_lat_%04d.png', 
            out_direct + '/u_lat.mp4', fps=12)
    if 'rho' in source.variables:
        PlotTools.merge_to_mp4(tmp_direct + '/rho_%04d.png', 
                out_direct + '/rho.mp4', fps=12)
    if 'p' in source.variables:
        PlotTools.merge_to_mp4(tmp_direct + '/pressure_%04d.png', 
                out_direct + '/pressure.mp4', fps=12)
else:
    shutil.move(tmp_direct + '/KE_0000.png',    
            out_direct + '/KE.png')
    shutil.move(tmp_direct + '/u_lon_0000.png', 
            out_direct + '/u_lon.png')
    shutil.move(tmp_direct + '/u_lat_0000.png', 
            out_direct + '/u_lat.png')
    if 'rho' in source.variables:
        shutil.move(tmp_direct + '/rho_0000.png', 
                out_direct + '/rho.png')
    if 'p' in source.variables:
        shutil.move(tmp_direct + '/pressure_0000.png', 
                out_direct + '/pressure.png')

# Now delete the frames
shutil.rmtree(tmp_direct)
