import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import PlotTools, cmocean
import matpy as mp
import sys, os, shutil, datetime

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

rho0    = 1e3
R_earth = 6371e3
D2R     = np.pi / 180
R2D     = 180 / np.pi
eps     = 1e-10

# Get the grid from the first filter
latitude  = results.variables['latitude'][:] * R2D
longitude = results.variables['longitude'][:] * R2D
scales    = results.variables['scale'][:]
depth     = results.variables['depth'][:]
time      = results.variables['time'][:] * (60 * 60) # time was in hours
mask      = results.variables['mask'][:]
v_r       = results.variables['vort_r'  ][:]
v_lon     = results.variables['vort_lon'][:]
v_lat     = results.variables['vort_lat'][:]

if not('baroclinic_transfer' in results.variables):
    sys.exit()

Lambda    = results.variables['baroclinic_transfer'][:]

Lambda = Lambda[:-1,:,:,:,:]

# Do some time handling tp adjust the epochs
# appropriately
epoch = datetime.datetime(1950,1,1)   # the epoch of the time dimension
dt_epoch = datetime.datetime.fromtimestamp(0)  # the epoch used by datetime
epoch_delta = dt_epoch - epoch  # difference
time = time - epoch_delta.total_seconds()  # shift

Nscales, Ntime, Ndepth, Nlat, Nlon = v_r.shape

if (Ntime == 1):
    print("Not enough time points to do enstrophy flux")
    sys.exit()

dlat = (latitude[1]  - latitude[0] ) * D2R
dlon = (longitude[1] - longitude[0]) * D2R
LON, LAT = np.meshgrid(longitude * D2R, latitude * D2R)
if source.variables['latitude'].units == 'm':
    dAreas = dlat * dlon * np.ones(LAT.shape)
else:
    dAreas = R_earth**2 * np.cos(LAT) * dlat * dlon

dArea = np.tile((mask*dAreas).reshape(1,Nlat,Nlon), (Ndepth,1,1)) 
mask  = np.tile(mask.reshape(1,Nlat,Nlon), (Ndepth, 1, 1))

if Ntime >= 5:
    order = 4
else:
    order = max(1, Ntime - 2)
Dt = mp.FiniteDiff(time, order, spb=False)

net_EN_flux = np.zeros((Ntime, Nscales-1))
net_lambda  = np.zeros((Ntime, Nscales-1))

## Scatter of fluxes vs Lambda
for iS in range(Nscales - 1):
    # First, sort out data
    EN_from_vort = 0.5 * rho0 * (  np.sum(v_r[  iS+1:,:,:,:,:], axis=0)**2 
                                 + np.sum(v_lat[iS+1:,:,:,:,:], axis=0)**2 
                                 + np.sum(v_lon[iS+1:,:,:,:,:], axis=0)**2 )
    EN_from_vort = EN_from_vort.reshape(Ntime, Ndepth*Nlat*Nlon)
    EN_flux = np.matmul(Dt, EN_from_vort)

    for iT in range(Ntime):

        timestamp = datetime.datetime.fromtimestamp(time[iT])
        sup_title = "{0:02d} - {1:02d} - {2:04d} ( {3:02d}:{4:02d} )".format(
                timestamp.day, timestamp.month, timestamp.year, 
                timestamp.hour, timestamp.minute)

        net_EN_flux[iT,iS] = np.sum(EN_flux[iT,:] * dArea.ravel() * mask.ravel())
        net_lambda[ iT,iS] = np.sum(Lambda[iS,iT,:] * dArea         * mask)

        EN_sel     = EN_flux[iT,:].ravel()[mask.ravel() == 1]
        Lambda_sel = Lambda[iS,iT,:,:,:].ravel()[mask.ravel() == 1]
    
        # Then plot
        # Initialize figure
        fig, axes = plt.subplots(2, 2, 
                gridspec_kw = dict(left = 0.1, right = 0.95, bottom = 0.1, top = 0.95,
                    hspace=0.02, wspace=0.02),
                figsize=(7.5, 6) )
        
        PlotTools.SignedLogScatter_hist(Lambda_sel, EN_sel, axes,
                force_equal = True, nbins_x = 200, nbins_y = 200)
        
        for ax in axes.ravel():
            xlim = ax.get_xlim()
            ylim = ax.get_ylim()
            ax.plot(xlim, xlim,'--c', label='$1:1$')
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
        
        axes[0,1].legend(loc='best')
        
        axes[0,0].set_xticklabels([])
        axes[0,1].set_xticklabels([])
        axes[0,1].set_yticklabels([])
        axes[1,1].set_yticklabels([])

        fig.suptitle(sup_title)

        # xlabel
        mid_x = 0.5 * ( axes[0,0].get_position().x0 + axes[1,1].get_position().x1 )
        plt.figtext(mid_x, 0.05, '$\Lambda^m$ $(\mathrm{W}^{\omega}\cdot\mathrm{m}^{-3})$',
             horizontalalignment='center', verticalalignment='top', rotation='horizontal', fontsize=16)

        # ylabel
        mid_y = 0.5 * ( axes[0,0].get_position().y0 + axes[1,1].get_position().y1 )
        plt.figtext(0.05, mid_y, '$\\frac{d}{dt}\left(\\frac{1}{2}\overline{\omega}\cdot\overline{\omega}\\right)$ $(\mathrm{W}^{\omega}\cdot\mathrm{m}^{-3})$',
           horizontalalignment='right', verticalalignment='center', rotation='vertical', fontsize=16)
        
        plt.savefig(tmp_direct + '/{0:.4g}_EN_fluxes_{1:04d}.png'.format(scales[iS]/1e3,iT), dpi=dpi)
        plt.close()
    
    # Now plot the space-integrated version
    colours = plt.rcParams['axes.prop_cycle'].by_key()['color']
    fig, axes = plt.subplots(1, 1, squeeze=False,
            gridspec_kw = dict(left = 0.15, right = 0.95, bottom = 0.1, top = 0.95,
            hspace=0.1))

    to_plot = net_EN_flux[:,iS] 
    label='$\int_{\Omega}\\frac{d}{dt}\left(\\frac{1}{2}\overline{\omega}\cdot\overline{\omega}\\right)$'
    axes[0,0].plot(np.ma.masked_where(to_plot<0, time),  np.ma.masked_where(to_plot<0, np.abs(to_plot)), '-',  color=colours[0], label=label)
    axes[0,0].plot(np.ma.masked_where(to_plot>0, time),  np.ma.masked_where(to_plot>0, np.abs(to_plot)), '--', color=colours[0])#, label=label)

    to_plot = net_lambda[:,iS] 
    label='$\int_{\Omega}\Lambda^m$'
    axes[0,0].plot(np.ma.masked_where(to_plot<0, time),  np.ma.masked_where(to_plot<0, np.abs(to_plot)), '-',  color=colours[1], label=label)
    axes[0,0].plot(np.ma.masked_where(to_plot>0, time),  np.ma.masked_where(to_plot>0, np.abs(to_plot)), '--', color=colours[1])#, label=label)

    axes[0,0].set_yscale('log')
    axes[0,0].legend(loc='best')
    axes[0,0].set_ylabel('$\mathrm{W}^{\omega}$')
    plt.savefig(out_direct + '/{0:.4g}km/EN_fluxes_net.pdf'.format(scales[iS]/1e3))
    plt.close()

# If more than one time point, create mp4s
if Ntime > 1:
    for iS in range(Nscales-1):
        PlotTools.merge_to_mp4(
                tmp_direct + '/{0:.4g}_EN_fluxes_%04d.png'.format(scales[iS]/1e3),
                out_direct + '/{0:.4g}km/EN_fluxes.mp4'.format(scales[iS]/1e3),
                fps=12)
else:
    for iS in range(Nscales-1):
        shutilmove(tmp_direct + '/{0:.4g}_EN_fluxes_%04d.png'.format(scales[iS]/1e3),
                out_direct + '/{0:.4g}km/EN_fluxes.mp4'.format(scales[iS]/1e3))

# Now delete the frames
shutil.rmtree(tmp_direct)
