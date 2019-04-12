import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import PlotTools, cmocean
import matpy as mp
import scipy.integrate as spi
import sys, os, shutil, datetime, glob

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

source  = Dataset('input.nc', 'r')
dAreas = PlotTools.getAreas(source)

rho0    = 1e3
eps     = 1e-10

# Loop through filters
for fp in files:
    with Dataset(fp, 'r') as results:
        scale = results.filter_scale
        if (rank == 0):
            print('  {0:.4g}km'.format(scale/1e3))

        time  = results.variables['time'][:] * (60*60) # convert hours to seconds
        Ntime = len(time)

        if (Ntime == 1):
            if (rank == 0):
                print("Not enough time points to plot KE flux")
            sys.exit()

        # Do some time handling tp adjust the epochs
        # appropriately
        epoch       = datetime.datetime(1950,1,1)   # the epoch of the time dimension
        dt_epoch    = datetime.datetime.fromtimestamp(0)  # the epoch used by datetime
        epoch_delta = dt_epoch - epoch  # difference
        time        = time - epoch_delta.total_seconds()  # shift

        if Ntime >= 10:
            order = 8
        elif Ntime >= 5:
            order = 4
        else:
            order = max(1, Ntime - 2)
        Dt = mp.FiniteDiff(time, order, spb=False)

        do_PEtoKE = ('PEtoKE' in results.variables)
        do_divJ   = ('div_Jtransport' in results.variables)

        depth  = results.variables['depth']
        Ndepth = len(depth)

        lat  = results.variables['latitude']
        Nlat = len(lat)

        lon  = results.variables['longitude']
        Nlon = len(lon)

        dArea = np.tile(dAreas.reshape(1,Nlat,Nlon), (Ndepth,1,1)) 

        u_r       = results.variables['coarse_u_r'  ][:]
        u_lon     = results.variables['coarse_u_lon'][:]
        u_lat     = results.variables['coarse_u_lat'][:]

        KE_from_vel = 0.5 * rho0 * (  u_r**2 + u_lat**2 + u_lon**2 )
        KE_from_vel = KE_from_vel.reshape(Ntime, Ndepth*Nlat*Nlon)
        KE_flux = np.dot(Dt, KE_from_vel).reshape(Ntime, Ndepth, Nlat, Nlon)
        KE_from_vel = KE_from_vel.reshape(Ntime, Ndepth, Nlat, Nlon)

        for iT in range(rank, Ntime, num_procs):

            timestamp = datetime.datetime.fromtimestamp(time[iT])
            sup_title = "{0:02d} - {1:02d} - {2:04d} ( {3:02d}:{4:02d} )".format(
                    timestamp.day, timestamp.month, timestamp.year, 
                    timestamp.hour, timestamp.minute)
            Pi        = results.variables['energy_transfer'][iT,:,:,:]

            KE_sel =  KE_flux[iT,:][~KE_flux[iT,:].mask].ravel()
            Pi_sel = -Pi[~Pi.mask].ravel()
            if do_PEtoKE:
                PEtoKE  = results.variables['PEtoKE'][iT,:,:,:]
                Pi_sel += PEtoKE[~PEtoKE.mask].ravel()
            if do_divJ:
                div_J   = results.variables['div_Jtransport'][iT,:,:,:]
                Pi_sel += -div_J[~div_J.mask].ravel()
        
            label = '$-\Pi$' 
            if do_divJ:
                label = label + '$-\\nabla\cdot J_{\mathrm{transport}}$'
            if do_PEtoKE:
                label = label + '$+\overline{\\rho}g\overline{u}_r$'
            label = label + ' $(\mathrm{W}\cdot\mathrm{m}^{-3})$'
    
            # Then plot
            # Initialize figure
            fig, axes = plt.subplots(2, 2, 
                    gridspec_kw = dict(left = 0.15, right = 0.95, bottom = 0.1, top = 0.9,
                        hspace=0.02, wspace=0.02),
                    figsize=(7.5, 6) )
        
            PlotTools.SignedLogScatter_hist(Pi_sel, KE_sel, axes,
                    force_equal = True, nbins_x = 100, nbins_y = 100)
        
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
            plt.figtext(mid_x, 0.05, label, horizontalalignment='center', 
                    verticalalignment='top', rotation='horizontal', fontsize=16)

            # ylabel
            mid_y = 0.5 * ( axes[0,0].get_position().y0 + axes[1,1].get_position().y1 )
            plt.figtext(0.05, mid_y, 
                '$\\frac{d}{dt}\left( \\frac{\\rho_0}{2}\overline{u}\cdot\overline{u} \\right)$ $(\mathrm{W}\cdot\mathrm{m}^{-3})$',
                horizontalalignment='right', verticalalignment='center', 
                rotation='vertical', fontsize=16)
        
            plt.savefig(tmp_direct + '/{0:.4g}_KE_fluxes_{1:04d}.png'.format(scale/1e3,iT), dpi=dpi)
            plt.close()

if (rank > 0):
    sys.exit()

# If more than one time point, create mp4s
if Ntime > 1:
    for fp in files:
        with Dataset(fp, 'r') as results:
            scale = results.filter_scale
            PlotTools.merge_to_mp4(
                tmp_direct + '/{0:.4g}_KE_fluxes_%04d.png'.format(scale/1e3),
                out_direct + '/{0:.4g}km/KE_fluxes.mp4'.format(scale/1e3),
                fps=12)
else:
    for fp in files:
        with Dataset(fp, 'r') as results:
            scale = results.filter_scale
            shutilmove(tmp_direct + '/{0:.4g}_KE_fluxes_%04d.png'.format(scale/1e3),
                out_direct + '/{0:.4g}km/KE_fluxes.mp4'.format(scale/1e3))

# Now delete the frames
shutil.rmtree(tmp_direct)
