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

        net_EN      = np.zeros((Ntime,))
        net_EN_proc = np.zeros((Ntime,))
        net_Lambda      = np.zeros((Ntime,))
        net_Lambda_proc = np.zeros((Ntime,))

        depth  = results.variables['depth']
        Ndepth = len(depth)

        lat  = results.variables['latitude']
        Nlat = len(lat)

        lon  = results.variables['longitude']
        Nlon = len(lon)

        dArea = np.tile(dAreas.reshape(1,Nlat,Nlon), (Ndepth,1,1)) 

        vort_r   = results.variables['coarse_vort_r'  ][:]
        #vort_lon = results.variables['coarse_vort_lon'][:]
        #vort_lat = results.variables['coarse_vort_lat'][:]

        EN_from_vort = 0.5 * rho0 * (  vort_r**2 )
        #EN_from_vort = 0.5 * rho0 * (  vort_r**2 + vort_lat**2 + vort_lon**2 )
        EN_from_vort = EN_from_vort.reshape(Ntime, Ndepth*Nlat*Nlon)
        EN_flux = np.dot(Dt, EN_from_vort).reshape(Ntime, Ndepth, Nlat, Nlon)
        EN_from_vort = EN_from_vort.reshape(Ntime, Ndepth, Nlat, Nlon)

        for iT in range(rank, Ntime, num_procs):

            timestamp = datetime.datetime.fromtimestamp(time[iT])
            sup_title = "{0:02d} - {1:02d} - {2:04d} ( {3:02d}:{4:02d} )".format(
                    timestamp.day, timestamp.month, timestamp.year, 
                    timestamp.hour, timestamp.minute)
            Lambda = results.variables['Lambda_m'][iT,:,:,:]

            net_EN_proc[iT] = np.sum(EN_from_vort[iT,:,:,:] * dArea)
            net_Lambda_proc[iT] = np.sum(Lambda * dArea)
    
            EN_sel     =  EN_flux[iT,:][~EN_flux[iT,:].mask].ravel()
            Lambda_sel = -Lambda[~Lambda.mask].ravel()
        
            label = '$\Lambda^m$' 
            label = label + ' $(\mathrm{W}^{\omega}\cdot\mathrm{m}^{-3})$'
    
            # Then plot
            # Initialize figure
            fig, axes = plt.subplots(2, 2, 
                    gridspec_kw = dict(left = 0.15, right = 0.95, bottom = 0.1, top = 0.9,
                        hspace=0.02, wspace=0.02),
                    figsize=(7.5, 6) )
        
            PlotTools.SignedLogScatter_hist(Lambda_sel, EN_sel, axes,
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
                '$\\frac{d}{dt}\left( \\frac{\\rho_0}{2}\overline{\omega}\cdot\overline{\omega} \\right)$ $(\mathrm{W}\cdot\mathrm{m}^{-3})$',
                horizontalalignment='right', verticalalignment='center', 
                rotation='vertical', fontsize=16)
        
            plt.savefig(tmp_direct + '/{0:.4g}_EN_fluxes_{1:04d}.png'.format(scale/1e3,iT), dpi=dpi)
            plt.close()

    
        comm.Allreduce(net_EN_proc, net_EN, op=MPI.SUM)
        comm.Allreduce(net_Lambda_proc, net_Lambda, op=MPI.SUM)
        if rank == 0:
            # Now plot the space-integrated version
            colours = plt.rcParams['axes.prop_cycle'].by_key()['color']
            fig = plt.figure()
            ax1    = fig.add_axes([0.1, 0.25, 0.8, 0.7])
            ax2    = fig.add_axes([0.1, 0.15, 0.8, 0.1])
            leg_ax = fig.add_axes([0.1, 0.05, 0.8, 0.1])

            net_EN_flux = np.dot(Dt, net_EN)

            to_plot = net_EN_flux
            label='$\int_{\Omega}\\frac{d}{dt}\left( \\frac{\\rho_0}{2}\overline{\omega}\cdot\overline{\omega} \\right)\mathrm{dA}$'
            ax1.plot(np.ma.masked_where(to_plot<0, time),  np.ma.masked_where(to_plot<0, np.abs(to_plot)), '-',  color=colours[0], label=label)
            ax1.plot(np.ma.masked_where(to_plot>0, time),  np.ma.masked_where(to_plot>0, np.abs(to_plot)), '--', color=colours[0])#, label=label)

            to_plot = net_Lambda
            label='$\int_{\Omega}\Lambda^m\mathrm{dA}$'
            ax1.plot(np.ma.masked_where(to_plot<0, time),  np.ma.masked_where(to_plot<0, np.abs(to_plot)), '-',  color=colours[1], label=label)
            ax1.plot(np.ma.masked_where(to_plot>0, time),  np.ma.masked_where(to_plot>0, np.abs(to_plot)), '--', color=colours[1])#, label=label)

            ax1.set_yscale('log')
            ax1.set_ylabel('$\mathrm{W}^{\omega}$')

            ax2.set_xlim(ax1.get_xlim())

            ax1.set_xticks([])
            ax2.set_xticks([])
            ax2.set_yticks([])

            PlotTools.LabelTimeAxis(ax2, time)

            # Add legend
            handles, labels = ax1.get_legend_handles_labels()
            leg_ax.legend(handles, labels, bbox_to_anchor=(0., 0., 1., 1.), ncol=3, mode='expand',
                    frameon = True, borderaxespad=0.)
            leg_ax.set_xticks([])
            leg_ax.set_yticks([])
            for pos in ['left', 'right', 'top', 'bottom']:
                leg_ax.spines[pos].set_visible(False)

            plt.savefig(out_direct + '/{0:.4g}km/EN_fluxes_net.pdf'.format(scale/1e3))
            plt.close()

            # Also plot the time-integrated version
            colours = plt.rcParams['axes.prop_cycle'].by_key()['color']
            fig = plt.figure()
            ax1    = fig.add_axes([0.1, 0.25, 0.8, 0.7])
            ax2    = fig.add_axes([0.1, 0.15, 0.8, 0.1])
            leg_ax = fig.add_axes([0.1, 0.05, 0.8, 0.1])

            label='$\int_{\Omega}\left( \\frac{\\rho_0}{2}\overline{\omega}\cdot\overline{\omega} \\right)\mathrm{dA}$'
            to_plot = net_EN
            ax1.plot(np.ma.masked_where(to_plot<0, time),  np.ma.masked_where(to_plot<0, np.abs(to_plot)), '-',  color=colours[0], label=label)
            ax1.plot(np.ma.masked_where(to_plot>0, time),  np.ma.masked_where(to_plot>0, np.abs(to_plot)), '--', color=colours[0])#, label=label)

            label='$\int^t\int_{\Omega}\Lambda^m\mathrm{dA}$'
            to_plot = spi.cumtrapz(net_Lambda, time)
            ax1.plot(np.ma.masked_where(to_plot<0, time[1:]),  np.ma.masked_where(to_plot<0, np.abs(to_plot)), '-',  color=colours[1], label=label)
            ax1.plot(np.ma.masked_where(to_plot>0, time[1:]),  np.ma.masked_where(to_plot>0, np.abs(to_plot)), '--', color=colours[1])#, label=label)

            ax1.set_yscale('log')
            ax1.set_ylabel('$\mathrm{J}^{\omega}$')

            ax2.set_xlim(ax1.get_xlim())

            ax1.set_xticks([])
            ax2.set_xticks([])
            ax2.set_yticks([])

            PlotTools.LabelTimeAxis(ax2, time)

            # Add legend
            handles, labels = ax1.get_legend_handles_labels()
            leg_ax.legend(handles, labels, bbox_to_anchor=(0., 0., 1., 1.), ncol=3, mode='expand',
                    frameon = True, borderaxespad=0.)
            leg_ax.set_xticks([])
            leg_ax.set_yticks([])
            for pos in ['left', 'right', 'top', 'bottom']:
                leg_ax.spines[pos].set_visible(False)

            plt.savefig(out_direct + '/{0:.4g}km/EN_net.pdf'.format(scale/1e3))
            plt.close()

if (rank > 0):
    sys.exit()

# If more than one time point, create mp4s
if Ntime > 1:
    for fp in files:
        with Dataset(fp, 'r') as results:
            scale = results.filter_scale
            PlotTools.merge_to_mp4(
                tmp_direct + '/{0:.4g}_EN_fluxes_%04d.png'.format(scale/1e3),
                out_direct + '/{0:.4g}km/EN_fluxes.mp4'.format(scale/1e3),
                fps=12)
else:
    for fp in files:
        with Dataset(fp, 'r') as results:
            scale = results.filter_scale
            shutilmove(tmp_direct + '/{0:.4g}_EN_fluxes_%04d.png'.format(scale/1e3),
                out_direct + '/{0:.4g}km/EN_fluxes.mp4'.format(scale/1e3))

# Now delete the frames
shutil.rmtree(tmp_direct)
