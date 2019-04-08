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

        net_KE      = np.zeros((Ntime,))
        net_KE_proc = np.zeros((Ntime,))
        net_Pi      = np.zeros((Ntime,))
        net_Pi_proc = np.zeros((Ntime,))
        net_div_vel      = np.zeros((Ntime,))
        net_div_vel_proc = np.zeros((Ntime,))
        if do_PEtoKE:
            net_PEtoKE      = np.zeros((Ntime,))
            net_PEtoKE_proc = np.zeros((Ntime,))
        if do_divJ:
            net_divJ      = np.zeros((Ntime,))
            net_divJ_proc = np.zeros((Ntime,))

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
            div_vel   = results.variables['full_vel_div'][iT,:,:,:]

            net_KE_proc[iT] = np.sum(KE_from_vel[iT,:,:,:] * dArea)
            net_Pi_proc[iT] = np.sum(- Pi * dArea)
            if do_PEtoKE:
                PEtoKE = results.variables['PEtoKE'][iT,:,:,:]
                net_PEtoKE_proc[iT] =  np.sum(PEtoKE * dArea)
            if do_divJ:
                div_J = results.variables['div_Jtransport'][iT,:,:,:]
                net_divJ_proc[iT] = -np.sum(div_J * dArea)
    
            KE_sel =  KE_flux[iT,:][~KE_flux[iT,:].mask].ravel()
            Pi_sel = -Pi[~Pi.mask].ravel()
            if do_PEtoKE:
                Pi_sel += PEtoKE[~PEtoKE.mask].ravel()
            if do_divJ:
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

    
        comm.Allreduce(net_KE_proc, net_KE, op=MPI.SUM)
        comm.Allreduce(net_Pi_proc, net_Pi, op=MPI.SUM)
        comm.Allreduce(net_div_vel_proc, net_div_vel, op=MPI.SUM)
        if do_PEtoKE:
            comm.Allreduce(net_PEtoKE_proc, net_PEtoKE, op=MPI.SUM)
        if do_divJ:
            comm.Allreduce(net_divJ_proc, net_divJ, op=MPI.SUM)
        if rank == 0:
            # Now plot the space-integrated version
            colours = plt.rcParams['axes.prop_cycle'].by_key()['color']
            fig = plt.figure()
            ax1    = fig.add_axes([0.1, 0.25, 0.8, 0.7])
            ax2    = fig.add_axes([0.1, 0.15, 0.8, 0.1])
            leg_ax = fig.add_axes([0.1, 0.05, 0.8, 0.1])

            net_KE_flux = np.dot(Dt, net_KE)

            to_plot = net_KE_flux
            label='$\int_{\Omega}\\frac{d}{dt}\left( \\frac{\\rho_0}{2}\overline{u}\cdot\overline{u} \\right)\mathrm{dA}$'
            ax1.plot(np.ma.masked_where(to_plot<0, time),  np.ma.masked_where(to_plot<0, np.abs(to_plot)), '-',  color=colours[0], label=label)
            ax1.plot(np.ma.masked_where(to_plot>0, time),  np.ma.masked_where(to_plot>0, np.abs(to_plot)), '--', color=colours[0])#, label=label)

            to_plot = net_Pi 
            label='$-\int_{\Omega}\Pi\mathrm{dA}$'
            ax1.plot(np.ma.masked_where(to_plot<0, time),  np.ma.masked_where(to_plot<0, np.abs(to_plot)), '-',  color=colours[1], label=label)
            ax1.plot(np.ma.masked_where(to_plot>0, time),  np.ma.masked_where(to_plot>0, np.abs(to_plot)), '--', color=colours[1])#, label=label)

            if do_PEtoKE:
                label='$\int_{\Omega}\overline{\\rho}g\overline{u}_r\mathrm{dA}$'
                to_plot = net_PEtoKE
                ax1.plot(np.ma.masked_where(to_plot<0, time),  np.ma.masked_where(to_plot<0, np.abs(to_plot)), '-',  color=colours[2], label=label)
                ax1.plot(np.ma.masked_where(to_plot>0, time),  np.ma.masked_where(to_plot>0, np.abs(to_plot)), '--', color=colours[2])#, label=label)

            if do_divJ:
                label='$-\int_{\Omega}\\nabla\cdot J\mathrm{dA}$'
                to_plot = net_divJ
                ax1.plot(np.ma.masked_where(to_plot<0, time),  np.ma.masked_where(to_plot<0, np.abs(to_plot)), '-',  color=colours[3], label=label)
                ax1.plot(np.ma.masked_where(to_plot>0, time),  np.ma.masked_where(to_plot>0, np.abs(to_plot)), '--', color=colours[3])#, label=label)

            ax1.set_yscale('log')
            ax1.set_ylabel('$\mathrm{W}$')

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

            plt.savefig(out_direct + '/{0:.4g}km/KE_fluxes_net.pdf'.format(scale/1e3))
            plt.close()

            # Also plot the time-integrated version
            colours = plt.rcParams['axes.prop_cycle'].by_key()['color']
            fig = plt.figure()
            ax1    = fig.add_axes([0.1, 0.25, 0.8, 0.7])
            ax2    = fig.add_axes([0.1, 0.15, 0.8, 0.1])
            leg_ax = fig.add_axes([0.1, 0.05, 0.8, 0.1])

            label='$\int_{\Omega}\left( \\frac{\\rho_0}{2}\overline{u}\cdot\overline{u} \\right)\mathrm{dA}$'
            to_plot = net_KE
            ax1.plot(np.ma.masked_where(to_plot<0, time),  np.ma.masked_where(to_plot<0, np.abs(to_plot)), '-',  color=colours[0], label=label)
            ax1.plot(np.ma.masked_where(to_plot>0, time),  np.ma.masked_where(to_plot>0, np.abs(to_plot)), '--', color=colours[0])#, label=label)

            label='$-\int^t\int_{\Omega}\Pi\mathrm{dA}$'
            to_plot = spi.cumtrapz(net_Pi, time)
            ax1.plot(np.ma.masked_where(to_plot<0, time[1:]),  np.ma.masked_where(to_plot<0, np.abs(to_plot)), '-',  color=colours[1], label=label)
            ax1.plot(np.ma.masked_where(to_plot>0, time[1:]),  np.ma.masked_where(to_plot>0, np.abs(to_plot)), '--', color=colours[1])#, label=label)

            if do_PEtoKE:
                label='$\int^t\int_{\Omega}\overline{\\rho}g\overline{u}_r\mathrm{dA}$'
                to_plot = spi.cumtrapz(net_PEtoKE, time)
                ax1.plot(np.ma.masked_where(to_plot<0, time[1:]),  np.ma.masked_where(to_plot<0, np.abs(to_plot)), '-',  color=colours[2], label=label)
                ax1.plot(np.ma.masked_where(to_plot>0, time[1:]),  np.ma.masked_where(to_plot>0, np.abs(to_plot)), '--', color=colours[2])#, label=label)

            if do_divJ:
                label='$-\int^t\int_{\Omega}\\nabla\cdot J\mathrm{dA}$'
                to_plot = spi.cumtrapz(net_divJ, time)
                ax1.plot(np.ma.masked_where(to_plot<0, time[1:]),  np.ma.masked_where(to_plot<0, np.abs(to_plot)), '-',  color=colours[3], label=label)
                ax1.plot(np.ma.masked_where(to_plot>0, time[1:]),  np.ma.masked_where(to_plot>0, np.abs(to_plot)), '--', color=colours[3])#, label=label)

            ax1.set_yscale('log')
            ax1.set_ylabel('$\mathrm{J}$')

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

            plt.savefig(out_direct + '/{0:.4g}km/KE_net.pdf'.format(scale/1e3))
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
