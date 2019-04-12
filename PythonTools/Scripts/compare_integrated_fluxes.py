import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import PlotTools, cmocean
import matpy as mp
import scipy.integrate as spi
import sys, os, shutil, datetime, glob

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

    if not(os.path.exists(out_direct)):
        os.makedirs(out_direct)

source  = Dataset('input.nc', 'r')
dAreas = PlotTools.getAreas(source)

# Sort by filter length
scales = np.zeros(len(files))
files2 = []
for (fp,ii) in zip(files, range(len(files))):
    with Dataset(fp, 'r') as results:
        scales[ii] = results.filter_scale

inds = np.argsort(scales)
for ii in range(len(files)):
    files2 += [files[inds[ii]]]
files = files2

rho0    = 1025.
eps     = 1e-10

for var_name in variables:
    if var_name in results.variables:
        if var_name in Dataset(files[0], 'r').variables:
            fig = plt.figure()
            ax1    = fig.add_axes([0.15, 0.25, 0.8, 0.7])
            ax2    = fig.add_axes([0.15, 0.15, 0.8, 0.1])
            leg_ax = fig.add_axes([0.15, 0.05, 0.8, 0.1])

            if var_name == "energy_transfer":
                label = '$\int_{\Omega}-\Pi\mathrm{dA}$ $\quad\left(\mathrm{J}\\right)$'
                scale_factor = -1
            elif var_name == "Lambda_m":
                label = '$\int_{\Omega}\Lambda^m\mathrm{dA}$'
                scale_factor = 1
            elif var_name == "PEtoKE":
                label = '$\int_{\Omega}\overline{\\rho}g\overline{u}_r\mathrm{dA}$'
                scale_factor = 1
            elif var_name == "div_Jtransport":
                label = '$\int_{\Omega}-\\nabla\cdot J\mathrm{dA}$'
                scale_factor = -1

            label = '$\int^t$' + label

            # Loop through filters
            cnt = 0
            for fp in files:
                with Dataset(fp, 'r') as results:
        
                    scale = results.filter_scale
                    if (rank == 0):
                        print('  {0} - {1:.4g}km'.format(var_name,scale/1e3))
    
                    time  = results.variables['time'][:] * (60*60) # convert hours to seconds
                    Ntime = len(time)
            
                    if (Ntime == 1):
                        if (rank == 0):
                            print("Not enough time points to compare integrated fluxes")
                        sys.exit()
            
                    # Do some time handling to adjust the epochs
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
            
                    net      = np.zeros((Ntime,))
                    net_proc = np.zeros((Ntime,))
            
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
            
                        data         = results.variables[var_name][iT,:,:,:]
                        net_proc[iT] = np.sum( scale_factor * data * dArea )
                
                    comm.Allreduce(net_proc, net, op=MPI.SUM)
            
                    ##
                    ## Also plot the time-integrated version
                    ##
            
                    colours = plt.rcParams['axes.prop_cycle'].by_key()['color']
            
                    to_plot = spi.cumtrapz(net, time)
                    ax1.plot(time[1:], to_plot, '-', color=colours[cnt], 
                        label='{0:.3g}km'.format(scale/1e3))
                    #ax1.plot(np.ma.masked_where(to_plot>0, time[1:]),
                    #         np.ma.masked_where(to_plot>0, np.abs(to_plot)),
                    #        '-', color=colours[cnt], label='{0:.3g}km'.format(scale/1e3))
                    #ax1.plot(np.ma.masked_where(to_plot<0, time[1:]),
                    #         np.ma.masked_where(to_plot<0, np.abs(to_plot)),
                    #        '--', color=colours[cnt])
                cnt += 1
            
            #ax1.set_yscale('log')
            ax1.set_ylabel(label)
            ax1.plot(time, 0*time, '--k')
            
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
            
            plt.savefig(out_direct + '/Comparison_{0}.pdf'.format(var_name))
            plt.close()
