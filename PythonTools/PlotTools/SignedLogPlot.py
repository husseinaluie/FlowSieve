from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import cmocean
import numpy as np

def SignedLogPlot(XX, YY, ZZ, ax, fig,
        num_ords=6, xscale = 'linear', yscale = 'linear',
        xlabel = '', ylabel = '', title = ''):
    # figsize:  size of the figure (passed to figsize flag in matplotlib)
    # num_ords: number of orders of magnitude to show on colour bar
    # dpi:      resolution of output file
    # xlabel, ylabel, title: corresponding matplotlib labels

    ## Extract nice colour maps from cmocean 'balance'
    Np = 256
    bal = cmocean.tools.get_dict(cmocean.cm.balance, N=Np)

    bal_red   = np.array([x[2] for x in bal['red']  ])
    bal_blue  = np.array([x[2] for x in bal['blue'] ])
    bal_green = np.array([x[2] for x in bal['green']])

    just_red  = np.vstack([ bal_red[ Np//2:],       bal_green[ Np//2:],       bal_blue[ Np//2:]       ]).T
    just_blue = np.vstack([ bal_red[:Np//2 ][::-1], bal_green[:Np//2 ][::-1], bal_blue[:Np//2 ][::-1] ]).T

    red_map  = cmocean.tools.cmap(just_red,  N = Np//2 )
    blue_map = cmocean.tools.cmap(just_blue, N = Np//2 )

    red_map.set_under( 'w')
    blue_map.set_under('w')

    ## Plot
    vmax = np.max(np.abs(ZZ))
    vmin = 10**(np.ceil(np.log10(vmax)-num_ords)-1.0)

    # Plot the positives in red
    q1 = ax.pcolormesh(XX, YY, np.abs(ZZ)*(ZZ > 0), cmap=red_map, \
            norm=LogNorm(vmin=vmin, vmax=vmax), linewidth=0, rasterized=True)

    # Plot the negatives in blue
    q2 = ax.pcolormesh(XX, YY, np.abs(ZZ)*(ZZ < 0), cmap=blue_map, \
            norm=LogNorm(vmin=vmin, vmax=vmax), linewidth=0, rasterized=True)

    ax.yaxis.get_major_formatter().set_powerlimits((0, 4))

    # Set axis labels
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)

    # Set axis scales
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)

    ax_pos = ax.get_position()
    x0 = ax_pos.x0
    y0 = ax_pos.y0
    width  = ax_pos.x1 - x0
    height = ax_pos.y1 - y0

    ax.set_position([x0, y0, 0.85*width, height])

    x0 += 0.9*width

    # Axes for colour bars
    cax1 = fig.add_axes([x0, y0 + height*0.525, 0.05*width, height*0.475])
    cax2 = fig.add_axes([x0, y0,               0.05*width, height*0.475])

    # Colour bars
    cbar1 = plt.colorbar(q1, cax=cax1, extend='min')
    cbar1.solids.set_rasterized(True)
    cbar1.solids.set_edgecolor("face")

    cbar2 = plt.colorbar(q2, cax=cax2, extend='min')
    cbar2.solids.set_rasterized(True)
    cbar2.solids.set_edgecolor("face")

    # Add negative sign to negative side
    #labels = [st[:14] + '-' + st[14:] for st in 
    #                [tick.get_text() for tick in 
    #                    cax2.get_yticklabels()]]
    #print(labels)
    #cax2.set_yticklabels(labels)

    cax2.invert_yaxis()

    return (q1, q2), (cbar1, cbar2)
