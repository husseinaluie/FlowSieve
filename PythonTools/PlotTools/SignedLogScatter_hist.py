import matplotlib.pyplot as plt
import cmocean
import numpy as np
from matplotlib.colors import LogNorm

def SignedLogScatter_hist(
        XX, YY, axes, s = 0.1, cmap = 'cmo.amp',
        force_equal = False, nbins_x = 100, nbins_y = 100, 
        log_cbar = True):

    x_max = np.nanmax(np.abs(XX))
    y_max = np.nanmax(np.abs(YY))

    #x_min = np.nanmin(np.abs(XX)[np.abs(XX)>0])
    #y_min = np.nanmin(np.abs(YY)[np.abs(YY)>0])
    x_min = np.percentile(np.abs(XX)[np.abs(XX)>0], 0.01)
    y_min = np.percentile(np.abs(YY)[np.abs(YY)>0], 0.01)

    if force_equal:
        x_ub = max(x_max, y_max)
        y_ub = max(x_max, y_max)

        x_lb = min(x_min, y_min)
        y_lb = min(x_min, y_min)
    else:
        x_ub = x_max
        y_ub = y_max

        x_lb = x_min
        y_lb = y_min

    x_edges = np.logspace(np.log10(x_lb), np.log10(x_ub), nbins_x+1)
    y_edges = np.logspace(np.log10(y_lb), np.log10(y_ub), nbins_y+1)

    xe_min = np.log10(x_edges.min())
    xe_max = np.log10(x_edges.max())
    xe_del = np.log10(x_edges[1]) - np.log10(x_edges[0])

    ye_min = np.log10(y_edges.min())
    ye_max = np.log10(y_edges.max())
    ye_del = np.log10(y_edges[1]) - np.log10(y_edges[0])

    dA = xe_del * ye_del

    ax_pp = axes[0,1]
    ax_nn = axes[1,0]
    ax_np = axes[0,0]
    ax_pn = axes[1,1]

    vmax = 0.
    vmin = np.inf

    # First, plot the all positives
    ax = ax_pp
    sel =  (XX>0) * (YY>0)
    x_sub = XX[sel == 1]
    y_sub = YY[sel == 1]

    counts = np.zeros((nbins_y, nbins_x))
    for II in range(len(x_sub)):
        Ix = int(np.floor((np.log10(x_sub[II]) - xe_min) / xe_del))
        Iy = int(np.floor((np.log10(y_sub[II]) - ye_min) / ye_del))
        if Ix >= nbins_x:
            Ix = nbins_x-1
        if Iy >= nbins_y:
            Iy = nbins_y-1
        counts[Iy,Ix] += 1
    num_counts = np.sum(counts)
    counts = counts / num_counts / dA

    vmax = max(vmax, counts.max())
    vmin = min(vmin, counts.min())

    qm_pp = ax.pcolormesh(x_edges, y_edges, counts, cmap=cmap)

    # Next, plot the all negatives
    ax = ax_nn
    sel =  (XX<0) * (YY<0)
    x_sub = -XX[sel == 1]
    y_sub = -YY[sel == 1]

    counts = np.zeros((nbins_y, nbins_x))
    for II in range(len(x_sub)):
        Ix = int(np.floor((np.log10(x_sub[II]) - xe_min) / xe_del))
        Iy = int(np.floor((np.log10(y_sub[II]) - ye_min) / ye_del))
        if Ix >= nbins_x:
            Ix = nbins_x-1
        if Iy >= nbins_y:
            Iy = nbins_y-1
        counts[Iy,Ix] += 1
    num_counts = np.sum(counts)
    counts = counts / num_counts / dA

    vmax = max(vmax, counts.max())
    vmin = min(vmin, counts.min())

    qm_nn = ax.pcolormesh(x_edges, y_edges, counts, cmap=cmap)

    # Next, plot X > 0, Y < 0
    ax = ax_pn
    sel =  (XX>0) * (YY<0)
    x_sub =  XX[sel == 1]
    y_sub = -YY[sel == 1]

    counts = np.zeros((nbins_y, nbins_x))
    for II in range(len(x_sub)):
        Ix = int(np.floor((np.log10(x_sub[II]) - xe_min) / xe_del))
        Iy = int(np.floor((np.log10(y_sub[II]) - ye_min) / ye_del))
        if Ix >= nbins_x:
            Ix = nbins_x-1
        if Iy >= nbins_y:
            Iy = nbins_y-1
        counts[Iy,Ix] += 1
    num_counts = np.sum(counts)
    counts = counts / num_counts / dA

    vmax = max(vmax, counts.max())
    vmin = min(vmin, counts.min())

    qm_pn = ax.pcolormesh(x_edges, y_edges, counts, cmap=cmap)

    # Next, plot X < 0, Y > 0
    ax = ax_np
    sel =  (XX<0) * (YY>0)
    x_sub = -XX[sel == 1]
    y_sub =  YY[sel == 1]

    counts = np.zeros((nbins_y, nbins_x))
    for II in range(len(x_sub)):
        Ix = int(np.floor((np.log10(x_sub[II]) - xe_min) / xe_del))
        Iy = int(np.floor((np.log10(y_sub[II]) - ye_min) / ye_del))
        if Ix >= nbins_x:
            Ix = nbins_x-1
        if Iy >= nbins_y:
            Iy = nbins_y-1
        counts[Iy,Ix] += 1
    num_counts = np.sum(counts)
    counts = counts / num_counts / dA

    vmax = max(vmax, counts.max())
    vmin = min(vmin, counts.min())

    qm_np = ax.pcolormesh(x_edges, y_edges, counts, cmap=cmap)

    for ax in [ax_pp, ax_np, ax_pn, ax_nn]:
        ax.set_xlim(x_lb, x_ub)
        ax.set_ylim(y_lb, y_ub)
        ax.set_xscale('log')
        ax.set_yscale('log')
    for ax in [ax_np, ax_nn]:
        ax.invert_xaxis()
    for ax in [ax_pn, ax_nn]:
        ax.invert_yaxis()

    # Explicitely add ticks. 
    #  Unfortunately, with order of rendering with log scales,
    #  this does seem to be necessary
    x_ub_exp = int(np.log10(x_ub))
    x_lb_exp = int(np.log10(x_lb))
    x_step   = int(np.ceil((x_ub_exp - x_lb_exp) / 5))
    xticks   = [np.power(10.,exp) for exp in np.arange(x_ub_exp, x_lb_exp, -x_step)][::-1]
    for ax in [ax_pp, ax_np, ax_pn, ax_nn]:
        ax.set_xticks(xticks)

    y_ub_exp = int(np.log10(y_ub))
    y_lb_exp = int(np.log10(y_lb))
    y_step   = int(np.ceil((y_ub_exp - y_lb_exp) / 5))
    yticks   = [np.power(10.,exp) for exp in np.arange(y_ub_exp, y_lb_exp, -y_step)][::-1]
    for ax in [ax_pp, ax_np, ax_pn, ax_nn]:
        ax.set_yticks(yticks)

    # Add negative signs to ax_nn ticklabels
    labels = ['$-10^{' + '{0:d}'.format(int(np.log10(tickval))) + '}$' for tickval in ax_nn.get_yticks()]
    ax_nn.set_yticklabels(labels)

    labels = ['$-10^{' + '{0:d}'.format(int(np.log10(tickval))) + '}$' for tickval in ax_nn.get_xticks()]
    ax_nn.set_xticklabels(labels)

    # Update norms / vmin/vmax
    for qm in [qm_pp, qm_np, qm_pn, qm_nn]:
        if log_cbar:
            vmin = vmax / 1e3
            qm.set_norm(LogNorm(vmin=vmin, vmax=vmax))
        else:
            qm.set_clim(vmin, vmax)


    # Add colour bar
    cbar = plt.colorbar(qm_np, ax=axes)
    cbar.ax.set_ylabel('Count Density')

