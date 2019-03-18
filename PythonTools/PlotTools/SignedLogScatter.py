import matplotlib.pyplot as plt
import cmocean
import numpy as np

def SignedLogScatter(
        XX, YY, CC, axes, s = 0.1, cmap = 'cmo.thermal',
        force_equal = False,
        scatter_kws = None, cbar_label=''):

    x_cv = np.max(np.abs(XX))
    y_cv = np.max(np.abs(YY))

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

    ax_pp = axes[0,1]
    ax_nn = axes[1,0]
    ax_np = axes[0,0]
    ax_pn = axes[1,1]

    # First, plot the all positives
    ax = ax_pp
    sel =  (XX>0) * (YY>0)
    x_sub = XX[sel == 1]
    y_sub = YY[sel == 1]
    c_sub = CC[sel == 1]

    ax.scatter(x_sub, y_sub, s = s, c = c_sub, **scatter_kws)

    # Next, plot the all negatives
    ax = ax_nn
    sel =  (XX<0) * (YY<0)
    x_sub = -XX[sel == 1]
    y_sub = -YY[sel == 1]
    c_sub =  CC[sel == 1]

    ax.scatter(x_sub, y_sub, s = s, c = c_sub, **scatter_kws)

    # Next, plot X > 0, Y < 0
    ax = ax_pn
    sel =  (XX>0) * (YY<0)
    x_sub =  XX[sel == 1]
    y_sub = -YY[sel == 1]
    c_sub =  CC[sel == 1]

    ax.scatter(x_sub, y_sub, s = s, c = c_sub, **scatter_kws)

    # Next, plot X < 0, Y > 0
    ax = ax_np
    sel =  (XX<0) * (YY>0)
    x_sub = -XX[sel == 1]
    y_sub =  YY[sel == 1]
    c_sub =  CC[sel == 1]

    sc = ax.scatter(x_sub, y_sub, s = s, c = c_sub, **scatter_kws)

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

    # Add colour bar
    cbar = plt.colorbar(sc, ax=axes)
    cbar.ax.set_ylabel(cbar_label)

