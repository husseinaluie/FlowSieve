import matplotlib.pyplot as plt
import cmocean
import numpy as np

def SignedLogScatter(
        XX, YY, CC, axes, s = 0.2, cmap = 'cmo.thermal',
        num_ords_x = 12, num_ords_y = 12, force_equal = False,
        scatter_kws = None):

    x_cv = np.max(np.abs(XX))
    y_cv = np.max(np.abs(YY))

    if force_equal:
        num_ords = max(num_ords_x, num_ords_y)
        x_ub = max(x_cv, y_cv)
        y_ub = max(x_cv, y_cv)

        x_lb = 10**(np.log10(min(x_cv, y_cv)) - num_ords)
        y_lb = 10**(np.log10(min(x_cv, y_cv)) - num_ords)
    else:
        x_ub = x_cv
        y_ub = y_cv

        x_lb = 10**(np.log10(x_cv) - num_ords_x)
        y_lb = 10**(np.log10(y_cv) - num_ords_y)

    c_min = np.min(CC)
    c_max = np.max(CC)

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

    ax.scatter(x_sub, y_sub, s = s, c = c_sub, vmin=c_min, vmax=c_max, **scatter_kws)
    ax.set_xlim(x_lb, x_ub)
    ax.set_ylim(y_lb, y_ub)
    ax.set_xscale('log')
    ax.set_yscale('log')

    # Next, plot the all negatives
    ax = ax_nn
    sel =  (XX<0) * (YY<0)
    x_sub = -XX[sel == 1]
    y_sub = -YY[sel == 1]
    c_sub =  CC[sel == 1]

    ax.scatter(x_sub, y_sub, s = s, c = c_sub, vmin=c_min, vmax=c_max, **scatter_kws)
    ax.set_xlim(x_lb, x_ub)
    ax.set_ylim(y_lb, y_ub)
    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.invert_xaxis()
    ax.invert_yaxis()

    # Next, plot X > 0, Y < 0
    ax = ax_pn
    sel =  (XX>0) * (YY<0)
    x_sub =  XX[sel == 1]
    y_sub = -YY[sel == 1]
    c_sub =  CC[sel == 1]

    ax.scatter(x_sub, y_sub, s = s, c = c_sub, vmin=c_min, vmax=c_max, **scatter_kws)
    ax.set_xlim(x_lb, x_ub)
    ax.set_ylim(y_lb, y_ub)
    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.invert_yaxis()

    # Next, plot X < 0, Y > 0
    ax = ax_np
    sel =  (XX<0) * (YY>0)
    x_sub = -XX[sel == 1]
    y_sub =  YY[sel == 1]
    c_sub =  CC[sel == 1]

    sc = ax.scatter(x_sub, y_sub, s = s, c = c_sub, vmin=c_min, vmax=c_max, **scatter_kws)
    ax.set_xlim(x_lb, x_ub)
    ax.set_ylim(y_lb, y_ub)
    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.invert_xaxis()

    # Add negative signs to ax_nn ticklabels
    #labels = [st[:14] + '-' + st[14:] for st in 
    #                [tick.get_text() for tick in 
    #                    ax_nn.get_yticklabels()]]
    #print(labels)
    #ax_nn.set_yticklabels(labels)
    #labels = [st[:14] + '-' + st[14:] for st in 
    #                [tick.get_text() for tick in 
    #                    ax_nn.get_xticklabels()]]
    #print(labels)
    #ax_nn.set_xticklabels(labels)


    # Add colour bar
    plt.colorbar(sc, ax=axes)

