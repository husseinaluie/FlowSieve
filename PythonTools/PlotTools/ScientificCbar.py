import matplotlib as mpl
import numpy as np

def ScientificCbar(cbar, units='',
        orientation='vertical', centre=True):

    # If requested, centre the colour bar
    if centre:
        cv = np.max(np.abs(cbar.mappable.get_clim()))
        cbar.mappable.set_clim(-cv, cv)

    # Limit the number of ticks on the colour bar
    tick_locator = mpl.ticker.MaxNLocator(nbins=5)
    cbar.locator = tick_locator
    cbar.update_ticks()
    ticks = cbar.get_ticks()

    # Re-scale the values to avoid the poorly-placed exponent
    #   above the colour bar
    scale = np.log10(np.max(np.abs(ticks)))
    scale = np.floor(scale)

    # Instead, simply add a ylabel to the colour bar giving the scale.
    if orientation == 'vertical':
        if scale != 0.:
            cbar.ax.set_yticklabels(["{0:.2g}".format(tick/(10**scale)) for tick in ticks])
        cbar.ax.set_ylabel('$\\times10^{' + '{0:d}'.format(int(scale)) + '}$ ' + units,
                rotation = '-90', labelpad=10)
    elif orientation == 'horizontal':
        if scale != 0.:
            cbar.ax.set_xticklabels(["{0:.2g}".format(tick/(10**scale)) for tick in ticks])
        cbar.ax.set_xlabel('$\\times10^{' + '{0:d}'.format(int(scale)) + '}$ ' + units,
                rotation = '0', labelpad=10)
