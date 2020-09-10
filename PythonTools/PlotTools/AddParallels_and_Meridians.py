import matplotlib.pyplot as plt
import numpy as np

def AddParallels_and_Meridians(
        ax, proj, 
        parallels, meridians, 
        latitude, longitude,
        label_parallels = True,
        label_meridians = True,
        colour = 'g',
        linestyle = '--',
        linewidth = 1):

    for mer in meridians:
        xs, ys = proj(mer*np.ones(latitude.shape), latitude, inverse=False)
        ax.plot(xs, ys, linestyle = linestyle, linewidth = linewidth, color = colour)
    for par in parallels:
        xs, ys = proj(longitude, par*np.ones(longitude.shape), inverse=False)
        ax.plot(xs, ys, linestyle = linestyle, linewidth = linewidth, color = colour)

    xs, ys = proj(meridians, latitude.min()*np.ones(meridians.shape), inverse=False)
    ax.set_xticks(xs)
    if label_meridians:
        ax.set_xticklabels(['{0:d}'.format(int(mer)) for mer in meridians])
    else:
        ax.set_xticklabels([])

    xs, ys = proj(longitude.min()*np.ones(parallels.shape), parallels, inverse=False)
    ax.set_yticks(ys)
    if label_parallels:
        ax.set_yticklabels(['{0:d}'.format(int(par)) for par in parallels])
    else:
        ax.set_yticklabels([])
