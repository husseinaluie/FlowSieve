import numpy as np
import cmocean, sys, datetime
from netCDF4 import Dataset
import PlotTools, subprocess, shutil, os

source = Dataset('input.nc', 'r')

fp = 'filter_output.nc'
results = Dataset(fp, 'r')

# Get the grid from the first filter
scales    = results.variables['scale'][:]

# If the output directory doesn't exist, create it.
#    Same with the subdirectories for each scale
out_direct = os.getcwd() + '/Videos'
if not(os.path.exists(out_direct)):
    os.makedirs(out_direct)

for scale in scales[:-1]:
    direct = out_direct + '/{0:.4g}km'.format(scale/1e3)
    if not(os.path.exists(direct)):
        os.makedirs(direct)
