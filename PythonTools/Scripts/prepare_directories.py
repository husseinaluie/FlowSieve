import numpy as np
import cmocean, sys, datetime
from netCDF4 import Dataset
import PlotTools, subprocess, shutil, os, glob

# If the output directory doesn't exist, create it.
#    Same with the subdirectories for each scale
out_direct = os.getcwd() + '/Videos'
if not(os.path.exists(out_direct)):
    os.makedirs(out_direct)

# Get the available filter files
files = glob.glob('filter_*.nc')

# Loop through filters
for fp in files:
    with Dataset(fp, 'r') as results:
        scale = results.filter_scale
        direct = out_direct + '/{0:.4g}km'.format(scale/1e3)
        if not(os.path.exists(direct)):
            os.makedirs(direct)
