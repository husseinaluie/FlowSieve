import numpy as np
from pyproj import Proj
from netCDF4 import Dataset

def null_proj(x, y, inverse=False):
    return x, y

def MapProjection(longitude, latitude, R_earth = 6371e3):

    lon_0 = np.mean(longitude)
    lat_0 = np.mean(latitude)

    source = Dataset('input.nc', 'r')

    if source.variables['latitude'].units == 'm':
        proj = null_proj
    else:
        proj = Proj(proj='wag7', lon_0=lon_0, lat_0=lat_0)
        #proj = Proj(proj='eck6', lon_0=lon_0, lat_0=lat_0)

    return proj
