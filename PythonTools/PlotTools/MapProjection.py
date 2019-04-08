import numpy as np
from pyproj import Proj
from netCDF4 import Dataset

def null_proj(x, y, inverse=False):
    coef = ( np.pi / 180 ) / 1e3
    return x * coef, y * coef

def MapProjection(longitude, latitude, R_earth = 6371e3):

    lon_0 = np.mean(longitude)
    lat_0 = np.mean(latitude)

    source = Dataset('input.nc', 'r')
    try:
        unit = source.variables['latitude'].units
    except:
        unit = ''

    if unit == 'm':
        proj = null_proj
    else:
        proj = Proj(proj='wag7', lon_0=lon_0, lat_0=lat_0)
        #proj = Proj(proj='eck6', lon_0=lon_0, lat_0=lat_0)

    return proj
