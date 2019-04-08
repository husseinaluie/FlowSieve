import numpy as np

def getAreas(source):
    
    D2R = np.pi / 180

    R_earth = 6371e3

    latitude  = source.variables['latitude' ][:] * D2R
    longitude = source.variables['longitude'][:] * D2R

    dlat = (latitude[1]  - latitude[0] )
    dlon = (longitude[1] - longitude[0])
    LON, LAT = np.meshgrid(longitude, latitude)

    try:
        units = source.variables['latitude'].units
    except:
        units = ''

    if units == 'm':
        dAreas = dlat * dlon * np.ones(LAT.shape)
    else:
        dAreas = R_earth**2 * np.cos(LAT) * dlat * dlon

    return dAreas
