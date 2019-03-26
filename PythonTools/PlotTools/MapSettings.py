import numpy as np
from mpl_toolkits.basemap import Basemap

def MapSettings(longitude, latitude, R_earth = 6371e3):

    # Resolution options are (increasing to the right): 'c', 'l' , 'i', 'h', 'f'

    # Some parameters for plotting
    lon_0 = np.mean(longitude)
    lat_0 = np.mean(latitude)

    # Set the projection
    if np.abs(longitude.max() - longitude.min()) > 0.:
        print("Using Eckert IV Projection")
        map_settings = dict(
            projection = 'eck4',
            lon_0 = lon_0,
            resolution = 'l')  
    else:
        print("Using Gnomonic Projection")

        tmp_map = Basemap(projection = 'gnom',
            lon_0 = lon_0,
            lat_0 = lat_0,
            width = 20,
            height = 20,
            resolution = 'c')  

        LON, LAT = np.meshgrid(longitude, latitude)
        X, Y = tmp_map(LON, LAT, inverse=False)

        ll_lon, ll_lat = tmp_map(X.min(), Y.min(), inverse=True)
        ur_lon, ur_lat = tmp_map(X.max(), Y.max(), inverse=True)

        map_settings = dict(
            projection = 'gnom',
            lon_0 = lon_0,
            lat_0 = lat_0,
            llcrnrlon = ll_lon,
            llcrnrlat = ll_lat,
            urcrnrlon = ur_lon,
            urcrnrlat = ur_lat,
            resolution = 'l')  

    return map_settings
