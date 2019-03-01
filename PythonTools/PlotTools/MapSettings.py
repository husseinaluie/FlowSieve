import numpy as np

def MapSettings(longitude, latitude, R_earth = 6371e3):

    # Resolution options are (increasing to the right): 'c', 'l' , 'i', 'h', 'f'

    # Some parameters for plotting
    lon_0 = np.mean(longitude)
    lat_0 = np.mean(latitude)

    # Set the projection
    if np.abs(longitude.max() - longitude.min()) > 75.:
        print("Using Eckert IV Projection")
        map_settings = dict(
            projection = 'eck4',
            lon_0 = -50,#lon_0, 
            resolution = 'l')  
    else:
        print("Using Gnomonic Projection")

        height = 1.2 * R_earth * (latitude.max()  - latitude.min() )
        width  = 1.2 * R_earth * (longitude.max() - longitude.min())

        ll_corner_lon = np.min(longitude)
        ll_corner_lat = np.min(latitude)

        ur_corner_lon = np.max(longitude)
        ur_corner_lat = np.max(latitude)

        map_settings = dict(
            projection = 'gnom',
            lon_0 = lon_0,
            lat_0 = lat_0,
            llcrnrlon = ll_corner_lon,
            llcrnrlat = ll_corner_lat - 5,
            urcrnrlon = ur_corner_lon + 15,
            urcrnrlat = ur_corner_lat,
            #height = height,
            #width  = width,
            resolution = 'l')  

    return map_settings
