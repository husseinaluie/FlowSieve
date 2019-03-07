import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import cmocean, sys
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
import PlotTools

source = Dataset('input.nc', 'r')

fp = 'filter_output.nc'
results = Dataset(fp, 'r')

R_earth = 6371e3
D2R     = np.pi / 180
R2D     = 180 / np.pi
eps     = 1e-10

# Get the grid from the first filter
latitude  = results.variables['latitude'][:] * R2D
longitude = results.variables['longitude'][:] * R2D
scales    = results.variables['scale'][:]
depth     = results.variables['depth'][:]
time      = results.variables['time'][:]
mask      = results.variables['mask'][:]

num_scales = len(scales)-1

dlat = latitude[1]  - latitude[0]
dlon = longitude[1] - longitude[0]
LON, LAT = np.meshgrid(longitude * D2R, latitude * D2R)
dAreas = R_earth**2 * np.cos(LAT * np.pi / 180) * dlat * dlon

uo = source.variables['uo'][0, 0, :, :]
vo = source.variables['vo'][0, 0, :, :]
Full_KE = 0.5 * (uo**2 + vo**2)


# Some parameters for plotting
map_settings = PlotTools.MapSettings(longitude, latitude)

meridians = np.round(np.linspace(longitude.min(), longitude.max(), 5))
parallels = np.round(np.linspace(latitude.min(),  latitude.max(),  5))

cbar_props = dict(pad = 0.1, shrink = 0.85, orientation = 'vertical')
gridspec_props = dict(wspace = 0.05, hspace = 0.05, left = 0.02, right = 0.98, bottom = 0.02, top = 0.98)

    
## Full KE
plt.figure()
ax = plt.subplot(1,1,1)

to_plot = Full_KE.copy()
to_plot = np.ma.masked_where(mask==0, to_plot)

m  = Basemap(ax = ax, **map_settings)

CV = np.max(np.abs(to_plot))
qm = m.pcolormesh(LON*R2D, LAT*R2D, to_plot, cmap='cmo.amp', vmin = 0, vmax = CV, latlon = True)
    
cbar = plt.colorbar(qm, ax = ax, **cbar_props)
PlotTools.ScientificCbar(cbar, units='')

# Add coastlines and lat/lon lines
m.drawcoastlines(linewidth=0.1)
m.drawparallels(parallels, linewidth=0.5, labels=[0,1,0,0], color='g')
m.drawmeridians(meridians, linewidth=0.5, labels=[0,0,1,0], color='g')
m.contourf(LON*R2D, LAT*R2D, mask, [-0.5, 0.5], colors='gray', hatches=['','///\\\\\\'], latlon=True)

plt.savefig('Figures/KE_full.png', dpi=500)
plt.close()

    
## Full uo
plt.figure()
ax = plt.subplot(1,1,1)

to_plot = uo.copy()
to_plot = np.ma.masked_where(mask==0, to_plot)

m  = Basemap(ax = ax, **map_settings)

CV = np.max(np.abs(to_plot))
qm = m.pcolormesh(LON*R2D, LAT*R2D, to_plot, cmap='cmo.balance', vmin = -CV, vmax = CV, latlon = True)
    
cbar = plt.colorbar(qm, ax = ax, **cbar_props)
PlotTools.ScientificCbar(cbar, units='')

# Add coastlines and lat/lon lines
m.drawcoastlines(linewidth=0.1)
m.drawparallels(parallels, linewidth=0.5, labels=[0,1,0,0], color='g')
m.drawmeridians(meridians, linewidth=0.5, labels=[0,0,1,0], color='g')
m.contourf(LON*R2D, LAT*R2D, mask, [-0.5, 0.5], colors='gray', hatches=['','///\\\\\\'], latlon=True)

plt.savefig('Figures/uo_full.png', dpi=500)
plt.close()

    
## Full vo
plt.figure()
ax = plt.subplot(1,1,1)

to_plot = vo.copy()
to_plot = np.ma.masked_where(mask==0, to_plot)

m  = Basemap(ax = ax, **map_settings)

CV = np.max(np.abs(to_plot))
qm = m.pcolormesh(LON*R2D, LAT*R2D, to_plot, cmap='cmo.balance', vmin = -CV, vmax = CV, latlon = True)
    
cbar = plt.colorbar(qm, ax = ax, **cbar_props)
PlotTools.ScientificCbar(cbar, units='')

# Add coastlines and lat/lon lines
m.drawcoastlines(linewidth=0.1)
m.drawparallels(parallels, linewidth=0.5, labels=[0,1,0,0], color='g')
m.drawmeridians(meridians, linewidth=0.5, labels=[0,0,1,0], color='g')
m.contourf(LON*R2D, LAT*R2D, mask, [-0.5, 0.5], colors='gray', hatches=['','///\\\\\\'], latlon=True)

plt.savefig('Figures/vo_full.png', dpi=500)
plt.close()

