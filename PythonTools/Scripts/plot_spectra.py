import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import cmocean, sys
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
import matpy as mp

fp = 'filter_output.nc'
results = Dataset(fp, 'r')

R_earth = 6371e3
rho0    = 1e3
D2R = np.pi / 180
R2D = 180 / np.pi
eps = 1e-10

# Get the grid from the first filter
latitude  = results.variables['latitude'][:] * R2D
longitude = results.variables['longitude'][:] * R2D
scales    = results.variables['scale'][:]
depth     = results.variables['depth'][:]
time      = results.variables['time'][:]
mask      = results.variables['mask'][:]
transfer  = results.variables['energy_transfer'][:, 0, 0, :, :] * 1e6  # unit conversion
u_r_set   = results.variables['u_r'  ][:, 0, 0, :, :]
u_lat_set = results.variables['u_lat'][:, 0, 0, :, :]
u_lon_set = results.variables['u_lon'][:, 0, 0, :, :]

L_lon = R_earth * (longitude.max() - longitude.min())
L_lat = R_earth * ( latitude.max() -  latitude.min())

LON, LAT = np.meshgrid(longitude * D2R, latitude * D2R)

dlat = (latitude[ 1] - latitude[ 0]) * np.pi / 180
dlon = (longitude[1] - longitude[0]) * np.pi / 180
dAreas = R_earth**2 * np.cos(LAT) * dlat * dlon

scales = scales[:-1]  # get rid of zero
num_scales = len(scales)

KE_tot = 0.5 * (   np.sum(u_r_set, axis=0)**2
                 + np.sum(u_lat_set, axis=0)**2
                 + np.sum(u_lon_set, axis=0)**2
                )
KE_tot = np.ma.masked_where(mask==0, KE_tot)
KE_tot = rho0 * np.ma.sum(KE_tot * dAreas * mask)

## First, compute domain-integrated KE
cumul_spec   = np.zeros(num_scales)
net_transfer = np.zeros(num_scales)

# Plot each band
for ii in range(num_scales):

    u_r   = u_r_set[  ii:,:,:]
    u_lat = u_lat_set[ii:,:,:]
    u_lon = u_lon_set[ii:,:,:]

    if (len(u_r.shape) == 3):
        u_r   = np.sum(u_r, axis=0)
        u_lat = np.sum(u_lat, axis=0)
        u_lon = np.sum(u_lon, axis=0)

    u_r   = np.ma.masked_where(mask==0, u_r)
    u_lat = np.ma.masked_where(mask==0, u_lat)
    u_lon = np.ma.masked_where(mask==0, u_lon)

    ener_transf = transfer[ii,:,:]
    ener_transf = np.ma.masked_where(mask==0, ener_transf)

    KE = 0.5 * ( u_r**2 + u_lat**2 + u_lon**2 )

    cumul_spec[ii] = rho0 * np.sum(KE * dAreas * mask)

    net_transfer[ii] = np.sum(transfer[ii,:,:] * dAreas * mask)

L = min(L_lon, L_lat)
k_l = L/scales
deriv_k_l     = mp.FiniteDiff(         k_l,  2, uniform=False)
deriv_k_l_log = mp.FiniteDiff(np.log10(k_l), 2, uniform=False)

spectrum        = deriv_k_l.dot(cumul_spec)
spectral_slopes = deriv_k_l_log.dot(np.log10(spectrum))

fig, axes = plt.subplots(4, 1, sharex=True, figsize=(5,8))

# For consistency with publications
horiz_scale = 1e4 / (scales / 1e3)

axes[0].plot(horiz_scale, cumul_spec / KE_tot, '-o')
axes[1].plot(horiz_scale, spectrum / KE_tot, '-o')
axes[2].plot(horiz_scale, net_transfer / KE_tot, '-o')
axes[3].plot(horiz_scale, spectral_slopes, '-o')
axes[3].set_xlabel('$10^4\cdot l^{-1} \,(\mathrm{km}^{-1})$')

#axes[2].set_ylabel('\Pi/\mathrm{KE}_{\mathrm{total}}\nW \cdot\mathrm{km}^{-2}\cdot\mathrm{m}^{-1}')
axes[2].plot([horiz_scale[0], horiz_scale[-1]], [0,0], '--k')
axes[3].plot([horiz_scale[0], horiz_scale[-1]], [0,0], '--k')

for ax in axes:
    #ax.xaxis.get_major_formatter().set_powerlimits((0, 4))
    ax.yaxis.get_major_formatter().set_powerlimits((0, 4))
    ax.set_xscale('log')

plt.tight_layout(True)
plt.savefig('Figures/spectra.png', dpi=500)
