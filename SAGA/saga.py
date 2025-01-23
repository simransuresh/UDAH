import netCDF4 as nc
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import cartopy.crs as ccrs
import numpy as np
# from helpers import *
from datetime import date, timedelta

def convert_time(time):
    reference_date = date(2011, 1, 1)
    # time = days since 1st January 0000 at 00:00:00
    offset_days = [day - time[0] for day in time]
    # from ref, compute days in y-m format
    time = [ (reference_date+timedelta(days=day) ).strftime('%Y-%m') for day in offset_days]
    return time

####### monthly sla data
file_path = 'SAGA__SSHA_GEOVEL__60N_88N__2011_2020_rev_ext_alongtrack.nc'
ncfile = nc.Dataset(file_path, 'r')

# days since 1st January 0000 at 00:00:00
time = ncfile.variables['time'][:]
# change this to readable year-month format
time = convert_time(time)
# print(time)

lats = ncfile.variables['lat'][:][:,0]
lons = ncfile.variables['lon'][:][0,:]

mdt = ncfile.variables['mdt20112020'][:, :]
sla = ncfile.variables['sla'][:, :, :]
# dot is sum of constant mdt and changing sla
dot = mdt+sla

# using geostrophic equations - g,f, and grad of dot wrt x,y
ug = ncfile.variables['ug'][:, :, :]
vg = ncfile.variables['vg'][:, :, :]


###### plot of dot, ug, vg mean over whole period 2011-2020
# dot = np.mean(dot, axis=0)
# ug = np.mean(ug, axis=0)
# vg = np.mean(vg, axis=0)

# # Form Lambert equal spacing grid using lats lons
# m = Basemap(projection='nplaea', boundinglat=60, lon_0=0, resolution='l', round=True, llcrnrlat=60, urcrnrlat=90,
#             llcrnrlon=-180, urcrnrlon=180)

# fig1 = plt.figure(figsize=(6, 6))
# ax = fig1.add_axes([0.1,0.1,0.8,0.8])
# m.drawcoastlines(linewidth=1.5)
# m.drawparallels(range(60, 91, 10), labels=[1, 1, 1])
# m.drawmeridians(range(-180, 181, 60), labels=[1, 1, 1, 1])

# longitudes, latitudes = np.meshgrid(lons, lats)
# im = m.pcolormesh(longitudes,latitudes,dot,cmap=plt.cm.YlGnBu, latlon=True)
# cbar = plt.colorbar(im, ax=ax, orientation='vertical', shrink=0.7)
# cbar.set_label('DOT (m)')
# im.set_clim(0, 0.8)

# # value 32 is nicer with spacing is 212km, if set to 133, spacing is close to 50km
# uproj, vproj, X, Y = m.transform_vector(ug, vg, lons, lats, 50, 50, returnxy=True, masked=True)
# print(lons.shape, lats.shape, ug.shape, vg.shape)
# print(X.shape, Y.shape, uproj.shape, vproj.shape)
# print('Spacing:', np.diff(X, axis=1).mean(), np.diff(Y, axis=0).mean())

# # calculate scale based on ug vg values for quiver
# scale = np.nanmax(np.sqrt(uproj**2 + vproj**2)) / 0.1
# print(scale)

# Q = m.quiver(X, Y, uproj, vproj, scale=scale, width=0.004)
# qk = plt.quiverkey(Q, 0.1, 0.1, 0.1, '10 cm/s', labelpos='W')

# plt.show()
# # plt.savefig('saga_2012_annmean_50km.png', dpi=300)
