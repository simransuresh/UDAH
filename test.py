
# from compute_gsc import lats, lons, dot, ug, vg
# from scipy.interpolate import interp2d
# print(lats.shape, lons.shape, dot.shape, ug.shape, vg.shape)
# lons = lons+180

# f = interp2d(lons, lats, dot)
# interp_dot = f(longitudes, latitudes)
# print(interp_dot.shape)
# f = interp2d(lons, lats, ug)
# interp_nug = f(longitudes, latitudes)
# print(interp_nug.shape)
# f = interp2d(lons, lats, vg)
# interp_nvg = f(longitudes, latitudes)
# print(interp_nvg.shape)

# import numpy as np
# import matplotlib.pyplot as plt
# import datetime
# from mpl_toolkits.basemap import Basemap, shiftgrid
# from netCDF4 import Dataset

# # specify date to plot.
# yyyy=1993; mm=3; dd=14; hh=0
# date = datetime.datetime(yyyy,mm,dd,hh)
# # set OpenDAP server URL.
# URLbase="https://www.ncei.noaa.gov/thredds/dodsC/model-cfs_reanl_6h_pgb/"
# URL=URLbase+"%04i/%04i%02i/%04i%02i%02i/pgbh00.gdas.%04i%02i%02i%02i.grb2" %\
#              (yyyy,yyyy,mm,yyyy,mm,dd,yyyy,mm,dd,hh)
# data = Dataset(URL)

# latitudes = data.variables['lat'][::-1]
# longitudes = data.variables['lon'][:].tolist() 
# # latitudes = [lat for lat in latitudes if float(lat)>=60 and float(lat)<=88]
# print(len(latitudes), latitudes)
# print(len(longitudes), longitudes)

# slpin = 0.01*data.variables['Pressure_msl'][:].squeeze()
# uin = data.variables['u-component_of_wind_height_above_ground'][:].squeeze()
# vin = data.variables['v-component_of_wind_height_above_ground'][:].squeeze()
# print(slpin.shape, uin.shape, vin.shape)
# print(uin)
# print(vin)

# # add cyclic points manually (could use addcyclic function)
# slp = np.zeros((slpin.shape[0],slpin.shape[1]+1),np.float64)
# slp[:,0:-1] = slpin[::-1]; slp[:,-1] = slpin[::-1,0]
# u = np.zeros((uin.shape[0],uin.shape[1]+1),np.float64)
# u[:,0:-1] = uin[::-1]; u[:,-1] = uin[::-1,0]
# v = np.zeros((vin.shape[0],vin.shape[1]+1),np.float64)
# v[:,0:-1] = vin[::-1]; v[:,-1] = vin[::-1,0]
# longitudes.append(360.); longitudes = np.array(longitudes)
# lons, lats = np.meshgrid(longitudes,latitudes)

# m = Basemap(projection='nplaea',boundinglat=60,lon_0=0,resolution='l')
# fig1 = plt.figure(figsize=(6,6))
# ax = fig1.add_axes([0.1,0.1,0.8,0.8])
# x, y = m(lons, lats)

# ugrid,newlons = shiftgrid(180.,u,longitudes,start=False)
# vgrid,newlons = shiftgrid(180.,v,longitudes,start=False)

# uproj,vproj,xx,yy = \
# m.transform_vector(ugrid,vgrid,newlons,latitudes,31,31,returnxy=True,masked=True)
# Q = m.quiver(xx,yy,uproj,vproj,scale=700)
# qk = plt.quiverkey(Q, 0.1, 0.1, 20, '20 m/s', labelpos='W')
# m.drawcoastlines(linewidth=1.5)
# plt.show()

# # ############# plotting animation of annual mean
# fig, ax = plt.subplots(1, 2, figsize=(10, 8))

# m = Basemap(projection='npstere', lon_0=0, lat_0=90, boundinglat=60, ax=ax[0])
# x, y = m(lon, lat)
# def animate1(frame):
#     m.pcolormesh(x, y, dot_annual_mean[frame][:,:], cmap='viridis', ax=ax[0])
#     m.colorbar(ax=ax[0], label='DOT [m]')
#     ax[0].quiver(x, y, ug_annual_mean[frame][:,:], vg_annual_mean[frame][:,:], scale=10, color='black', alpha=0.5)
#     ax[0].set_title(frame)
#     return fig,
# animation1 = FuncAnimation(fig, animate1, frames=dot_annual_mean.keys(), interval=100) 


# ############# plotting animation of biannual mean
# m = Basemap(projection='npstere', lon_0=0, lat_0=90, boundinglat=60, ax=ax[1])
# x, y = m(lon, lat)
# def animate2(frame):
#     m.pcolormesh(x, y, dot_biannual_mean[frame][:,:], cmap='viridis', ax=ax[1])
#     m.colorbar(ax=ax[1], label='DOT [m]')
#     ax[1].quiver(x, y, ug_biannual_mean[frame][:,:], vg_biannual_mean[frame][:,:], scale=10, color='black', alpha=0.5)
#     ax[1].set_title(frame)
#     return fig,
# animation2 = FuncAnimation(fig, animate2, frames=dot_biannual_mean.keys(), interval=200) 

# m.drawcoastlines(ax=ax[0])
# m.drawcoastlines(ax=ax[1])
# ax[0].set_aspect('equal')
# ax[1].set_aspect('equal')
# # plt.tight_layout()
# plt.show()

import pandas as pd
from scipy.spatial import KDTree
import time

print(time.time())
# Load the CSV file into a DataFrame
data = pd.read_csv("data/landsea_04.csv")

# Extract coordinates and create a KDTree
coordinates = list(zip(data['Latitude'], data['Longitude']))
kdtree = KDTree(coordinates)

def is_land_or_ocean(lat, lon):
    _, idx = kdtree.query((lat, lon))
    return 'ocean' if int(data.iloc[idx]['Bottom_Standard_level'])==0 else 'land'

# Test the function
lat, lon = 60, 40.125  # Example input coordinates
result = is_land_or_ocean(lat, lon)
print(result)

# if result == 1:
#     print("The point is on land.")
# else:
#     print("The point is in the ocean.")

print(time.time())