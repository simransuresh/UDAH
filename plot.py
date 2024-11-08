import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from scipy.spatial import KDTree

###### NOTE plotting gridded UDAH data
# landmask_data = pd.read_csv("data/landsea_04.csv")
# landmask_coordinates = list(zip(landmask_data['Latitude'], landmask_data['Longitude']))
# kdtree = KDTree(landmask_coordinates)

# def is_land_or_ocean(lat, lon):
#     _, idx = kdtree.query((lat, lon))
#     depth = landmask_data.iloc[idx]['Bottom_Standard_level']
#     return 'land' if depth == 1 else 'ocean'

# main_data1 = pd.read_csv("results/grid_dh_20110901_1.csv")
# main_data1['Land_Or_Ocean'] = main_data1.apply(lambda row: is_land_or_ocean(row['Latitude'], row['Longitude']), axis=1)
# main_data2 = pd.read_csv("results/grid_dh_20110901_2.csv")
# main_data2['Land_Or_Ocean'] = main_data2.apply(lambda row: is_land_or_ocean(row['Latitude'], row['Longitude']), axis=1)

main_data1 = pd.read_csv("results/gridded_dh_2012_01.csv")
lat = main_data1['Latitude']
lon = main_data1['Longitude']
dh = main_data1['DH_error']
# lat = np.concatenate((main_data1['Latitude'], main_data2['Latitude'])).reshape(121, 480)[:113, :]
# lon = np.concatenate((main_data1['Longitude'], main_data2['Longitude'])).reshape(121, 480)[:113, :]
# dh = np.concatenate((main_data1['Dynamic_Height'], main_data2['Dynamic_Height'])).reshape(121, 480)[:113, :]

# from read_sla_currents import mdt, temp_mean_sla, sla, ug, vg
# dot = mdt + sla + temp_mean_sla
# dot = dot[12, :, :]  # sept 2011 month index in SAGA
# ug = ug[12, :, :]
# vg = vg[12, :, :]

# final = np.abs(dh - dot)
           
fig, ax = plt.subplots(figsize=(10, 10), subplot_kw={'projection': ccrs.NorthPolarStereo()})
# ax.set_extent([-160, -120, 65, 85], ccrs.PlateCarree())
ax.set_extent([-180, 180, 60, 90], ccrs.PlateCarree())

ax.coastlines()
ax.gridlines(draw_labels=True)

scatter = ax.scatter(lon, lat, c=dh, cmap='YlGnBu', s=10, transform=ccrs.PlateCarree(), vmin=0, vmax=1.2)
cbar = plt.colorbar(scatter, ax=ax, orientation='vertical', label='Difference in DH and DOT [m]')

# ax.set_title("Dynamic Height of Arctic ocean from hydrographic observations 2011-2018")
plt.show()
# plt.savefig('figures/dh_dot_diff.png')


########## NOTE plotting all dynamic height from UDAH
# import pandas as pd
# import matplotlib.pyplot as plt
# import cartopy.crs as ccrs
# import cartopy.feature as cfeature

# # Load the CSV file
# file_path = 'results/dh_2011_2018.csv'  # Update this path
# data = pd.read_csv(file_path)

# # Convert Datetime to datetime format
# data['Datetime'] = pd.to_datetime(data['Datetime'])

# # Filter data for the year 2011
# data_2011 = data[data['Datetime'].dt.year == 2011]

# # Plotting the dynamic height for 2011 in a North Polar Stereographic projection
# fig, ax = plt.subplots(figsize=(10, 10), subplot_kw={'projection': ccrs.NorthPolarStereo()})

# # Step 5: Set map boundaries and limits
# ax.set_extent([-180, 180, 60, 90], ccrs.PlateCarree())

# # Add coastlines, gridlines, and labels
# ax.coastlines()
# ax.gridlines(draw_labels=True)

# # Step 6: Plot the dynamic height data
# scatter = ax.scatter(data_2011['Longitude'], data_2011['Latitude'], c=data_2011['Dynamic_height'], cmap='YlGnBu', s=10, transform=ccrs.PlateCarree(), vmin=0, vmax=1.2)

# # Add a color bar to represent dynamic height values
# cbar = plt.colorbar(scatter, ax=ax, orientation='vertical', label='Dynamic Height [m]')

# # Title and display the plot
# ax.set_title("Dynamic Height of Arctic ocean 2011")
# plt.savefig('figures/dh_hydr_obs_2011.png')
# # plt.show()