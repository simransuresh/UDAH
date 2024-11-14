import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from scipy.spatial import KDTree

###### NOTE plotting gridded UDAH data
# main_data1 = pd.read_csv("results/grd_dh_2012_01.csv")
df = pd.read_csv("results/grd_dh_2011_01_final.csv")
dh = df['Dynamic_Height']

# for idx, row in  enumerate(main_data1['Latitude']):
#     if main_data1['Longitude'][idx]>-20 and main_data1['Longitude'][idx]<98 and main_data1['Latitude'][idx]<82: # remove Fram strait from mapping
#         main_data1['Dynamic_Height'][idx] = np.nan
#         main_data1['DH_error'][idx] = np.nan

# for idx, row in  enumerate(df['Latitude']):
#     if df['Longitude'][idx]>-20 and df['Longitude'][idx]<80 and df['Latitude'][idx]<82: # remove Fram strait from mapping
#         df['Depth'][idx] = np.nan
#     if df['Longitude'][idx]>-90 and df['Longitude'][idx]<0 and df['Latitude'][idx]<80:
#         df['Depth'][idx] = np.nan

# # df.to_csv('filt_depth.csv', index=False)
# # Drop rows with any NaN values and save to a new CSV file
# df.dropna().to_csv('data_points.csv', index=False)

# # Count rows where all columns have non-NaN values
# count_no_nan = (df.notna().all(axis=1)).sum()
# print(f"Number of rows with no NaN values: {count_no_nan}")

# lat = np.array(df['Latitude']).reshape(121, 480)[:113, :]
# lon = np.array(df['Longitude']).reshape(121, 480)[:113, :]
# dh = np.array(df['Dynamic_Height']).reshape(121, 480)[:113, :]
# dhe = np.array(df['DH_error']).reshape(121, 480)[:113, :]
   
# from read_sla_currents import mdt, temp_mean_sla, sla, ug, vg
# dot = mdt + sla + temp_mean_sla
# dot = dot[8, :, :]  # sept 2011 month index in SAGA
# ug = ug[8, :, :]
# vg = vg[8, :, :]   # TODO compare by populating with np.nan mat and filling with val at filt grd points

# final = np.abs(dh - dot)
                   
print(np.nanmin(dh), np.nanmax(dh))

fig, ax = plt.subplots(figsize=(10, 10), subplot_kw={'projection': ccrs.NorthPolarStereo()})
# ax.set_extent([-160, -120, 65, 85], ccrs.PlateCarree())
ax.set_extent([-180, 180, 70, 90], ccrs.PlateCarree())

ax.coastlines()
ax.gridlines(draw_labels=True)

# scatter = ax.scatter(lon, lat, c=final, cmap='YlGnBu', s=10, transform=ccrs.PlateCarree(), vmin=0, vmax=1.2)
scatter = ax.scatter(df['Longitude'], df['Latitude'], c=dh, cmap='gist_rainbow', s=10, transform=ccrs.PlateCarree(), vmin=0, vmax=1)
cbar = plt.colorbar(scatter, ax=ax, orientation='vertical', label='Difference in DH and DOT [m]')

# ax.set_title("Dynamic Height of Arctic ocean from hydrographic observations 2011-2018")
plt.show()
# plt.savefig('figures/grd_dh_201101_final.png')


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