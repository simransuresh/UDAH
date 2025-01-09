import numpy as np
import pandas as pd
# from helpers import *
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

# Example DataFrame
# df = pd.read_csv('grd_dh_2012_01.csv')
# df['ug'] = np.full(len(df), np.nan)
# df['vg'] = np.full(len(df), np.nan)

###### checking nearest neigh function
# target_lat = 77.81912879787636
# target_lon = -141.38977105221358
# lat = round(target_lat, 3)
# lon = round(target_lon,3)

# Find the 4 nearest coordinates within the 49-51 km range
# nearest_coords = find_nearest_coords(target_lat, target_lon, df)
# print(nearest_coords)

# Print results
# idx = df.loc[(df['Latitude']==target_lat ) & (df['Longitude']==target_lon)].index[0]
# print(idx, df['X_meters'][idx], df['Y_meters'][idx], df['Surf_DH'][idx])

# nidx = coord[0]
# print(nidx, df['X_meters'][nidx], df['Y_meters'][nidx], coord[3], df['Surf_DH'][nidx])

# dx_plus = df.loc[(df['X_meters']==df['X_meters'][idx]+50000 ) & (df['Y_meters']==df['Y_meters'][idx])].index[0]
# dh_plus = df['Surf_DH'][dx_plus]
# dx_minus = df.loc[(df['X_meters']==df['X_meters'][idx]-50000 ) & (df['Y_meters']==df['Y_meters'][idx])].index[0]
# dh_minus = df['Surf_DH'][dx_minus]

# grad_x = (dh_plus-dh_minus)/(2*50000)

# dy_plus = df.loc[(df['X_meters']==df['X_meters'][idx] ) & (df['Y_meters']==df['Y_meters'][idx]+50000)].index[0]
# dh_plus = df['Surf_DH'][dy_plus]
# dy_minus = df.loc[(df['X_meters']==df['X_meters'][idx] ) & (df['Y_meters']==df['Y_meters'][idx]-50000)].index[0]
# dh_minus = df['Surf_DH'][dy_minus]

# grad_y = (dh_plus-dh_minus)/(2*50000)

# g = get_g(target_lat)
# f = coriolis(target_lat)
# print(grad_x, grad_y)

# ug = -g*grad_y/f
# vg = g*grad_x/f
# print(ug, vg)
    

##### combine different mapping vars to single file
# df = pd.read_csv('grd_dh_2012_01.csv')

# hfw_df = pd.read_csv('grd_dh_2012_01_hfw.csv')
# df['hFW'] = hfw_df['hFW']
# df['hFW_err'] = hfw_df['hFW_err']

# disoh_df = pd.read_csv('grd_dh_2012_01_disoh.csv')
# df['D_Siso'] = disoh_df['D_Siso']
# df['D_Siso_err'] = disoh_df['D_Siso_err']

# # print(df.tail)
# df.to_csv('grd_dh_2012_01.csv', index=False)


###### combine all months in a year to one - 2012
import os
# print(os.listdir('./'))
csv_files = [file for file in os.listdir('./') if file.startswith('grd_dh_2015_') and file not in ['data_500m.csv', 'grid_50km_nplaea.csv']]
print(csv_files)

# Initialize an empty list to hold DataFrames
dataframes = []

# Loop through the CSV files and read each into a DataFrame
for csv_file in csv_files:
    file_path = os.path.join('./', csv_file)
    df = pd.read_csv(file_path)
    dataframes.append(df)

merged_df = pd.concat(dataframes, ignore_index=True)
merged_df = merged_df.dropna(subset=['ug', 'vg'])
print(merged_df.head)

output_file = "grd_dh_2015.csv"
merged_df.to_csv(output_file, index=False)

# print(f"All CSV files merged successfully into {output_file}")

# import pandas as pd
# import numpy as np
# import matplotlib.pyplot as plt
# import cartopy.crs as ccrs
# import cartopy.feature as cfeature

# # Load the CSV file
# file_path = 'grd_2011_2018.csv'  # Replace with your actual file path
# data = pd.read_csv(file_path)
# # Ensure the 'Datetime' column is parsed as a datetime object
# data['Datetime'] = pd.to_datetime(data['Datetime'])

# # Extract the year and filter rows for 2011, 2012, 2013
# filtered_data = data[data['Datetime'].dt.year.isin([2014, 2015, 2016, 2017, 2018])]
# # print(filtered_data.head)
# filtered_data.to_csv(file_path)

# # Filter data for the specified region
# lat_min, lat_max = 70.0, 82.0
# lon_min, lon_max = -170.0, -120.0
# region_data = data[(data['Latitude'] >= lat_min) & (data['Latitude'] <= lat_max) &
#                    (data['Longitude'] >= lon_min) & (data['Longitude'] <= lon_max)]
# region_data['Datetime'] = pd.to_datetime(region_data['Datetime'])
# region_data['Year'] = region_data['Datetime'].dt.year

# # Calculate global min and max DH for all years to set fixed colorbar limits
# dh_vmin = region_data['Surf_DH'].min()
# dh_vmax = region_data['Surf_DH'].max()

# # Create figure and axes for 8 years (1 row, 8 columns)
# fig, axes = plt.subplots(1, 8, figsize=(36, 6), subplot_kw={'projection': ccrs.NorthPolarStereo()})
# axes = axes.flatten()

# # Loop through each year (2011 to 2018)
# for i, year in enumerate(range(2011, 2019)):
#     # Filter data for the current year
#     region_data_year = region_data[region_data['Year'] == year]
#     region_mean = region_data_year.groupby(['Latitude', 'Longitude'])[['Surf_DH', 'ug', 'vg']].mean().reset_index()

#     # Extract variables
#     lons = region_mean['Longitude'].values
#     lats = region_mean['Latitude'].values
#     dhs = region_mean['Surf_DH'].values
#     ug = region_mean['ug'].values
#     vg = region_mean['vg'].values

#     # Current axis
#     ax = axes[i]
#     ax.set_extent([-180, -120, 72, 82], crs=ccrs.PlateCarree())

#     # Add map features
#     ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
#     ax.add_feature(cfeature.BORDERS, linestyle=':', linewidth=0.5)
#     ax.add_feature(cfeature.LAND, color='lightgrey')

#     # Scatter plot for Dynamic Height (DH) with fixed colorbar limits
#     sc = ax.scatter(lons, lats, c=dhs, cmap='YlGnBu', s=130, transform=ccrs.PlateCarree())

#     # Quiver plot for geostrophic currents
#     scale = np.nanmax(np.sqrt(ug**2 + vg**2)) / 0.1
#     quiver = ax.quiver(lons, lats, ug, vg, scale=scale, width=0.003, color='black', transform=ccrs.PlateCarree())

#     # Title for each subplot
#     ax.set_title(f'{year}', fontsize=10)

# print(dh_vmin, dh_vmax)

# # Add a single colorbar outside the loop
# # cbar_ax = fig.add_axes([0.92, 0.2, 0.01, 0.6])  # Position for the colorbar
# cbar_ax = fig.add_axes([0.4, 0.08, 0.2, 0.03]) # Position for the colorbar
# cbar = plt.colorbar(sc, cax=cbar_ax, orientation='horizontal')
# cbar.set_label('Dynamic Height (m)', fontsize=10)
# # cbar.set_ticks(np.linspace(0.3, 1, 5))  # Optional: Set specific tick values

# # Super title for the entire figure
# # fig.suptitle('Dynamic Height and Geostrophic Currents (70°-82°N, -170° to -120°E)', fontsize=14)

# # Adjust spacing between subplots
# # plt.subplots_adjust(wspace=0.1)
# plt.tight_layout()
# # Show the plot
# plt.savefig('annual_bg.png', dpi=300)
