import cartopy.crs as ccrs
from mpl_toolkits.basemap import Basemap
import cartopy.feature as cfeature
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
# from bathymetric_depth import *
from ao import *

# file_path = 'grd_2011_2018.csv'  # Replace with your actual file path
file_path = 'grd_dh_2015.csv'  # Replace with your actual file path
# file_path = 'results/gsc/grd_2015.csv'  # Replace with your actual file path
data = pd.read_csv(file_path)

# Filter data for the specified region
lat_min, lat_max = 70.0, 82.0
lon_min, lon_max = -170.0, -120.0
region_data = data[(data['Latitude'] >= lat_min) & (data['Latitude'] <= lat_max) &
                   (data['Longitude'] >= lon_min) & (data['Longitude'] <= lon_max) ]
# region_data['Datetime'] = pd.to_datetime(region_data['Datetime'])
# region_data['Year'] = region_data['Datetime'].dt.year
# region_data_year = region_data[region_data['Year'] == 2018]

# Compute mean DH and geostrophic currents
# region_mean = data.groupby(['X_meters', 'Y_meters'])[['Surf_DH', 'ug', 'vg']].mean().reset_index()
# region_mean = data.groupby(['Longitude', 'Latitude'])[['Surf_DH', 'ug', 'vg']].mean().reset_index()
# region_mean = region_data.groupby(['Latitude', 'Longitude'])[['Surf_DH', 'ug', 'vg']].mean().reset_index()
region_mean = region_data.groupby(['X_meters', 'Y_meters'])[['Surf_DH', 'ug', 'vg']].mean().reset_index()

# Extract coordinates and data
# x = region_mean['Longitude'].values
# y = region_mean['Latitude'].values
x = region_mean['X_meters'].values
y = region_mean['Y_meters'].values
dhs = region_mean['Surf_DH'].values
ug = region_mean['ug'].values
vg = region_mean['vg'].values


# Plot dynamic height as a scatter plot
plt.figure(figsize=(10, 10))
# sc = plt.scatter(x, y, c=dhs, cmap="YlGnBu", s=90)
sc = plt.scatter(x, y, c=dhs, cmap="YlGnBu", s=280)
cb = plt.colorbar(sc, orientation="vertical", label="Dynamic Height (m)")

# Add geostrophic currents as vectors
scale = np.nanmax(np.sqrt(ug**2 + vg**2)) / 0.05
plt.quiver(x, y, ug, vg, scale=scale, color='black', width=0.002)

# Customize the plot
plt.title("Dynamic Height and Geostrophic Currents", fontsize=14)
plt.xlabel("X (meters)")
plt.ylabel("Y (meters)")
plt.axis("equal")  # Ensure equal scaling for x and y
plt.grid(True)

plt.show()

####### mean dh plot with mean gsc on map
# fig = plt.figure(figsize=(6, 4))    # change proj and plot
# projection = ccrs.NorthPolarStereo()
# ax = plt.axes(projection=projection)

# # Set the extent to focus on the region of interest
# # ax.set_extent([-180, 180, 70, 90], crs=ccrs.PlateCarree())
# ax.set_extent([-170, -120, 68, 82], crs=ccrs.PlateCarree())

# Add map features
# ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
# ax.add_feature(cfeature.BORDERS, linestyle=':', linewidth=0.5)
# ax.add_feature(cfeature.LAND, color='lightgrey')

# # Scatter plot for mean DH
# sc = plt.scatter(x, y, c=dhs, cmap='YlGnBu', s=180, transform=ccrs.PlateCarree(), label='Mean DH')
# cb = plt.colorbar(sc, orientation='vertical', pad=0.05, shrink=0.7)
# cb.set_label('Mean Dynamic Height (m)', fontsize=12)

# # Quiver plot for geostrophic currents (u_g, v_g)
# scale = np.nanmax(np.sqrt(ug**2 + vg**2)) / 0.1
# quiver = ax.quiver(x, y, ug, vg, transform=ccrs.PlateCarree(), scale=scale, color='black', width=0.003)

# # Add labels and title
# plt.title('Mean Dynamic Height and Geostrophic Currents (70°-82°N, 170°-120°W)', fontsize=14)
# plt.legend(loc='upper right')
# plt.show()
# plt.savefig('bg_map_2011_og.png', dpi=300)


###### annual dh and gsc maps
# region_data['Datetime'] = pd.to_datetime(region_data['Datetime'])
# region_data['Year'] = region_data['Datetime'].dt.year
# fig, axes = plt.subplots(1, 8, figsize=(10, 8), subplot_kw={'projection': ccrs.NorthPolarStereo()})
# fig.subplots_adjust(wspace=0.3)

# # Flatten the axes for easy iteration
# axes = axes.flatten()

# # Loop over each year (2011-2018)
# for i, year in enumerate(range(2011, 2019)):
#     # Filter data for the current year
#     region_data_year = region_data[region_data['Year'] == year]
    
#     # Extract the relevant variables
#     lons = region_data_year['Longitude'].values
#     lats = region_data_year['Latitude'].values
#     dhs = region_data_year['Surf_DH'].values  # Dynamic Height
#     # ug = region_data_year['ug'].values        # Geostrophic current u-component
#     # vg = region_data_year['vg'].values        # Geostrophic current v-component

#     # Get the current axis
#     ax = axes[i]
    
#     ax.set_extent([-180, -120, 68, 85], crs=ccrs.PlateCarree())

#     # Add map features
#     ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
#     ax.add_feature(cfeature.BORDERS, linestyle=':', linewidth=0.5)
#     ax.add_feature(cfeature.LAND, color='lightgrey')

#     # Plot the mean dynamic height (DH) as a scatter plot
#     sc = ax.scatter(lons, lats, c=dhs, cmap='YlGnBu', s=20, transform=ccrs.PlateCarree())

#     # Plot the geostrophic surface currents (GSC) as quiver vectors
#     # scale = np.nanmax(np.sqrt(ug**2 + vg**2)) / 0.05  # Adjust scaling for quiver vectors
#     # quiver = ax.quiver(lons, lats, ug, vg, scale=scale, width=0.003, color='black', transform=ccrs.PlateCarree())
#     # ax.quiverkey(quiver, 0.1, 0.1, 0.1, '0.1 m/s', labelpos='E', color='black')

#     # Set the title for each subplot
#     ax.set_title(f'{year}', fontsize=12)

#     # Remove axes box (for each subplot)
#     ax.set_frame_on(False)
    
# # Add a super title for the whole figure
# fig.suptitle('Mean Dynamic Height and Geostrophic Currents (70°-82°N, -170° to -120°E)', fontsize=16)

# # Plot the colorbar for the dynamic height (DH)
# cbar = plt.colorbar(sc, ax=axes, orientation='vertical', shrink=0.5)
# cbar.set_label('Mean Dynamic Height (m)', fontsize=10)

# # Show the plot
# plt.show()


######## time series of annual DH, gsc
# Convert 'Datetime' column to datetime type
# region_data['Datetime'] = pd.to_datetime(region_data['Datetime'])
# region_data['Year'] = region_data['Datetime'].dt.year

# # print(region_data.head)

# # Compute the annual mean DH
# annual_mean_dh = region_data.groupby('Year')['Surf_DH'].mean()
# annual_mean_ug = region_data.groupby('Year')['ug'].mean()
# annual_mean_vg = region_data.groupby('Year')['vg'].mean()

# fig, ax1 = plt.subplots(figsize=(10, 6))
# ax1.plot(annual_mean_dh.index, annual_mean_dh.values, marker='o', linestyle='-', color='red', label='DH')
# ax1.set_xlabel('Year', fontsize=12)
# ax1.set_ylabel('Mean Dynamic Height (m)', fontsize=12)
# ax1.tick_params(axis='y')
# ax1.grid(True, linestyle='--', alpha=0.6)
# ax1.set_ylim([0.5,0.8])
# ax2 = ax1.twinx()
# # Plot geostrophic velocities
# ax2.plot(annual_mean_ug.index, annual_mean_ug.values*100, marker='o', linestyle='-.', color='blue', label='Zonal velocity (ug)')
# ax2.plot(annual_mean_vg.index, annual_mean_vg.values*100, marker='o', linestyle='-.', color='green', label='Meridional velocity (vg)')
# ax2.set_ylabel('Mean Velocity (cm/s)', fontsize=12)
# ax2.tick_params(axis='y')

# # Add legends
# lines1, labels1 = ax1.get_legend_handles_labels()
# lines2, labels2 = ax2.get_legend_handles_labels()
# ax2.legend(lines1 + lines2, labels1 + labels2, loc='upper right', fontsize=10)

# plt.title('Annual Mean DH and Geostrophic Velocities (70°-82°N, -170° to -120°E)', fontsize=14)
# plt.tight_layout()
# # plt.show()
# plt.savefig('bg_dot_gsc.png', dpi=300)


###### section across BG - TPD - dot, gsc, surface/volume transport
# fig, ax1 = plt.subplots(figsize=(10, 6))
# ax1.plot(df['Year'], df['AO'], c='black', marker='.', linestyle='-', label='Annual AO')
# ax1.set_ylim([-0.7, 0.7])
# ax1.set_xlabel('Year', fontsize=12)
# ax1.set_ylabel('AO index', fontsize=12)
# ax1.tick_params(axis='y')
# ax1.grid(True, linestyle='--', alpha=0.6)

# for i, year in enumerate(df['Year']):
#     if year==2013: continue
#     if df['AO'][i] > 0:
#         ax1.axvspan(year - 0.1, year + 0.1, color='blue', alpha=0.2)  # Shaded region
#     if df['AO'][i] < 0:
#         ax1.axvspan(year - 0.1, year + 0.1, color='red', alpha=0.2)  # Shaded region

# ax2 = ax1.twinx()
# # Plot transport 

# ###### transport
# ## ALONG -170E
# lon_min = -171
# lon_max = -169
# lat_min, lat_max = 70.0, 82.0
# data_filtered = data[(data['Longitude'] >=lon_min) & (data['Longitude'] <=lon_max) &
#                      (data['Latitude'] >= lat_min) & 
#                      (data['Latitude'] <= lat_max)]

# data_filtered['Datetime'] = pd.to_datetime(data_filtered['Datetime'])
# data_filtered['Year'] = data_filtered['Datetime'].dt.year

# # Compute transport per segment (v_g * Δy)
# data_filtered['Transport'] = data_filtered['vg'] * 50000 # for volume * 400m

# # Group by year and compute total transport (in Sverdrups)
# # annual_transport170 = np.abs(data_filtered.groupby('Year')['Transport'].sum() / 1e6)  # Convert to Sv
# annual_transport170 = data_filtered.groupby('Year')['Transport'].sum() / 1e6  # Convert to Sv
# print(annual_transport170.values)

# print(np.corrcoef(annual_mean_dh.values, annual_transport170.values)[0][1])
# corr = np.corrcoef(df['AO'].values, annual_transport170.values)[0][1]
# ax2.plot(annual_mean_ug.index, annual_transport170.values, marker='o', linestyle='-.', color='green', label=f'-170E, corr={corr:.2f}')

# ### Along -160E
# lon_min = -161
# lon_max = -159
# lat_min, lat_max = 70.0, 82.0
# data_filtered = data[(data['Longitude'] >=lon_min) & (data['Longitude'] <=lon_max) &
#                      (data['Latitude'] >= lat_min) & 
#                      (data['Latitude'] <= lat_max)]

# data_filtered['Datetime'] = pd.to_datetime(data_filtered['Datetime'])
# data_filtered['Year'] = data_filtered['Datetime'].dt.year

# # Compute transport per segment (v_g * Δy)
# data_filtered['Transport'] = data_filtered['vg'] * 50000 # for volume * 400m

# # Group by year and compute total transport (in Sverdrups)
# # annual_transport170 = np.abs(data_filtered.groupby('Year')['Transport'].sum() / 1e6)  # Convert to Sv
# annual_transport170 = data_filtered.groupby('Year')['Transport'].sum() / 1e6  # Convert to Sv
# print(annual_transport170.values)

# # print(annual_transport.head)
# print(np.corrcoef(annual_mean_dh.values, annual_transport170.values)[0][1])
# corr = np.corrcoef(df['AO'].values, annual_transport170.values)[0][1]
# ax2.plot(annual_mean_ug.index, annual_transport170.values, marker='o', linestyle='-.', color='blue', label=f'-160E, corr={corr:.2f}')

# ##### Along -150E
# lon_min = -151
# lon_max = -149
# lat_min, lat_max = 70.0, 82.0
# data_filtered = data[(data['Longitude'] >=lon_min) & (data['Longitude'] <=lon_max) &
#                      (data['Latitude'] >= lat_min) & 
#                      (data['Latitude'] <= lat_max)]

# data_filtered['Datetime'] = pd.to_datetime(data_filtered['Datetime'])
# data_filtered['Year'] = data_filtered['Datetime'].dt.year

# # Compute transport per segment (v_g * Δy)
# data_filtered['Transport'] = data_filtered['vg'] * 50000 # for volume * 400m

# # Group by year and compute total transport (in Sverdrups)
# # annual_transport170 = np.abs(data_filtered.groupby('Year')['Transport'].sum() / 1e6)  # Convert to Sv
# annual_transport170 = data_filtered.groupby('Year')['Transport'].sum() / 1e6  # Convert to Sv
# print(annual_transport170.values)

# # print(annual_transport.head)
# print(np.corrcoef(annual_mean_dh.values, annual_transport170.values)[0][1])
# corr = np.corrcoef(df['AO'].values, annual_transport170.values)[0][1]
# ax2.plot(annual_mean_ug.index, annual_transport170.values, marker='o', linestyle='-.', color='cyan', label=f'-150E, corr={corr:.2f}')

# ##### Along -140E
# lon_min = -141
# lon_max = -139
# lat_min, lat_max = 70.0, 82.0
# data_filtered = data[(data['Longitude'] >=lon_min) & (data['Longitude'] <=lon_max) &
#                      (data['Latitude'] >= lat_min) & 
#                      (data['Latitude'] <= lat_max)]

# data_filtered['Datetime'] = pd.to_datetime(data_filtered['Datetime'])
# data_filtered['Year'] = data_filtered['Datetime'].dt.year

# # Compute transport per segment (v_g * Δy)
# data_filtered['Transport'] = data_filtered['vg'] * 50000 # for volume * 400m

# # Group by year and compute total transport (in Sverdrups)
# # annual_transport170 = np.abs(data_filtered.groupby('Year')['Transport'].sum() / 1e6)  # Convert to Sv
# annual_transport170 = data_filtered.groupby('Year')['Transport'].sum() / 1e6  # Convert to Sv
# print(annual_transport170.values)

# # print(annual_transport.head)
# print(np.corrcoef(annual_mean_dh.values, annual_transport170.values)[0][1])
# corr = np.corrcoef(df['AO'].values, annual_transport170.values)[0][1]
# ax2.plot(annual_mean_ug.index, annual_transport170.values, marker='o', linestyle='-.', color='magenta', label=f'-140E, corr={corr:.2f}')


# ##### along -130E
# lon_min = -131
# lon_max = -129
# lat_min, lat_max = 70.0, 82.0
# data_filtered = data[(data['Longitude'] >=lon_min) & (data['Longitude'] <=lon_max) &
#                      (data['Latitude'] >= lat_min) & 
#                      (data['Latitude'] <= lat_max)]

# data_filtered['Datetime'] = pd.to_datetime(data_filtered['Datetime'])
# data_filtered['Year'] = data_filtered['Datetime'].dt.year

# # Compute transport per segment (v_g * Δy)
# data_filtered['Transport'] = data_filtered['vg'] * 50000 # for volume * 400m

# # Group by year and compute total transport (in Sverdrups)
# # annual_transport170 = np.abs(data_filtered.groupby('Year')['Transport'].sum() / 1e6)  # Convert to Sv
# annual_transport170 = data_filtered.groupby('Year')['Transport'].sum() / 1e6  # Convert to Sv
# print(annual_transport170.values)

# # print(annual_transport.head)
# print(np.corrcoef(annual_mean_dh.values, annual_transport170.values)[0][1])
# corr = np.corrcoef(df['AO'].values, annual_transport170.values)[0][1]
# ax2.plot(annual_mean_ug.index, annual_transport170.values, marker='o', linestyle='-.', color='red', label=f'-130E, corr={corr:.2f}')


# # ax2.set_ylabel('Transport (Sv)', fontsize=12)
# ax2.tick_params(axis='y')
# # ax2.set_ylim([-1, -0.6])

# # Add legends
# lines1, labels1 = ax1.get_legend_handles_labels()
# lines2, labels2 = ax2.get_legend_handles_labels()
# ax2.legend(lines1 + lines2, labels1 + labels2, loc='lower right', fontsize=8)

# ax2.set_ylabel(r'$ Transport[Sv]\ \longleftarrow southward\ northward\ \longrightarrow$', fontsize=12)

# plt.title('AO index and Transport along BG section', fontsize=14)
# plt.tight_layout()
# # plt.show()
# plt.savefig('bg_ao_transport_ts.png', dpi=300)




