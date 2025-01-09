import netCDF4 as nc
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
from helpers import *
import csv
from pyproj import Proj, transform, Transformer
import cartopy.crs as ccrs
import cartopy.feature as cfeature

####### monthly SLA data
from saga import time, lats, lons, dot, ug, vg
# time = time[:96]  # Slice the 1D array time (no additional dimensions required)

# Compute the mean over the first 96 time steps for 3D arrays
# dot = np.mean(dot[:96, :, :], axis=0)  # Average over time dimension (axis 0)
# ug = np.mean(ug[:96, :, :], axis=0)
# vg = np.mean(vg[:96, :, :], axis=0)
# dot = dot[0,:,:]
# ug = ug[0,:,:]
# vg = vg[0,:,:]

# print(time)  # Print the new time array after removing the last 24 steps

# Mask the data to only include latitudes between 70 and 82, and longitudes between -170 and -120
lat_min, lat_max = 70.0, 82.0
lon_min, lon_max = -170.0, -120.0

# Create masks for latitude and longitude
lat_mask = (lats >= lat_min) & (lats <= lat_max)
lon_mask = (lons >= lon_min) & (lons <= lon_max)

lats = lats[lat_mask]
lons = lons[lon_mask]
dot = dot[np.ix_(lat_mask, lon_mask)]
ug = ug[np.ix_(lat_mask, lon_mask)]
vg = vg[np.ix_(lat_mask, lon_mask)] # Slice dot for time and latitudes, keeping all longitudes

# Check the shapes of the resulting arrays
print(f"dot shape: {dot.shape}")
print(f"ug shape: {ug.shape}")
print(f"vg shape: {vg.shape}")

lons, lats = np.meshgrid(lons, lats)
lons = np.array(lons).flatten()
lats = np.array(lats).flatten()
dhs = np.array(dot).flatten()
ug = np.array(ug).flatten()
vg = np.array(vg).flatten()

####### mean dh plot with mean gsc on map
fig = plt.figure(figsize=(6, 4))
projection = ccrs.NorthPolarStereo()
ax = plt.axes(projection=projection)

# Set the extent to focus on the region of interest
ax.set_extent([-180, -120, 68, 85], crs=ccrs.PlateCarree())

# Add map features
ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
ax.add_feature(cfeature.BORDERS, linestyle=':', linewidth=0.5)
ax.add_feature(cfeature.LAND, color='lightgrey')

# Scatter plot for mean DH
sc = ax.scatter(lons, lats, c=dhs, cmap='YlGnBu', s=180, transform=ccrs.PlateCarree(), label='Mean DH')
cb = plt.colorbar(sc, orientation='vertical', pad=0.05, shrink=0.7)
cb.set_label('Mean Dynamic Height (m)', fontsize=12)

# Quiver plot for geostrophic currents (u_g, v_g)
scale = np.nanmax(np.sqrt(ug**2 + vg**2)) / 0.03
quiver = ax.quiver(lons, lats, ug, vg, transform=ccrs.PlateCarree(), scale=scale, color='black', width=0.003)

# Add labels and title
plt.title('Mean Dynamic Height and Geostrophic Currents (70°-82°N, -170° to -120°E)', fontsize=14)
plt.savefig('bg_saga_2015.png', dpi=300)
# plt.show()

####### Plot of dot, ug, vg time series
# from saga import time, lats, lons, dot, ug, vg

# # Define region of interest
# lat_min, lat_max = 70.0, 82.0
# lon_min, lon_max = -170.0, -120.0

# # Monthly time series: slice only the first 96 time steps (8 years)
# time = time[:96]
# dot = dot[:96, :, :]
# ug = ug[:96, :, :]
# vg = vg[:96, :, :]

# # Create masks for latitude and longitude
# lat_mask = (lats >= lat_min) & (lats <= lat_max)
# lon_mask = (lons >= lon_min) & (lons <= lon_max)

# # Compute annual means for dot, ug, and vg
# dot_annual_mean = []
# ug_annual_mean = []
# vg_annual_mean = []

# years = np.arange(2011, 2019)  # Years from 2011 to 2018
# mon_idx = 0  # Initialize monthly index

# for year in years:
#     # Slice data for the current year (12 months)
#     dot_year = dot[mon_idx:mon_idx+12, :, :]
#     ug_year = ug[mon_idx:mon_idx+12, :, :]
#     vg_year = vg[mon_idx:mon_idx+12, :, :]
    
#     dot_region_mean = np.nanmean(dot_year[:, lat_mask, :][:, :, lon_mask])
#     ug_region_mean = np.nanmean(ug_year[:, lat_mask, :][:, :, lon_mask])
#     vg_region_mean = np.nanmean(vg_year[:, lat_mask, :][:, :, lon_mask])
    
#     # Append the annual means
#     dot_annual_mean.append(dot_region_mean)
#     ug_annual_mean.append(ug_region_mean)
#     vg_annual_mean.append(vg_region_mean)
    
#     # Increment the monthly index by 12
#     mon_idx += 12

# # Convert lists to numpy arrays
# dot_annual_mean = np.array(dot_annual_mean)
# ug_annual_mean = np.array(ug_annual_mean)
# vg_annual_mean = np.array(vg_annual_mean)

# # Print the annual means
# print("Annual DOT Means:", dot_annual_mean)
# print("Annual UG Means:", ug_annual_mean)
# print("Annual VG Means:", vg_annual_mean)

# # Plot the annual time series of DOT
# fig, ax1 = plt.subplots(figsize=(10, 6))
# ax1.plot(years, dot_annual_mean, marker='o', linestyle='-', color='red', label='DOT')
# ax1.set_xlabel('Year', fontsize=12)
# ax1.set_ylabel('Mean Dynamic Height (m)', fontsize=12)
# ax1.tick_params(axis='y')
# ax1.grid(True, linestyle='--', alpha=0.6)

# ax2 = ax1.twinx()
# # Plot geostrophic velocities
# ax2.plot(years, ug_annual_mean*100, marker='o', linestyle='-', color='blue', label='Zonal velocity (ug)')
# ax2.plot(years, vg_annual_mean*100, marker='o', linestyle='-', color='green', label='Meridional velocity (vg)')
# ax2.set_ylabel('Mean Velocity (cm/s)', fontsize=12)
# ax2.tick_params(axis='y')

# # Add legends
# lines1, labels1 = ax1.get_legend_handles_labels()
# lines2, labels2 = ax2.get_legend_handles_labels()
# ax2.legend(lines1 + lines2, labels1 + labels2, loc='upper right', fontsize=10)
# # ax1.set_ylim([0.5, 0.7])
# # ax2.set_ylim([-0.4, 0.4])

# plt.title('Annual Mean DH and Geostrophic Velocities (70°-82°N, -170° to -120°E)', fontsize=14)
# plt.tight_layout()
# # plt.savefig('bg_dot_gsc_saga.png', dpi=300)
# plt.show()

####### ao
# from ao import *
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

# ##### transport 
# ## ALONG -170E
# lon_min = -171
# lon_max = -169
# lat_min, lat_max = 70.0, 82.0
# lat_mask = (lats >= lat_min) & (lats <= lat_max)
# lon_mask = (lons >= lon_min) & (lons <= lon_max)
# dot_region_mean = np.nanmean(dot_year[:, lat_mask, :][:, :, lon_mask], axis=[1,2])
# ug_region_mean = np.nanmean(ug_year[:, lat_mask, :][:, :, lon_mask], axis=[1,2])
# vg_region_mean = np.nanmean(vg_year[:, lat_mask, :][:, :, lon_mask], axis=[1,2])

# Compute transport per segment (v_g * Δy)
# data_filtered['Transport'] = ug['vg'] * 50000 # for volume * 400m

# # Group by year and compute total transport (in Sverdrups)
# # annual_transport170 = np.abs(data_filtered.groupby('Year')['Transport'].sum() / 1e6)  # Convert to Sv
# annual_transport170 = data_filtered.groupby('Year')['Transport'].sum() / 1e6  # Convert to Sv
# print(annual_transport170.values)

# print(np.corrcoef(annual_mean_dh.values, annual_transport170.values)[0][1])
# corr = np.corrcoef(df['AO'].values, annual_transport170.values)[0][1]
# ax2.plot(annual_mean_ug.index, annual_transport170.values, marker='o', linestyle='-.', color='green', label=f'-170E, corr={corr:.2f}')

### Along -160E
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
# plt.show()