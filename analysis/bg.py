import cartopy.crs as ccrs
from mpl_toolkits.basemap import Basemap
import cartopy.feature as cfeature
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
# from bathymetric_depth import *
from ao import *
from pyproj import Proj, Transformer

file_path = 'grd_gsc_1980_2018.csv'  # Replace with your actual file path
data = pd.read_csv(file_path)

# Filter data for the specified region
# lat_min, lat_max = 70.0, 81.0
# lon_min, lon_max = -180.0, -130.0

# FS
# lat_min, lat_max = 86, 87.0
# lon_min, lon_max = -45.0, 10.0

# LR start
lat_min, lat_max = 81.0, 82.0
lon_min, lon_max = 120.0, 160.0
region_data = data[(data['Latitude'] >= lat_min) & (data['Latitude'] <= lat_max) &
                   (data['Longitude'] >= lon_min) & (data['Longitude'] <= lon_max) ]
# region_mean = region_data.groupby(['X_meters', 'Y_meters'])[['Surf_DH', 'ug', 'vg']].mean().reset_index()
# region_mean = data.groupby(['X_meters', 'Y_meters'])[['Surf_DH', 'ug', 'vg']].mean().reset_index()

# Extract coordinates and data
# x = region_mean['X_meters'].values
# y = region_mean['Y_meters'].values
# dhs = region_mean['Surf_DH'].values
# ug = region_mean['ug'].values
# vg = region_mean['vg'].values

###### x, y, plot
# plt.figure(figsize=(10, 10))
# # sc = plt.scatter(x, y, c=dhs, cmap="YlGnBu", s=90)
# sc = plt.scatter(x, y, c=dhs, cmap="YlGnBu", s=350)
# cb = plt.colorbar(sc, orientation="vertical", label="Dynamic Height (m)")

# # Add geostrophic currents as vectors
# scale = np.nanmax(np.sqrt(ug**2 + vg**2)) / 0.1
# plt.quiver(x, y, ug, vg, scale=scale, color='black', width=0.003)   # add quiver legend

# # plt.title('Beaufort gyre 1980-2018 Mean DH and Gsc')
# plt.axis("equal")  # Ensure equal scaling for x and y
# plt.grid(True)
# plt.show()
# plt.savefig('bg_1980_2018_mean_xy.png', dpi=300)

# lons = region_mean['Longitude'].values
# lats = region_mean['Latitude'].values

# # Create a polar stereographic projection
# proj = ccrs.NorthPolarStereo()

# # Initialize the figure and axes
# plt.figure(figsize=(12, 12))
# ax = plt.axes(projection=proj)
# ax.set_extent([-180, 180, 68, 90], crs=ccrs.PlateCarree())  # Define the map extent

# # Add map features
# ax.add_feature(cfeature.LAND, facecolor='lightgray')
# ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
# ax.gridlines(draw_labels=True)

# # Transform coordinates to the projection
# transform = ccrs.PlateCarree()
# x, y = lons, lats  # Assuming lon and lat are arrays of longitudes and latitudes

# # Scatter plot for dynamic height
# sc = plt.scatter(x, y, c=dhs, cmap="YlGnBu", s=350, transform=transform)
# cb = plt.colorbar(sc, orientation="vertical", label="Dynamic Height (m)")

# # Add geostrophic currents as vectors
# scale = np.nanmax(np.sqrt(ug**2 + vg**2)) / 0.1
# plt.quiver(x, y, ug, vg, scale=scale, color='black', width=0.003, transform=transform)

# # Set title
# plt.title('Dynamic Height and Geostrophic Currents on North Polar Stereographic Projection')
# plt.show()


####### annnual maps of dot # TODO compute trends 2011-2018
# region_data['Datetime'] = pd.to_datetime(region_data['Datetime'])
# region_data['Year'] = region_data['Datetime'].dt.year
# fig, axes = plt.subplots(2, 4, figsize=(8, 8))

# # Flatten the axes for easy iteration
# axes = axes.flatten()

# # Loop over each year (2011-2018)
# for i, year in enumerate(range(2011, 2019)):
#     # Filter data for the current year
#     region_data_year = region_data[region_data['Year'] == year]
    
#     # Extract the relevant variables
#     x = region_data_year['X_meters'].values
#     y = region_data_year['Y_meters'].values
#     dhs = region_data_year['Surf_DH'].values  # Dynamic Height
#     ug = region_data_year['ug'].values
#     vg = region_data_year['vg'].values

#     ax = axes[i]
    
#     sc = ax.scatter(x, y, c=dhs, cmap="YlGnBu", s=90)
#     # cb = ax.colorbar(sc, orientation="vertical", label="Dynamic Height (m)")

#     # Add geostrophic currents as vectors
#     scale = np.nanmax(np.sqrt(ug**2 + vg**2)) / 0.1
#     quiver = ax.quiver(x, y, ug, vg, scale=scale, color='black', width=0.002)   # add quiver legend
#     ax.quiverkey(quiver, 0.1, 0.1, 0.1, '0.1 m/s', labelpos='E', color='black')
    
#     ax.set_title(f'{year}', fontsize=12)
#     # plt.axis("equal")  # Ensure equal scaling for x and y
#     # plt.grid(True)

# cbar = plt.colorbar(sc, ax=axes, orientation='vertical', shrink=0.5)
# cbar.set_label('Mean Dynamic Height (m)', fontsize=10)
# plt.show()


######## time series of annual DH, gsc SAGA
region_data['Datetime'] = pd.to_datetime(region_data['Datetime'])
region_data['Year'] = region_data['Datetime'].dt.year
region_data['Month'] = region_data['Datetime'].dt.month
# print(region_data['Month'])

# region_data = region_data[(region_data['Month']==7) | (region_data['Month']==8) | (region_data['Month']==9) | (region_data['Month']==10)]
# region_data = region_data[(region_data['Month']==11) | (region_data['Month']==12) | (region_data['Month']==1) | 
#                           (region_data['Month']==2) | (region_data['Month']==3) | (region_data['Month']==4) |
#                           (region_data['Month']==5) | (region_data['Month']==6)]
# Compute the annual mean DH
annual_mean_dh = region_data.groupby('Year')['Surf_DH'].mean()
annual_mean_ug = region_data.groupby('Year')['ug'].mean()
annual_mean_vg = region_data.groupby('Year')['vg'].mean()

print(max(annual_mean_dh))  # SUMMER=0.659368508697068, ANNUAL=0.6655937294554898, WINTER=0.674810271346563

years = range(2011, 2019)
# bg from SAGA
# annual
# dot = [0.57778235,0.56724194,0.57190907,0.57744078,0.58295251,0.59840932,0.58362414,0.59542531,0.59836627,0.58921805][0:8]
# print(max(dot))
# # summer
# dot = [0.5947648348405712, 0.5732655636484042, 0.5825672944503133, 0.5954438863750228, 0.6063547770995144, 0.6060834147943638, 0.6021691799720162, 0.6054839308218327, 0.6229220802438532, 0.5978694310515648][0:8]
# print(max(dot))
# # # winter 
# dot = [0.5692911108039135, 0.564230133520331, 0.5665799565382659, 0.5684392228582397, 0.5712513821758182, 0.5945722665921782, 0.574351624399514, 0.590395997402988, 0.586088360960378, 0.5848923545688148][0:8]
# print(max(dot)) # annual=0.59840932, summer=0.6063547770995144, winter=0.5945722665921782

# dot = 
# ug = [i**2 for i in [-0.00166289, 0.00097516, -0.00171596, -0.00297301, -0.00284071, -0.00065336, 0.00092955, -0.00027355]]
# vg = [i**2 for i in [0.0069653, 0.00612324, 0.00506053, 0.00575097, 0.00511311, 0.00567128, 0.00836546, 0.00817728]]

# FS from SAGA
# annual 
# dot = [0.23133023,0.19609327,0.22711721,0.21471036,0.20398841,0.24041543,0.22145704,0.22920151,0.23966722,0.23550441][0:8]
# print(max(dot))

# summer
# dot = [0.2530981630237924, 0.21462254002690315, 0.2508420189943265, 0.2409440826171556, 0.2477739525304453, 0.26545400730762125, 0.23981466028033882, 0.24192172050778124, 0.26220591229864876, 0.25340970569365734][0:8]
# print(max(dot))

# winter
# dot = [0.22044626912605522, 0.1868286419659853, 0.2152548055617592, 0.20159349768372203, 0.18209564412339918, 0.2278961441064304, 0.21227823237108218, 0.2228414084938531, 0.22839786705082735, 0.22655176163424512][0:8]
# print(max(dot))

# LR start from SAGA
# annual 
# dot = [0.2266353,0.20565696,0.21339384,0.21704378,0.19703093,0.20158889,0.17114796,0.19286856,0.18545996,0.18372723][0:8]
# print(max(dot))
# summer 
dot = [0.25390159081391717, 0.21584939222644878, 0.22472804998634038, 0.2452059229042519, 0.23210521613558133, 0.2096662828008886, 0.1766658088537278, 0.19899534476507041, 0.2094545527395827, 0.20474409173691163][0:8]
print(max(dot))

# # winter
# dot = [0.21300215742516296, 0.20056075085023486, 0.20772673088746765, 0.20296270107953912, 0.1794937835989037, 0.19755019523824255, 0.168389034545463, 0.1898051726166159, 0.173462661286747, 0.1732187978474906][0:8]
# print(max(dot))

# ug = [i**2 for i in [ 0.00502478,-0.001866,-0.00196271,-0.00574145,-0.00136061,-0.00077899,-0.00227324,-0.00179609]]
# vg = [i**2 for i in [0.023437,0.02278752,0.02041138,0.0266506,0.02799238,0.02613746,0.02126416,0.02323556]]
# currents = np.add(np.array(ug), np.array(vg))

# from transport import annual_transport_df
from scipy.stats import pearsonr

fig, ax1 = plt.subplots(figsize=(10, 6))
ax1.plot(annual_mean_dh.index, annual_mean_dh.values, marker='o', linestyle='-', color='red', label='Mean DH rel.to 400m (UDAH)')
ax1.plot(years, dot, marker='o', linestyle='-', color='blue', label='DOT (SAGA/Cryosat2)')

correlation, p_value = pearsonr(dot, annual_mean_dh.values[-8:])    # annual=0.776, summer=0.230, winter=0.830
print(f"Pearson Correlation Coefficient: {correlation:.3f}")
print(f"P-value: {p_value:.3e}")

ax1.set_xlabel('Year', fontsize=12)
ax1.set_ylabel('Mean Dynamic Height (m)', fontsize=12)
ax1.tick_params(axis='y')

# ax2 = ax1.twinx()
# ax2.plot(range(1988, 2019), annual_transport_df.values*100, linestyle='-', color='black', label='weighted mean velocity')
# ax2.plot(annual_mean_dh.index, np.sqrt(annual_mean_ug.values**2+annual_mean_vg.values**2)*100, marker='o', linestyle='-', color='black', label='Mean GSC (UDAH)')
# ax2.plot(years, np.sqrt(np.array(currents)), marker='o', linestyle='-', color='cyan', label='Mean GSC (SAGA/Cryosat2)')
# ax2.set_ylabel('Mean Velocity (cm/s)', fontsize=12)
# ax2.tick_params(axis='y')

# lines1, labels1 = ax1.get_legend_handles_labels()
# lines2, labels2 = ax2.get_legend_handles_labels()
# ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper left', fontsize=10)

# # ax1.set_ylim([0.4, 0.7])
# # ax2.set_ylim([-0.4, 0.4])

ax1.grid(True, linestyle='--', alpha=0.6)
# # plt.tight_layout()
plt.legend()
# plt.title('Transport drift (LR start) mean DH and GSC 1980-2018')
plt.title('Transport drift (end of LR-Fram side) mean DH 1980-2018')
# plt.title('BG SUMMER mean DH and GSC 1980-2018')
# plt.savefig('tpd_1980_2018_dot_gsc.png', dpi=300)
plt.show()
