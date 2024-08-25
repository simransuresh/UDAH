import pandas as pd
import gsw
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from datetime import datetime
from scipy import interpolate

# Read the data from the text file
data = pd.read_csv('2011.csv')

print('Reading all hydrographic data...')

# Extract columns
dt = data['Datetime'].values
temperature = data['Temperature'].values
salinity = data['Salinity'].values
pressure = data['Pressure'].values
latitude = data['Latitude'].values  
longitude = data['Longitude'].values  
depth = data['Depth'].values  

print('Done reading all hydrographic data...')

# For each lat lon, get the temp, sali and pressure profile for all depth and save as key
hydr_data = {}

for idx in range(len(latitude)):

    ts = datetime.strptime(dt[idx], "%Y-%m-%d %H:%M:%S")

    # take only summer months measurements to remove bias
    # if ts.month == 6 or ts.month == 7 or ts.month == 8:
        
    if (latitude[idx],longitude[idx]) not in list(hydr_data.keys()):
        hydr_data[(latitude[idx], longitude[idx])] = {'temp':[], 'sali':[], 'pres':[], 'dept':[], 'ssh':np.nan}
        
    hydr_data[(latitude[idx], longitude[idx])]['temp'].append(temperature[idx])
    hydr_data[(latitude[idx], longitude[idx])]['sali'].append(salinity[idx])
    hydr_data[(latitude[idx], longitude[idx])]['pres'].append(pressure[idx])
    hydr_data[(latitude[idx], longitude[idx])]['dept'].append(depth[idx])
        
print('Set the data with key as lat lon and val as T/S/P/D/SSH...')


############## using specific volume anomaly method 
# for each latlon compute the dynamic height using t, s, p and absolute ssh 
i=0
for latlon,data in hydr_data.items():
    
    # interpolate upto 400m depth for each 2m step
    D = [d for d in data['dept'] if d<=400]
    T = [data['temp'][idx] for idx in range(len(D))]
    S = [data['sali'][idx] for idx in range(len(D))]
    P = [data['pres'][idx] for idx in range(len(D))]
    print('Retrieving depth upto 400m ... ')
    
    if len(D) > 0:
        xdept = np.linspace(0,400,201)
        ytemp = interpolate.interp1d(D, T, fill_value='extrapolate')(xdept)
        ysali = interpolate.interp1d(D, S, fill_value='extrapolate')(xdept)
        ypres = interpolate.interp1d(D, P, fill_value='extrapolate')(xdept)
        print('Interpolating for every 2m upto 400m...') 
        
        # Compute Absolute Salinity from practical salinity 
        SA = gsw.SA_from_SP(ysali, ypres, latlon[1], latlon[0])
        print('Computing absolute salinity...')
        
        # Compute Conservative Temperature from SA and temperature
        CT = gsw.CT_from_t(SA, ytemp, ypres)
        print('Computing conservative temperature...')
        
        try:  
            # Integrate specific volume anomaly to compute dynamic height 
            # REF: https://www.teos-10.org/pubs/gsw/html/gsw_geo_strf_dyn_height.html
            dyn_height_anom = gsw.geostrophy.geo_strf_dyn_height(SA, CT, ypres, p_ref=400)
            print('Computing dynamic height anomaly...')
            
            # REF: https://www.teos-10.org/pubs/gsw/html/gsw_geo_strf_steric_height.html
            ster_height = dyn_height_anom[0] / 9.8 # TODO change g wrt lat
            
            # for pref = 400 and positive SH and set ssh as SH[0] => min, max -> 1.9120461764456176 -4.042504447640991
            # for pref = 0 and negative and and set ssh as SH[-1] => same
            
            # DOT or absolute ssh is the dynamic height at the surface pressure 0dbar
            hydr_data[latlon]['ssh'] = ster_height
            print('Setting steric height...', hydr_data[latlon]['ssh'])
            
        except ValueError as err:
            print('Setting DH to nan due to error...')
            hydr_data[latlon]['ssh'] = np.nan
            i=i+1
        
print(i)    # 5 profiles throw error of pressure not monotonically increasing

            
####### plot dh of hydrographic data
# fig = plt.figure(figsize=(8, 8))
# ax = plt.axes(projection=ccrs.NorthPolarStereo())

# # Add land and coastlines
# ax.add_feature(cfeature.LAND, zorder=0, edgecolor='black')
# ax.add_feature(cfeature.COASTLINE, zorder=0, edgecolor='black')

# lats = np.array([latlon[0] for latlon in hydr_data.keys()])
# lons = np.array([latlon[1] for latlon in hydr_data.keys()])
# ssh = np.array([hydr_data[(lats[idx], lons[idx])]['ssh'] for idx in range(len(lats))])
# print('Max, Min:', max(ssh), min(ssh))

# # Plot the scatter 
# sc = ax.scatter(lons, lats, c=ssh, cmap='YlGnBu', s=10, transform=ccrs.PlateCarree(), vmin=0, vmax=1.2)

# # Add a colorbar
# cbar = plt.colorbar(sc, orientation='vertical', pad=0.1)
# cbar.set_label('Dynamic Height (m)')

# ax.set_extent([-180, 180, 60, 90], ccrs.PlateCarree())
# ax.gridlines(draw_labels=True)
# plt.title('Dynamic height in 2011')
# plt.show()

# plt.savefig('dh_2011_400m.png', dpi=300)


# ###### Histogram: plot ssh of all stations 
# # data to histogram
# dist1 = [val['ssh'] for key,val in hydr_data.items()]

# # Plot the color histogram
# plt.figure(figsize=(10, 6))
# n, bins, patches = plt.hist(dist1, bins=30, edgecolor='black', alpha=0.7)

# # Normalize the color scale (colormap)
# bin_centers = 0.5 * (bins[:-1] + bins[1:])
# colormap = plt.cm.viridis((bin_centers - min(bin_centers)) / (max(bin_centers) - min(bin_centers)))

# # Apply the colormap to each patch
# for patch, color in zip(patches, colormap):
#     patch.set_facecolor(color)

# # Add labels and title
# plt.xlabel('Sea Surface Height (m)')
# plt.ylabel('Frequency')
# plt.title('Color Histogram of Sea Surface Height')

# # Show the plot
# plt.show()




