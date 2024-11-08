import pandas as pd
import gsw
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from datetime import datetime
from scipy import interpolate
import csv

# print('Reading all hydrographic data...')
# hydr_data = {}

def read_data(fname, hydr_data):
    data = pd.read_csv(fname)

    dt = data['Datetime'].values
    temperature = data['Temperature'].values
    salinity = data['Salinity'].values
    pressure = data['Pressure'].values
    latitude = data['Latitude'].values  
    longitude = data['Longitude'].values  
    depth = data['Depth'].values  

    print('Set the data with key as lat lon and val as T/S/P/D/SSH...')
    for idx in range(len(latitude)):

        ts = datetime.strptime(dt[idx], "%Y-%m-%d %H:%M:%S")
        # take only summer months measurements to remove bias
        # if ts.month == 6 or ts.month == 7 or ts.month == 8:
            
        if (latitude[idx],longitude[idx]) not in list(hydr_data.keys()):
            hydr_data[(latitude[idx], longitude[idx])] = {'dt': ts.date(), 'temp':[], 'sali':[], 'pres':[], \
                                                          'dept':[], 'ssh':np.nan}
            
        hydr_data[(latitude[idx], longitude[idx])]['temp'].append(temperature[idx])
        hydr_data[(latitude[idx], longitude[idx])]['sali'].append(salinity[idx])
        hydr_data[(latitude[idx], longitude[idx])]['pres'].append(pressure[idx])
        hydr_data[(latitude[idx], longitude[idx])]['dept'].append(depth[idx])

    # print(hydr_data)
    return hydr_data

# hydr_data = read_data('data/2011.csv', hydr_data)
# hydr_data = read_data('data/2012.csv', hydr_data)
# hydr_data = read_data('data/2013.csv', hydr_data)
# hydr_data = read_data('data/2014.csv', hydr_data)
# hydr_data = read_data('data/2015.csv', hydr_data)
# hydr_data = read_data('data/2016.csv', hydr_data)
# hydr_data = read_data('data/2017.csv', hydr_data)
# hydr_data = read_data('data/2018.csv', hydr_data)
# print('Done reading all hydrographic data...')

# def get_g(lat):
#     # Define constants
#     g0 = 9.780327  # Standard gravity at the equator in m/s^2
#     phi = np.radians(lat) 
#     g = g0 * (1 + 0.0053024 * np.sin(phi)**2 - 0.0000058 * np.sin(2 * phi)**2)
#     return g

# fp = open('dh_2011_2018.csv', 'w+', newline='')
# writer = csv.writer(fp, delimiter=',')
# writer.writerow(['Datetime','Latitude','Longitude','Dynamic_height'])


# ############## using specific volume anomaly method 
# # for each latlon compute the dynamic height using t, s, p and absolute ssh 
# i=0
# for latlon,data in hydr_data.items():
    
#     if max(data['dept']) > 500:
#         # interpolate upto 400m depth for each 2m step
#         D = [d for d in data['dept'] if d<=400]
#         T = [data['temp'][idx] for idx in range(len(D))]
#         S = [data['sali'][idx] for idx in range(len(D))]
#         P = [data['pres'][idx] for idx in range(len(D))]
        
#         if len(D) > 0:
#             xdept = np.linspace(0,400,201)
#             ytemp = interpolate.interp1d(D, T, fill_value='extrapolate')(xdept)
#             ysali = interpolate.interp1d(D, S, fill_value='extrapolate')(xdept)
#             ypres = interpolate.interp1d(D, P, fill_value='extrapolate')(xdept)
            
#             SA = gsw.SA_from_SP(ysali, ypres, latlon[1], latlon[0])
#             CT = gsw.CT_from_t(SA, ytemp, ypres)
        
#             try:  
#                 # Integrate specific volume anomaly to compute dynamic height https://www.teos-10.org/pubs/gsw/html/gsw_geo_strf_dyn_height.html
#                 dyn_height_anom = gsw.geostrophy.geo_strf_dyn_height(SA, CT, ypres, p_ref=400)
#                 ster_height = dyn_height_anom[0] / get_g(latlon[0]) 
                
#                 # DOT or absolute ssh is the dynamic height at the surface pressure 0dbar
#                 hydr_data[latlon]['ssh'] = ster_height
#                 print('Setting steric height...', hydr_data[latlon]['ssh'])
                
#             except ValueError as err:
#                 print('Setting DH to nan due to error...')
#                 hydr_data[latlon]['ssh'] = None
#                 i=i+1
                
#         writer.writerow([hydr_data[latlon]['dt'],latlon[0],latlon[1],hydr_data[latlon]['ssh']])
        
# print(i) 
# fp.close()

def plot_dh(hydr_data, year):
    ####### plot dh of hydrographic data
    fig = plt.figure(figsize=(8, 8))
    ax = plt.axes(projection=ccrs.NorthPolarStereo())

    # Add land and coastlines
    ax.add_feature(cfeature.LAND, zorder=0, edgecolor='black')
    ax.add_feature(cfeature.COASTLINE, zorder=0, edgecolor='black')

    # lats = np.array([latlon[0] for latlon in hydr_data.keys()])
    # lons = np.array([latlon[1] for latlon in hydr_data.keys()])
    # ssh = np.array([hydr_data[(lats[idx], lons[idx])]['ssh'] for idx in range(len(lats))])
    lats = hydr_data['Latitude']
    lons = hydr_data['Longitude']
    ssh = hydr_data['Dynamic_height']
    print('Max, Min:', max(ssh), min(ssh))

    # Plot the scatter 
    sc = ax.scatter(lons, lats, c=ssh, cmap='YlGnBu', s=8, transform=ccrs.PlateCarree(), vmin=0, vmax=1.2)

    # Add a colorbar
    cbar = plt.colorbar(sc, orientation='vertical', pad=0.1)
    cbar.set_label('Dynamic Height (m)')

    ax.set_extent([-180, 180, 60, 90], ccrs.PlateCarree())
    ax.gridlines(draw_labels=True)
    # plt.title('Dynamic height in', year)
    # plt.show()

    plt.savefig(f'dh_{year}_deep.png', dpi=300)    
    
# plot_dh(hydr_data, '2011')      
# # plot_dh(hydr_data, '2012')      
# # plot_dh(hydr_data, '2013')      
# # plot_dh(hydr_data, '2014')      
# # plot_dh(hydr_data, '2015')      
# # plot_dh(hydr_data, '2016')      
# # plot_dh(hydr_data, '2017')      
# # plot_dh(hydr_data, '2018')    

hydr_data = pd.read_csv('results/dh_2011_2018.csv') 
plot_dh(hydr_data, '2011_2018')      
