import pandas as pd
import gsw
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from datetime import datetime
from scipy import interpolate
import csv
from helpers import get_g
from bathymetric_depth import nearest_depth

print('Reading all hydrographic data...')
hydr_data = {}

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
        
        latitude[idx], longitude[idx] = round(latitude[idx],5), round(longitude[idx],5)
        ts = datetime.strptime(dt[idx], "%Y-%m-%d %H:%M:%S")
        
        # take only summer months measurements to remove bias
        # if ts.month == 6 or ts.month == 7 or ts.month == 8:
            
        if (latitude[idx],longitude[idx]) not in list(hydr_data.keys()):
            hydr_data[(latitude[idx], longitude[idx])] = {'dt': ts.date(), 'temp':[], 'sali':[], 'pres':[], \
                                                          'dept':[]}
            
        hydr_data[(latitude[idx], longitude[idx])]['temp'].append(temperature[idx])
        hydr_data[(latitude[idx], longitude[idx])]['sali'].append(salinity[idx])
        hydr_data[(latitude[idx], longitude[idx])]['pres'].append(pressure[idx])
        hydr_data[(latitude[idx], longitude[idx])]['dept'].append(depth[idx])

    # print(hydr_data)
    return hydr_data

# hydr_data = read_data('data/2008.csv', hydr_data)
# hydr_data = read_data('data/2009.csv', hydr_data)
# hydr_data = read_data('data/2010.csv', hydr_data)
hydr_data = read_data('data/2011.csv', hydr_data)
# hydr_data = read_data('data/2012.csv', hydr_data)
# hydr_data = read_data('data/2013.csv', hydr_data)
# hydr_data = read_data('data/2014.csv', hydr_data)
# hydr_data = read_data('data/2015.csv', hydr_data)
# hydr_data = read_data('data/2016.csv', hydr_data)
# hydr_data = read_data('data/2017.csv', hydr_data)
# hydr_data = read_data('data/2018.csv', hydr_data)
print('Done reading all hydrographic data...')

# df = pd.read_csv('data_500m_extended.csv')
# df['Latitude'] = np.array([round(lat, 5) for lat in df['Latitude'].values])
# df['Longitude'] = np.array([round(lon, 5) for lon in df['Longitude'].values])
# df['Depth'] = np.array([round(d, 5) for d in df['Depth'].values])
# # df['Datetime'] = np.array([datetime(2012, 7, 1) if dt is None else dt for dt in df['Datetime'].values])
# df['Surf_DH'] = np.full(len(df), np.nan)
# print(df.head)

# fp = open('data_500m_extended.csv', 'w+', newline='')
# writer = csv.writer(fp, delimiter=',')
# writer.writerow(['Latitude','Longitude','Depth','Datetime','Surf_DH','D_Siso','hFW'])

Sref=35
Siso=34 # isohaline

############## using specific volume anomaly method 
# for each latlon compute the dynamic height using t, s, p and absolute ssh 
i=0
for latlon,data in hydr_data.items():
               
    if max(data['dept'])>400:

        # interpolate upto 400m depth for each 2m step
        D = [d for d in data['dept'] if d<=400]
        T = [data['temp'][idx] for idx in range(len(D))]
        S = [data['sali'][idx] for idx in range(len(D))]
        P = [data['pres'][idx] for idx in range(len(D))]
        
        if len(D) > 0:
            xdept = np.linspace(0,400,201)
            ytemp = interpolate.interp1d(D, T, fill_value='extrapolate')(xdept)
            ysali = interpolate.interp1d(D, S, fill_value='extrapolate')(xdept)
            ypres = interpolate.interp1d(D, P, fill_value='extrapolate')(xdept)
            
            SA = gsw.SA_from_SP(ysali, ypres, latlon[1], latlon[0])
            CT = gsw.CT_from_t(SA, ytemp, ypres)
            
            try:  
                # Integrate specific volume anomaly to compute dynamic height https://www.teos-10.org/pubs/gsw/html/gsw_geo_strf_dyn_height.html
                dyn_height_anom = gsw.geostrophy.geo_strf_dyn_height(SA, CT, ypres, p_ref=400)
                ster_height = dyn_height_anom[0] / get_g(latlon[0]) 
                    
                # DOT or absolute ssh is the dynamic height at the surface pressure 0dbar
                hydr_data[latlon]['Surf_DH'] = round(ster_height,5)
                print(ster_height)
                # df.loc[(df['Latitude'] == latlon[0]) & (df['Longitude'] == latlon[1]), 'Surf_DH'] = hydr_data[latlon]['Surf_DH']

                # writer.writerow([latlon[0],latlon[1],nearest_depth(latlon[0],latlon[1]), hydr_data[latlon]['dt'],hydr_data[latlon]['Surf_DH'], np.nan, np.nan])
            
            except ValueError as err:
                print('Setting DH to nan due to error...', latlon, ypres)
                i=i+1
                break

    
print(i) 
# print(df.head)
# df.to_csv('2012.csv', index=False)
# df.to_csv('data_500m_extended.csv', index=False)    # for tac 2012 month

# fp.close()
