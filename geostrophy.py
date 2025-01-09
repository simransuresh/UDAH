import pandas as pd
from helpers import *
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from datetime import datetime

### Compute gsc now
def geostrophy(df):
    for idx, row in df.iterrows():
        lat = row['Latitude']
        lon = row['Longitude']
        x = row['X_meters']
        y = row['Y_meters']
        
        if row['Surf_DH'] !=np.nan or row['Surf_DH']!='nan':
             
            dx_plus_idx = df.loc[ (df['X_meters']==x+50000) & (df['Y_meters']==y) ].index
            dx_minus_idx = df.loc[ (df['X_meters']==x-50000) & (df['Y_meters']==y) ].index

            dy_plus_idx = df.loc[ (df['X_meters']==x) & (df['Y_meters']==y+50000) ].index
            dy_minus_idx = df.loc[ (df['X_meters']==x) & (df['Y_meters']==y-50000) ].index

            if dx_plus_idx.size > 0 and dx_minus_idx.size > 0 and dy_plus_idx.size > 0 and dy_minus_idx.size > 0:

                # print('from conventional')
                grad_x = (df['Surf_DH'][dx_plus_idx[0]] - df['Surf_DH'][dx_minus_idx[0]]) / (2*50000)
                grad_y = (df['Surf_DH'][dy_plus_idx[0]] - df['Surf_DH'][dy_minus_idx[0]]) / (2*50000)

                g = get_g(lat)
                f = coriolis(lat)

                ug = -g*grad_y/f
                vg = g*grad_x/f
                
                df.loc[(df['X_meters']==x) & (df['Y_meters']==y), 'ug'] = ug
                df.loc[(df['X_meters']==x) & (df['Y_meters']==y), 'vg'] = vg
            
            # if there is no exact 50km distant neighbour for computing gsc, take nearest 4 neighbours
            else:
                # takes time also applies to some coordinates along north pole only so including it 
                nearest_coords = find_nearest_coords(lat, lon, df)
            # # print(lat, lon)
            
                if len(nearest_coords) == 4:
                    dx_plus = dx_minus = dy_plus = dy_minus = np.nan
                    # print('!! from nearest neighbours', idx) # boundary points

                    for neigh in range(len(nearest_coords)):
                        nidx = nearest_coords[neigh][0]
                        nx = df['X_meters'][nidx]
                        ny = df['Y_meters'][nidx]
                    
                        if np.abs(nx-x)==0: # grad along y  # TODO remove this 
                            # print('grad_y, x same:', x, y, nx, ny)
                            if ny-y > 49000 and ny-y <51000:
                                dy_plus = df['Surf_DH'][nidx]
                                # print('north', df['Latitude'][nidx], df['Longitude'][nidx])
                            elif ny-y > -51000 and ny-y < -49000:
                                dy_minus = df['Surf_DH'][nidx]
                                # print('south', df['Latitude'][nidx], df['Longitude'][nidx])

                        if np.abs(ny-y)==0: # grad along x
                            # print('grad_x, y same:', x, y, nx, ny)
                            if nx-x > 49000 and nx-x <51000:
                                dx_plus = df['Surf_DH'][nidx]
                                # print('east', df['Latitude'][nidx], df['Longitude'][nidx])
                            elif nx-x > -51000 and nx-x < -49000:
                                dx_minus = df['Surf_DH'][nidx]
                                # print('west', df['Latitude'][nidx], df['Longitude'][nidx])

                    # if neigh==3:
                    if dy_plus!=np.nan and dy_minus!=np.nan and dx_plus!=np.nan and dx_minus!=np.nan:
                        # print(dx_plus, dx_minus, dy_plus, dy_minus)
                        grad_x = (dx_plus - dx_minus) / (2*50000)
                        grad_y = (dy_plus - dy_minus) / (2*50000)
                        # print(grad_x, grad_y)
                        
                        g = get_g(lat)
                        f = coriolis(lat)
                        # print(grad_x, grad_y)

                        ug = -g*grad_y/f
                        vg = g*grad_x/f
                        # print(ug, vg)
                        
                        df.loc[(df['X_meters']==x) & (df['Y_meters']==y), 'ug'] = ug
                        df.loc[(df['X_meters']==x) & (df['Y_meters']==y), 'vg'] = vg
#         # print(lat, lon, row['X_meters'], row['Y_meters'], nearest_coords)

        
    return df
         
 
def gsc(data):
    dh_dx = np.full(len(data), np.nan)
    dh_dy = np.full(len(data), np.nan)

    # Compute gradients
    for i in range(1, len(data) - 1):
        if not np.isnan(data['Surf_DH'].iloc[i]):
            # dh/dx
            lon_plus = data[(data['Longitude'] == data['Longitude'].iloc[i] + 0.75) &
                            (data['Latitude'] == data['Latitude'].iloc[i])]
            lon_minus = data[(data['Longitude'] == data['Longitude'].iloc[i] - 0.75) &
                            (data['Latitude'] == data['Latitude'].iloc[i])]
            if not lon_plus.empty and not lon_minus.empty:
                dh_dx[i] = (lon_plus['Surf_DH'].values[0] -
                            lon_minus['Surf_DH'].values[0]) / (np.deg2rad(0.75 * 2))

            # dh/dy
            lat_plus = data[(data['Longitude'] == data['Longitude'].iloc[i]) &
                            (data['Latitude'] == data['Latitude'].iloc[i] + 0.25)]
            lat_minus = data[(data['Longitude'] == data['Longitude'].iloc[i]) &
                            (data['Latitude'] == data['Latitude'].iloc[i] - 0.25)]
            if not lat_plus.empty and not lat_minus.empty:
                dh_dy[i] = (lat_plus['Surf_DH'].values[0] -
                            lat_minus['Surf_DH'].values[0]) / (np.deg2rad(0.25 * 2))

    # Add gradients to DataFrame
    # data['dh_dx'] = dh_dx
    # data['dh_dy'] = dh_dy

    # Compute geostrophic currents
    lats = data['Latitude']
    f = 2 * 7.29e-5 * np.sin(np.deg2rad(lats))  # Coriolis parameter
    g = 9.81  # Gravity
    Re = 6371008  # Earth's radius
    data['ug'] = -(g / (f * Re)) * dh_dy
    data['vg'] = (g / (f * Re * np.cos(np.deg2rad(lats)))) * dh_dx
    
    return data

# print(k)

for i in range(1, 13):
    df = pd.read_csv('results/dh/grd_dh_2015.csv')
    # df = pd.read_csv('grd_dh_2011.csv')
    if len(str(i))==1: i='0'+str(i)
    df = df[df['Datetime']==f'2015-{i}-01']
    df['ug'] = np.full(len(df), np.nan)
    df['vg'] = np.full(len(df), np.nan)
    df = geostrophy(df)
    # df = gsc(df)
    df.dropna().to_csv(f'grd_dh_2015_{i}.csv', index=False)


# for year in range(2014, 2019):  # Update to your years of interest
#     file_path = f'./results/dh/grd_dh_{year}.csv'  # Replace with your file path
#     df = pd.read_csv(file_path)
    
#     # Convert Datetime column to datetime type
#     df['Datetime'] = pd.to_datetime(df['Datetime'])
    
#     # Extract unique months
#     months = df['Datetime'].dt.month.unique()
    
#     for month in months:
#         # Filter data for the specific month
#         df_month = df[df['Datetime'].dt.month == month]
        
#         # Initialize columns for geostrophic velocities
#         df_month['ug'] = np.full(len(df_month), np.nan)
#         df_month['vg'] = np.full(len(df_month), np.nan)
        
#         # Compute geostrophic velocities
#         df_month = geostrophy(df_month)
        
#         # Save the result to a new file
#         output_file = f'grd_dh_{year}_{month:02d}.csv'
#         df_month.to_csv(output_file, index=False)
#         print(f"Processed and saved: {output_file}")