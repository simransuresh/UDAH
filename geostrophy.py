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
            
                if len(nearest_coords) == 4:
                    dx_plus = dx_minus = dy_plus = dy_minus = np.nan

                    for neigh in range(len(nearest_coords)):
                        nidx = nearest_coords[neigh][0]
                        nx = df['X_meters'][nidx]
                        ny = df['Y_meters'][nidx]
                    
                        if np.abs(nx-x)==0: # grad along y  
                            if ny-y > 49000 and ny-y <51000:
                                dy_plus = df['Surf_DH'][nidx]
                            elif ny-y > -51000 and ny-y < -49000:
                                dy_minus = df['Surf_DH'][nidx]

                        if np.abs(ny-y)==0: # grad along x
                            if nx-x > 49000 and nx-x <51000:
                                dx_plus = df['Surf_DH'][nidx]
                            elif nx-x > -51000 and nx-x < -49000:
                                dx_minus = df['Surf_DH'][nidx]

                    if dy_plus!=np.nan and dy_minus!=np.nan and dx_plus!=np.nan and dx_minus!=np.nan:
                        grad_x = (dx_plus - dx_minus) / (2*50000)
                        grad_y = (dy_plus - dy_minus) / (2*50000)
                        
                        g = get_g(lat)
                        f = coriolis(lat)

                        ug = -g*grad_y/f
                        vg = g*grad_x/f
                        
                        df.loc[(df['X_meters']==x) & (df['Y_meters']==y), 'ug'] = ug
                        df.loc[(df['X_meters']==x) & (df['Y_meters']==y), 'vg'] = vg

        
    return df
         

for year in range(1993, 2011):
    for i in range(1, 13):
        df = pd.read_csv(f'results/dh/grd_dh_{year}.csv')
        if len(str(i))==1: i='0'+str(i)
        df = df[df['Datetime']==f'{year}-{i}-01']
        df['ug'] = np.full(len(df), np.nan)
        df['vg'] = np.full(len(df), np.nan)
        df = geostrophy(df)
        df.dropna().to_csv(f'interm/grd_gsc_{year}_{i}.csv', index=False)


