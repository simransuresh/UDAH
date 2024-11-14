import numpy as np
from helpers import *
import math
import random
import pandas as pd
import time
import csv

L1 = 600
phi1 = 1
L2 = 300
phi2 = 0.4
T = 60

# first all profiles from dynamic_height.py -> dh_2011_2018_deep.csv
# dp = pd.read_csv('data_points.csv') -> only points within CAO
dp = pd.read_csv('all_500m.csv')  # takes all >500m profiles for mapping grid points within CAO

hydr_data = {
    (row['Latitude'], row['Longitude']): {
        'depth': row['Depth'],
        'dt': row['Datetime'] if row['Datetime']!='' else np.nan, 
        'dh': float(row['Dynamic_height']) if row['Dynamic_height']!='' else np.nan
    }
    for _, row in dp.iterrows()
}

# outlier elimination
mean = np.nanmean([val['dh'] for val in hydr_data.values()])
sd = np.nanstd([val['dh'] for val in hydr_data.values()])
hydr_data = {key:val for key,val in hydr_data.items() if abs(val['dh'] - mean) < 2*sd}

def find_weights(subset, latg, long, tg=None):
    
    D = D_mat(subset, target=(latg, long))
    PV = PV_mat([ll[0] for ll in subset], [ll[1] for ll in subset], latg, long, depth_info)[:,0]
    
    weighting_func = (D/L1)**2 + (PV/phi1)**2
    
    # for stage 2 including t also
    if tg is not None:  
        Od_t = [hydr_data[item]['dt'] for item in subset]
        t = tdiff(Od_t, tg)
        weighting_func += (t/T)**2
       
    weight =  np.exp( -weighting_func )
    weight = {ll: weight[idx] for idx, ll in enumerate(subset)}   # idx of subset should be same of D,PV,t mat idx

    return weight

def objmap(latg, long, tg):

    # print('#### ----- objective mapping started for the grid.....', latg, long, tg)

    ####### get data points 
    subset1 = [ll for ll, val in hydr_data.items() # take points within L1 radius, within +-3yrs
               if D_mat(ll, target=(latg, long)) <= L1 and tdiff(val['dt'], tg) <= 1095 and \
                math.isnan(val['dh']) is False]   
    # print(len(subset1))
    
    ###### pick data points for mapping
    if len(subset1) > 60:
        # pick random 1/3 points here   
        subset = random.sample(subset1, 20)
        rem_subset = [point for point in subset1 if point not in subset]
        
        # 1/3 - from remaining get weights based on D, PV and sorted by highest to lowest and take first 20 points
        w1 = find_weights(rem_subset, latg, long, tg=None)    
        l1_subset = sorted(w1, key=w1.get, reverse=True)[0:20]
        rem_subset1 = [point for point in subset1 if point not in l1_subset and point not in subset]

        # 1/3 - same but used t weight also        
        w2 = find_weights(rem_subset1, latg, long, tg=tg)
        l2_subset = sorted(w2, key=w2.get, reverse=True)[0:20]
        rem_subset2 = [point for point in subset1 if point not in l2_subset and point not in l1_subset and point not in subset]
        
        subset = list(set(subset+l1_subset+l2_subset))
        
    else:
        subset = subset1

    # print('final subset', len(subset), time.time())
    
    #### get value of dh from chosen subset
    Od = np.array([hydr_data[ll]['dh'] for ll in subset])
    n = len(Od)
    
    if len(Od)==0 or Od is None:
        return np.nan, np.nan
    
    # print('######### STAGE 1 ############', L1, phi1)

    signal_variance = signal(Od, n)
    noise_variance = noise(subset, Od, n) 
    # print(signal_variance, noise_variance)
    
    D_dd = D_mat(subset)
    D_dg = D_mat(subset, target=(latg, long))
    # print(D_dd.shape, D_dg.shape)
    
    given_lats = [lat for lat, lon in subset]
    given_lons = [lon for lat, lon in subset]

    PV_dd = PV_mat(given_lats, given_lons, given_lats, given_lons, depth_info)
    PV_dg = PV_mat(given_lats, given_lons, latg, long, depth_info)[:,0] # get first column for one gp
    # print(PV_dd.shape, PV_dg.shape)
    
    Cdd = covar1(D_dd, PV_dd, signal_variance, L1, phi1)
    Cdg = covar1(D_dg, PV_dg, signal_variance, L1, phi1)
    # print(Cdd.shape, Cdg.shape)
    
    Cdd += noise_variance * np.eye(Cdd.shape[0])
    Od_mean = sum(Cdd@Od) / sum(sum(Cdd))
    
    ### mapping procedure
    mapped_height = np.linalg.solve(Cdd, Od-Od_mean) 
    Og1 = np.dot(Cdg.T, mapped_height) + Od_mean  
        
    Og1_error = np.sqrt(signal_variance - Cdg.T@Cdd@Cdg)
    # print("Og1:", Og1, Og1_error) 


    ######## stage2
    # print('######### STAGE 2 ############', L2, phi2, T)
    
    Od = Od-Og1
    
    signal_variance = signal(Od, len(Od))
    noise_variance = noise(subset, Od, len(Od)) 
    
    Od_t = [hydr_data[item]['dt'] for item in subset]
    tdiff_dd = tdiff(Od_t) # days between dd matrix 
    tdiff_dg = tdiff(Od_t, tg) #days between dg matrix
    
    Cdd = covar2(D_dd, PV_dd, tdiff_dd, signal_variance, L2, phi2, T)
    Cdg = covar2(D_dg, PV_dg, tdiff_dg, signal_variance, L2, phi2, T)
    # print(Cdd.shape, Cdg.shape)
    
    Cdd += noise_variance * np.eye(Cdd.shape[0])
    Od_mean = sum(Cdd@Od) / sum(sum(Cdd))

    mapped_height = np.linalg.solve(Cdd, Od-Od_mean)  
    Og2 = np.dot(Cdg.T, mapped_height) + Od_mean 
    
    Og2_error = np.sqrt(signal_variance - Cdg.T@Cdd@Cdg)
    # print("Og2:", Og2, Og2_error) 

    # ###### final 
    # print('######### FINAL STAGE ############')

    Og = Og1+Og2
    Og_error = Og2_error
    # print('Dynamic height', Og, Og_error)
    
    return Og, Og_error

# objmap(latg, long, tg)
# objmap(87.5, 0, '2011-09-01')

gp = pd.read_csv('grid_points.csv')
t = '2011-01-01' #datetime(2011, 1, 1)

# output_file = "results/grd_dh_2011_01_wnn.csv"  # another run with distance sorted poiints selected within L1 L2
output_file = "results/grd_dh_2011_01_final.csv"  # another run with weight function inclusive

fp = open(output_file, mode='w', newline='') 
writer = csv.writer(fp)
writer.writerow(["Datetime", "Latitude", "Longitude", "Depth", "Dynamic_Height", "DH_error"])    

for idx, row in gp.iterrows():
    lat, lon = row['Latitude'], row['Longitude']
    (Og, Og_err) = objmap(lat, lon, t)
    
    print(idx, row['Latitude'], row['Longitude'], time.time(), Og, Og_err)
    writer.writerow([t, lat, lon, depth_info[(lat, lon)]['depth'], Og, Og_err])
    
    # if idx==2:
    #     break

fp.close()
print(f"Results written to {output_file}")
