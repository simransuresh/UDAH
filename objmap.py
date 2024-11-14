import numpy as np
from helpers import *
import math
import random
import csv
import pandas as pd
from datetime import datetime
import time

L1 = 600
phi1 = 1
L2 = 300
phi2 = 0.4
T = 60

dp = pd.read_csv('data_points.csv')
# dp = pd.read_csv('results/depth_500m.csv')  # TODO add DH into this csv and use this inside of filt dp csv

hydr_data = {
    (row['Latitude'], row['Longitude']): {
        'depth': row['Depth'],
        'dt': row['Datetime'] if row['Datetime']!='' else np.nan, 
        'dh': float(row['Dynamic_height']) if row['Dynamic_height']!='' else np.nan
    }
    for _, row in dp.iterrows()
}
# print(hydr_data)

def objmap(latg, long, tg):

    # print('#### ----- objective mapping started for the grid.....', latg, long, tg)
    
    ####### get data points and eliminate outliers more than 2SD
    subset1 = [ll for ll, val in hydr_data.items() # take points within radius, within +-3yrs
               if D_mat(ll, target=(latg, long)) <= L1 and \
                tdiff(val['dt'], tg) <= 1095 and \
                math.isnan(val['dh']) is False]    
    Od1 = {ll:hydr_data[ll]['dh'] for ll in subset1}
    Od1 = {ll:val for ll,val in Od1.items() 
            if abs(val - np.mean(list(Od1.values()))) < 2*np.std(list(Od1.values()))}
    subset1 = list(Od1.keys())
    # print(len(subset1), time.time())

    subset2 = [ll for ll, val in hydr_data.items() 
               if D_mat(ll, target=(latg, long)) <= L2 and \
                tdiff(val['dt'], tg) <= 1095 and \
                math.isnan(val['dh']) is False]    
    Od2 = {ll:hydr_data[ll]['dh'] for ll in subset2}
    Od2 = {ll:val for ll,val in Od2.items() 
            if abs(val - np.mean(list(Od2.values()))) < 2*np.std(list(Od2.values()))}
    subset2 = list(Od2.keys())
    # print(len(subset2), time.time())
    # print(len(subset2))
    
    ###### pick data points for mapping
    if len(subset1) > 60:
        # pick random 1/3 points here   
        subset = random.sample(subset1, 20)

        # # pick 1/3 points from L=600km radius with weight to closest points, sort to closest dist and select 1/3 points
        dist1 = [(ll, D_mat(ll, target=(latg, long))) for ll in list(Od1.keys())]
        sorted1 = sorted(dist1, key=lambda x: x[1])
        sorted1 = [point for point, distance in sorted1 if point not in subset]
        subset = subset + sorted1[0:20]

        # # pick remaining 1/3 points from L=300km radius 
        dist2 = [(ll, D_mat(ll, target=(latg, long))) for ll in list(Od2.keys())]   # TODO <60 points, if less available OK?
        sorted2 = sorted(dist2, key=lambda x: x[1])
        sorted2 = [point for point, distance in sorted2 if point not in subset]
        subset = subset + sorted2[0:20]

        # print(len(subset))
    else:
        subset = subset1

    # print('final subset', len(subset), time.time())
    
    Od = np.array([hydr_data[ll]['dh'] for ll in subset])
    n = len(Od)
    
    if len(Od)==0 or Od is None:
        return np.nan, np.nan
    
    # print('######### STAGE 1 ############')
    # print('Lphi---------', L1, phi1)

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
    
    mapped_height = np.linalg.solve(Cdd, Od-Od_mean) 
    
    Og1 = np.dot(Cdg.T, mapped_height) + Od_mean  
        
    Og1_error = np.sqrt(signal_variance - Cdg.T@Cdd@Cdg)
    
    # print("Og1:", Og1, Og1_error) 


    ######## stage2
    # print('######### STAGE 2 ############')
    # print('Lphi---------', L2, phi2, T)
    
    Od_t = [hydr_data[item]['dt'] for item in subset]
    Od = Od-Og1
    signal_variance = signal(Od, len(Od))
    noise_variance = noise(subset, Od, len(Od)) 
    
    tdiff_dd = tdiff(Od_t) #days between dd matrix 
    tdiff_dg = np.array([[d] for d in tdiff(Od_t, tg)])[:,0] #days between dg matrix
    # print(tdiff_dd.shape, tdiff_dg.shape)
    
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

# gp = pd.read_csv('grid_points.csv')
# t = '2011-01-01' #datetime(2011, 1, 1)

# for idx, row in gp.iterrows():
#     (Og, Og_err) = objmap(row['Latitude'], row['Longitude'], '2011-01-01')
#     print(idx, time.time(), Og, Og_err)
