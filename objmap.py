import numpy as np
from helpers import signal, noise, dist_mat, dist_pt, PV, covar1, covar2, tdiff, nearest_depth
import math
import random
import csv
# from datetime import datetime

L1 = 600
phi1 = 1
L2 = 300
phi2 = 0.4
T = 60

print('Loading dynamic heights computed using T-S..... ')

hydr_data = {}

# Open and read the CSV file
# with open('results/dh_2011_2018.csv', mode='r') as file:  # Max, Min: 7253503.838336006 -13.164722531269735
with open('results/dh_2011_2018_deep.csv', mode='r') as file:    # Max, Min: 1.9064679533773496 -0.3754665568686381
        reader = csv.DictReader(file)

        # Iterate over each row in the CSV
        for row in reader:
                # Store the datetime and dynamic height in the dict with the key (lat, lon)
                # # nearest_depth(float(row['Latitude']), float(row['Longitude']))>=0.5 # removes 500m isobath/Shelf seas
                hydr_data[(float(row['Latitude']), float(row['Longitude']))] = {
                        'dt': row['Datetime'], 'ssh': float(row['Dynamic_height']) 
                        if row['Dynamic_height']!='' #and 
                        # float(row['Dynamic_height'])<=1.5 and 
                        # float(row['Dynamic_height'])>=0 
                        else np.nan
                }
           
file.close()

def objmap(latg, long, tg):

    print('#### ----- objective mapping started for the grid.....', latg, long, tg)
    
    ####### get data points and eliminate outliers more than 2SD
    subset1 = [ll for ll, val in hydr_data.items() # take points within radius, within +-3yrs
               if dist_pt(ll[0], ll[1], latg, long) <= L1 and \
                tdiff(val['dt'], tg) <= 1095 and \
                math.isnan(val['ssh']) is False]    
    Od1 = {ll:hydr_data[ll]['ssh'] for ll in subset1}
    Od1 = {ll:val for ll,val in Od1.items() 
            if abs(val - np.mean(list(Od1.values()))) < 2*np.std(list(Od1.values()))}
    subset1 = list(Od1.keys())
    # print(len(subset1))

    subset2 = [ll for ll, val in hydr_data.items() 
               if dist_pt(ll[0], ll[1], latg, long) <= L2 and \
                tdiff(val['dt'], tg) <= 1095 and \
                math.isnan(val['ssh']) is False]    
    Od2 = {ll:hydr_data[ll]['ssh'] for ll in subset2}
    Od2 = {ll:val for ll,val in Od2.items() 
            if abs(val - np.mean(list(Od2.values()))) < 2*np.std(list(Od2.values()))}
    subset2 = list(Od2.keys())
    # print(len(subset2))
    
    ###### pick data points for mapping
    if len(subset1) > 60:
        # pick random 1/3 points here   
        subset = random.sample(subset1, 20)

        # pick 1/3 points from L=600km radius with weight to closest points, sort to closest dist and select 1/3 points
        dist1 = [(ll, dist_pt(ll[0], ll[1], latg, long)) for ll in list(Od1.keys())]
        sorted1 = sorted(dist1, key=lambda x: x[1])
        sorted1 = [point for point, distance in sorted1 if point not in subset]
        subset = subset + sorted1[0:20]

        # pick remaining 1/3 points from L=300km radius 
        dist2 = [(ll, dist_pt(ll[0], ll[1], latg, long)) for ll in list(Od2.keys())]
        sorted2 = sorted(dist2, key=lambda x: x[1])
        sorted2 = [point for point, distance in sorted2 if point not in subset]
        subset = subset + sorted2[0:20]
        # print(len(subset))
    else:
        subset = subset1

    Od = np.array([hydr_data[ll]['ssh'] for ll in subset])
    n = len(Od)
    
    if len(Od)==0 or Od is None:
        return np.nan, np.nan
    
    # print('######### STAGE 1 ############')
    # print('Lphi---------', L1, phi1)

    signal_variance = signal(Od, n)
    noise_variance = noise(subset, Od, n) 

    distances_dd = dist_mat(subset, subset)
    distances_dg = dist_mat(subset, [(latg, long)])
    
    given_lats = [lat for lat, lon in subset]
    given_lons = [lon for lat, lon in subset]

    pv_values_dd = PV(given_lats, given_lons, given_lats, given_lons)
    pv_values_dg = PV(given_lats, given_lons, latg, long)
    
    Cdd = covar1(distances_dd, pv_values_dd, signal_variance, L1, phi1)
    Cdg = covar1(distances_dg, pv_values_dg, signal_variance, L1, phi1)
    
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
    Od = Od-Og1[0]
    signal_variance = signal(Od, len(Od))
    noise_variance = noise(subset, Od, len(Od)) 
    
    tdiff_dd = tdiff(Od_t) #days between dd matrix 
    tdiff_dg = np.array([[d] for d in tdiff(Od_t, tg)]) #days between dg matrix
    
    Cdd = covar2(distances_dd, pv_values_dd, tdiff_dd, signal_variance, L2, phi2, T)
    Cdg = covar2(distances_dg, pv_values_dg, tdiff_dg, signal_variance, L2, phi2, T)
    
    Cdd += noise_variance * np.eye(Cdd.shape[0])
    
    Od_mean = sum(Cdd@Od) / sum(sum(Cdd))

    mapped_height = np.linalg.solve(Cdd, Od-Od_mean)  
    
    Og2 = np.dot(Cdg.T, mapped_height) + Od_mean 
    
    Og2_error = np.sqrt(signal_variance - Cdg.T@Cdd@Cdg)
    
    # print("Og2:", Og2, Og2_error) 

    ###### final 
    # print('######### FINAL STAGE ############')

    Og = (Og1+Og2)[0]
    Og_error = Og2_error[0][0] 
    # print('Dynamic height', Og, Og_error)
    
    return Og, Og_error

# objmap(latg, long, tg)
# objmap(87.5, 0, '2011-09-01')
