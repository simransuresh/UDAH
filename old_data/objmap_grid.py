import numpy as np
from old.compute_dyn_heights import hydr_data
from old.helpers import signal, noise, distance, PV, covariance, dist
import math
import random

L1 = 600
phi1 = 1
L2 = 300
phi2 = 0.4
T = 60


def objmap(latg, long):

    print('#### ----- objective mapping started for the grid.....', latg, long)
    
    ####### get data points and eliminate outliers
    subset1 = [ll for ll, val in hydr_data.items() if dist(ll[0], ll[1], latg, long) <= L1 and \
            math.isnan(val['ssh']) is False]
    Od1 = {ll:hydr_data[ll]['ssh'] for ll in subset1}
    Od1 = {ll:val for ll,val in Od1.items() 
            if abs(val - np.mean(list(Od1.values()))) < 2*np.std(list(Od1.values()))}
    subset1 = list(Od1.keys())
    print(len(subset1))

    subset2 = [ll for ll, val in hydr_data.items() if dist(ll[0], ll[1], latg, long) <= L2 and \
            math.isnan(val['ssh']) is False]
    Od2 = {ll:hydr_data[ll]['ssh'] for ll in subset2}
    Od2 = {ll:val for ll,val in Od2.items() 
            if abs(val - np.mean(list(Od2.values()))) < 2*np.std(list(Od2.values()))}
    subset2 = list(Od2.keys())
    print(len(subset2))
    
    ###### pick data points for mapping
    if len(subset1) > 60:
        # pick random 1/3 points here   
        subset = random.sample(subset1, 20)

        # pick 1/3 points from L=600km radius with weight to closest points, sort to closest dist and select 1/3 points
        dist1 = [(ll, dist(ll[0], ll[1], latg, long)) for ll in list(Od1.keys())]
        sorted1 = sorted(dist1, key=lambda x: x[1])
        sorted1 = [point for point, distance in sorted1 if point not in subset]
        subset = subset + sorted1[0:20]

        # pick remaining 1/3 points from L=300km radius 
        dist2 = [(ll, dist(ll[0], ll[1], latg, long)) for ll in list(Od2.keys())]
        sorted2 = sorted(dist2, key=lambda x: x[1])
        sorted2 = [point for point, distance in sorted2 if point not in subset]
        subset = subset + sorted2[0:20]
        # print(len(subset))
    else:
        subset = subset1

    Od = np.array([hydr_data[ll]['ssh'] for ll in subset])
    n = len(Od)
    
    if len(Od)==0 or Od is None:
        return np.nan
    
    print('######### STAGE 1 ############')
    print('Lphi---------', L1, phi1)

    signal_variance = signal(Od, n)
    noise_variance = noise(subset, Od, n) 

    distances_dd = distance(subset, subset)
    distances_dg = distance(subset, [(latg, long)])

    given_lats = [lat for lat, lon in subset]
    given_lons = [lon for lat, lon in subset]

    pv_values_dd = PV(given_lats, given_lons, given_lats, given_lons)
    pv_values_dg = PV(given_lats, given_lons, latg, long)

    Cdd = covariance(distances_dd, pv_values_dd, signal_variance, L1, phi1)
    Cdg = covariance(distances_dg, pv_values_dg, signal_variance, L1, phi1)

    Cdd += noise_variance * np.eye(Cdd.shape[0])
    Od_mean = sum(Cdd@Od) / sum(sum(Cdd))

    mapped_height = np.linalg.solve(Cdd, Od-Od_mean) 
    Og1 = np.dot(Cdg.T, mapped_height) + Od_mean  

    print("Og1:", Og1) 


    ######## stage2
    print('######### STAGE 2 ############')
    print('Lphi---------', L2, phi2)

    Od = Od-Og1[0]
    signal_variance = signal(Od, len(Od))
    noise_variance = noise(subset, Od, len(Od)) 

    Cdd = covariance(distances_dd, pv_values_dd, signal_variance, L2, phi2)
    Cdg = covariance(distances_dg, pv_values_dg, signal_variance, L2, phi2)

    Cdd += noise_variance * np.eye(Cdd.shape[0])
    Od_mean = sum(Cdd@Od) / sum(sum(Cdd))

    mapped_height = np.linalg.solve(Cdd, Od-Od_mean)  
    Og2 = np.dot(Cdg.T, mapped_height) + Od_mean 

    print("Og2:", Og2)


    ###### final 
    print('######### FINAL STAGE ############')

    Og = Og1+Og2
    print('Dynamic height', Og[0])
    
    return Og[0]

