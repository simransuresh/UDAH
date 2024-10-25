import numpy as np
import geopy.distance
import math
# from ..SAGA_analysis.read_sla_currents import lat, lon
import random
from distance import dist
from get_data_points import get_data_points, outlier_elimination
from old.compute_dyn_heights import hydr_data
from pick_data_points import pick_data_points
from variances import signal_variance, noise_variance
from potential_vorticity import PV

# apply objective mapping to regular grid (eg: SAGA)
# reg_lats = np.linspace(60, 90, nlat)
# reg_lons = np.linspace(-180, 180, nlon)

####### Constants
# for each point in reg grid, from the keys of hydr_data find haversine distance 
# if within L, then take it for <Od>
latg = 75.25
long = -150.75
L1 = 600
phi1 = 1
L2 = 300
phi2 = 0.4
T = 60

subset1, Od1, subset2, Od2 = get_data_points(hydr_data, latg, long, L1, L2) 
# Od1, Od2 = outlier_elimination(Od1, Od2)

####### stage 1
print('######### STAGE 1 ############')
print('Lphi---------', L1, phi1)

subset, Od, n = pick_data_points(hydr_data, subset1, latg, long, Od1, Od2) 
print('Chosen points:', subset)

s2 = signal_variance(Od, n)
n2 = noise_variance(hydr_data, subset, Od, n)

####### covariance computation
def get_Cdg(subset, latg, long):
    # distance
    Dij = np.array([dist(ll[0], ll[1], latg, long) for ll in subset])
    PVij = np.array([PV(ll[0], ll[1], latg, long) for ll in subset])
    
    ####### Cdg
    Cdg = np.array(s2*np.exp(-(Dij**2/L1**2 + PVij**2/phi1**2))).T
    return Cdg

####### Cdd 
def get_Cdd(n, subset):
    Cdd = np.zeros((n,n))
    latd = [ll[0] for ll in subset]
    lond = [ll[1] for ll in subset]
    
    for i1 in range(n):
        for i2 in range(n):
            
            Dij = dist(latd[i1], lond[i1], latd[i2], lond[i2])
            PVij = PV(latd[i1], lond[i1], latd[i2], lond[i2])  
            Cdd[i1][i2] = s2*np.exp(-(Dij**2/L1**2 + PVij**2/phi1**2))  
              
    return Cdd

# Cdg = get_Cdg(subset, latg, long)
# print('Cdg:', Cdg, Cdg.shape)

# Cdd = get_Cdd(n, subset)     
# print('Cdd:', Cdd, Cdd.shape)

####### objective mapping procedure 
def obj_map(s2, n2, Od, Cdg, Cdd):

    Cdd = (Cdd + np.eye(n)*n2)**-1

    # weighted mean
    Od_mean = sum(Cdd@Od) / sum(sum(Cdd)) 
    print('Od_mean', Od_mean)
    print('Od', Od)
    print('mean-Od:', Od-Od_mean)

    # weight   
    W = Cdg @ Cdd
    # W = Cdg * Cdd
    print('W:', W, W.shape)
    
    print('weight after mean-Od:', W @(Od-Od_mean))
    
    # grid value 
    Og1 = Od_mean + W @ (Od-Od_mean)
    # Og1 = Od_mean + W * (Od-Od_mean)
    print('Og1', Og1)

    # grid error 
    Og1_err = np.sqrt( s2 - Cdg.T@Cdd@Cdg + (1 - sum(Cdg@Cdd))**2/sum(sum(Cdd)) )
    # Og1_err = s2 - W @ Cdg.T
    
    print('Og_err', Og1_err)
    
    return Og1, Og1_err

# Og1, Og1_err = obj_map(s2, n2, Od, Cdg, Cdd)
# obj_map(s2, n2, Cdg, Cdd, Od)

####### Fine tuning 
# iterate until the grid value is within max/min of the data values
# Nmax_iter = 10
# N_iter = 0

# while N_iter <= Nmax_iter and Og1 != np.nan and (Og1 > max(Od) or Og1 < min(Od)) :
#     N_iter = N_iter + 1
#     print('Iter#############', N_iter)
    
#     # add rand*10^-4 to lat, lon and depth
#     latg = latg + random.randint(1, 10)*10**-4
#     long = long + random.randint(1, 10)*10**-4
    
#     Cdg = get_Cdg(subset, latg, long)
#     # print('Cdg:', Cdg, Cdg.shape)

#     Og1, Og1_err = obj_map(s2, n2, Od, Cdg, Cdd)




