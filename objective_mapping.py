import numpy as np
import geopy.distance
import math
# from ..SAGA_analysis.read_sla_currents import lat, lon
from compute_dyn_heights import hydr_data
from bathymetric_depth import lat, lon, z
from scipy.spatial import cKDTree
import random

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

####### distance between the target grid point and all other points
def distance(lat1, lon1, lat2, lon2):
    return geopy.distance.geodesic((lat1, lon1), (lat2, lon2)).km

####### Finding potential data points
# find <Od> spatial mean of ssh over all data points within L, Od array and subset of latlon
subset1 = [ll for ll, val in hydr_data.items() 
                   if distance(ll[0], ll[1], latg, long) <= L1
                   and math.isnan(val['ssh']) is False]
Od1 = {ll:hydr_data[ll]['ssh'] for ll in subset1}
print(len(Od1))

subset2 = [ll for ll, val in hydr_data.items() 
                   if distance(ll[0], ll[1], latg, long) <= L2 
                   and math.isnan(val['ssh']) is False]
Od2 = {ll:hydr_data[ll]['ssh'] for ll in subset2}
print(len(Od2))
# print(subset, len(subset))

####### Outlier elimination
# for the potential data points, neglect those outside 2 SD
Od1 = {ll:val for ll,val in Od1.items() 
       if abs(val - np.mean(list(Od1.values()))) < 2*np.std(list(Od1.values()))}
Od2 = {ll:val for ll,val in Od2.items() 
       if abs(val - np.mean(list(Od2.values()))) < 2*np.std(list(Od2.values()))}

####### Subselection of data points     
# Experiments
# N=60 -> 20/20/20 -> Og = -26.04, 0.72m   => CHOSE THIS 1/3
# N=60 -> 10/25/25 -> Og = -2.98, 0.15m OK
# N=70 -> 10/30/30 -> Og = -8.22, 0.31m
# N=65 -> 5/30/30 -> Og = -2.81, 0.12m OK
# N=75 -> 5/35/35 -> Og = -7.87, 0.23m
# N=85 -> 5/40/40 -> Og = -3.27, 0.12m 

######### stage 1
print('######### STAGE 1 ############')
print('Lphi---------', L1, phi1)

# partition only when data points are more than 60 else take full data as it is
def get_data_points():
    if len(subset1) > 60:
        # pick random 1/3 points here   
        subset = random.sample(subset1, 20)
        # print(len(subset))

        # pick 1/3 points from L=600km radius with weight to closest points 
        # sort ll array based on distance closest and select 1/3 points
        dist1 = [(ll, distance(ll[0], ll[1], latg, long)) for ll in list(Od1.keys())]
        sorted1 = sorted(dist1, key=lambda x: x[1])
        sorted1 = [point for point, distance in sorted1 if point not in subset]
        subset = subset + sorted1[0:20]
        # print(len(subset))

        # pick remaining 1/3 points from L=300km radius 
        dist2 = [(ll, distance(ll[0], ll[1], latg, long)) for ll in list(Od2.keys())]
        sorted2 = sorted(dist2, key=lambda x: x[1])
        sorted2 = [point for point, distance in sorted2 if point not in subset]
        subset = subset + sorted2[0:20]
        print(len(subset))
    else:
        subset = subset1
    return subset

# value at data points
subset = get_data_points()
Od = np.array([hydr_data[ll]['ssh'] for ll in subset])
n = len(Od)
# print('Od:', Od, n)
print('n', n)
print('mean', np.mean(Od))

####### signal variance <s2>
def signal_variance(Od):
    s2 = np.sum([(dval - np.mean(Od))**2 for dval in Od])/n
    print('s2:', s2)
    return s2

s2 = signal_variance(Od)

####### noise variance <n2>
def noise_variance(subset, Od):
    # Finding noise between a data point with its nearest neighbour using k-d tree
    tree = cKDTree(np.array(subset))
    nn = []
    for i, point in enumerate(subset):
        dist, idx = tree.query(point, k=2)  # k=2 because the closest point to a point is itself
        nearest_idx = idx[1]  # idx[0] is the point itself, idx[1] is the nearest neighbor
        nearest_neighbor = subset[nearest_idx]
        nn.append(nearest_neighbor)
        
    n2 = np.sum([(Od[idx] - hydr_data[nn[idx]]['ssh'])**2 for idx in range(len(subset))])/(2*n)
    print('n2:', n2)
    return n2

n2 = noise_variance(subset, Od)

####### potential vorticity, planetary vorticity is f/H param
def potential_vorticity(f1, Z1, f2, Z2):
    return abs(f1/Z1 - f2/Z2) / math.sqrt(f1**2/Z1 + f2**2/Z2)

####### covariance computation
def get_Cdg(subset):
    # distance
    D = [distance(ll[0], ll[1], latg, long) for ll in subset]

    ####### coriolis force of data points and grid point
    latd = [ll[0] for ll in subset]
    lond = [ll[1] for ll in subset]
    fds = 2 * 7.29e-5 * np.sin(np.deg2rad(latd))
    fg = 2 * 7.29e-5 * np.sin(np.deg2rad(latg))

    ####### bathymetric depth given in meters from IBCAO dataset
    points = np.column_stack((lat.flatten(), lon.flatten()))
    tree = cKDTree(points)
    # finding closest depth in km to a grid point (negate depth as it is already negative or else sqrt(neg) will give err)
    def nearest_depth(target_lat, target_lon):
        _, idx = tree.query([target_lat, target_lon])
        return z.flatten()[idx]/1000

    Zg = -nearest_depth(latg, long)
    Zds = [-nearest_depth(latd[idx], lond[idx]) for idx in range(len(latd))]
                                        
    PV = [potential_vorticity(fds[idx], Zds[idx], fg, Zg) for idx in range(len(fds))]

    ####### Cdg
    Cdg = np.array([s2*np.exp(-(D[idx]**2/L1**2 + PV[idx]**2/phi1**2)) for idx in range(len(D))]).T
    # print('Cdg:', Cdg)
    return Cdg, latd, lond, fds, Zds

Cdg, latd, lond, fds, Zds = get_Cdg(subset)
print('Cdg:', Cdg, Cdg.shape)

####### Cdd 
def get_Cdd(n, latd, lond, fds, Zds):
    Cdd = np.zeros((n,n))

    for i1 in range(n):
        for i2 in range(n):
            
            if i1 == i2:
                # D and PV of same points will be zero so s2*exp(0)=s2*1=s2
                Cdd[i1][i2] = s2
            else:
                if i1 < i2:
                    Di1i2 = distance(latd[i1], lond[i1], latd[i2], lond[i2])
                    PVi1i2 = potential_vorticity(fds[i1], Zds[i1], fds[i2], Zds[i2])
                    Cdd[i1][i2] = Cdd[i2][i1] = s2*np.exp(-(Di1i2**2/L1**2 + PVi1i2**2/phi1**2))
    return Cdd
       
Cdd = get_Cdd(n, latd, lond, fds, Zds)     
print('Cdd:', Cdd, Cdd.shape)

####### objective mapping procedure 
def obj_map(s2, n2, Cdg, Cdd, Od):
    Cdd = (Cdd + np.eye(n)*n2)**-1

    # weighted mean
    Od_mean = np.sum(Cdd*Od) / np.sum(Cdd)
    print('Od_mean', Od_mean)

    # weight    TODO check the shape of Ben code
    # W = Cdg @ Cdd
    W = Cdg * Cdd
    print('W:', W, W.shape)
    
    # grid value
    # Og1 = Od_mean + W @ (Od-Od_mean)
    Og1 = Od_mean + W * (Od-Od_mean)
    print('Og1', Og1)

    # # grid error
    # Og1_err = np.sqrt( s2 - Cdg.T@Cdd@Cdg + (1 - np.sum(Cdg*Cdd))**2/np.sum(Cdd) )
    # print('Og_err', Og1_err)
    
    # return Og1, Og1_err
    return 

# Og1, Og1_err = obj_map(s2, n2, Cdg, Cdd, Od)
obj_map(s2, n2, Cdg, Cdd, Od)

####### Fine tuning
# # iterate until the grid value is within max/min of the data values
# Nmax_iter = 10
# N_iter = 0

# while N_iter <= Nmax_iter and Og1 != np.nan and (Og1 > max(Od) or Og1 < min(Od)) :
#     N_iter = N_iter + 1
#     print('Iter#############', N_iter)
    
#     subset = get_data_points()
#     Od = np.array([hydr_data[ll]['ssh'] for ll in subset])
#     n = len(Od)
#     print('n', n)
#     print('mean', np.mean(Od))
    
#     s2 = signal_variance(Od)
#     n2 = noise_variance(subset, Od)
    
#     Cdg, latd, lond, fds, Zds = get_Cdg(subset)
#     Cdd = get_Cdd(n, latd, lond, fds, Zds) 
    
#     Og1, Og1_err = obj_map(s2, n2, Cdg, Cdd, Od)
