# import pandas as pd 

# # data = pd.read_csv('60e_2015.csv')
# # lon = data['Longitude_[deg]']

# with open('60E_2015_50m.txt', 'r') as fp:
#     data = fp.readlines()
    
# pressure_ = []
# depth_ = []
# temp_ = []
# salin_ = []

# for d in data:
#     d = d[:-1]
#     d = d.split('\t')
#     # print(d)
    
#     lon = d[0]
#     lat = d[1]
#     pressure_.append(d[2])
#     depth_.append(d[3]) 
#     temp_.append(d[5]) 
#     salin_.append(d[7]) 
    
# # print(lon, lat, pressure_, depth_, temp_, salin_)
    
    
# def tsp2density(T, S, p):
#     T = float(T)
#     S = float(S)
#     p = float(p)
#     # % Density of pure water as a function of the temperature (kg/m3)
#     # % water salinity = 0 ppt
#     rhow = 999.842594 + (6.793952e-2 * T) - (9.095290e-3 * T**2) + (1.001685e-4 * T**3) - (1.120083e-6 * T**4) + (6.536332e-9 * T**5)
    
#     # % Water density (kg/m3) at one standard atmosphere, p = 0. 
#     rhost0 = rhow + (S * (0.824493 - (4.0899e-3 * T) + (7.6438e-5 * T**2) - (8.2467e-7 * T**3) + (5.3875e-9 * T**4))) + ( S**(3/2) * (-5.72466e-3 + (1.0227e-4 * T) - (1.6546e-6 * T**2))) + ( 4.8314e-4 * S**2)
#     # % Water pure secant bulk modulus
#     kw = 19652.21 + (148.4206 * T) - (2.327105 * T**2) + (1.360477e-2 * T**3) - (5.155288e-5 * T**4)
#     # % Water secant bulk modulus at one standard atmosphere (p = 0)
#     kst0 = kw + (S * (54.6746 - (0.603459 * T) + (1.09987e-2 * T**2) -(6.1670e-5 * T**3))) + (S **(3/2) * (7.944e-2 + (1.6483e-2 * T) -(5.3009e-4 * T**2)))
    
#     # % Water secant bulk modulus at pressure values, p
#     kstp = kst0 + \
#         (p*(3.239908 + (1.43713e-3 * T) + (1.16092e-4 * T**2) -(5.77905e-7* T**3))) + \
#         ((p*S) *(2.2838e-3 - (1.0981e-5 * T) - (1.6078e-6 * T**2))) + \
#         (1.91075e-4 * p * S**(3/2)) + \
#         (p**2 * (8.50935e-5 - (6.12293e-6 * T) + (5.2787e-8 * T**2))) + \
#         ((p**2 *S) * (-9.9348e-7 + (2.0816e-8 * T) + (9.1697e-10 * T**2)))
        
#     # % Water density at any pressure (kg/m3)
#     rho = rhost0/(1-(p/kstp))
#     return rho

# rho_ = []
# for i in range(len(pressure_)):
#     rho = tsp2density(temp_[i], salin_[i], pressure_[i])
#     rho_.append(rho)
    
    # print(rho)


# import pandas as pd
# from compute_dyn_heights import hydr_data
# import csv

# fp = open('ll.csv', 'w+', newline='')
# writer = csv.writer(fp, delimiter=',')

# for k, v in hydr_data.items():
#     lat, lon = k[0], k[1]
#     writer.writerow([lat, lon])
    
# fp.close()

# # for each latlon compute the dynamic height using t, s, p and absolute ssh 
# i=0
# for latlon,data in hydr_data.items():
    
#     # Compute Absolute Salinity from practical salinity 
#     SA = gsw.SA_from_SP(data['sali'], data['pres'], latlon[1], latlon[0])

#     # Compute Conservative Temperature from SA and temperature
#     CT = gsw.CT_from_t(SA, data['temp'], data['pres'])

#     try:

#         # Integrate specific volume anomaly to compute dynamic height
#         geopot_height_anom = gsw.geo_strf_dyn_height(SA, CT, data['pres'], p_ref=0)
#         dyn_height_anom = geopot_height_anom / 9.8
#         # print(dynamic_height_anom)
        
#         # DOT or absolute ssh is the dynamic height at the surface pressure 0dbar
#         # abs_ssh = dynamic_height[0]
#         hydr_data[latlon]['ssh'] = dyn_height_anom[-1]
        
#         print(latlon, hydr_data[latlon]['ssh'])
        
#     except ValueError as err:
#         # print(latlon, data['pres'])
#         hydr_data[latlon]['ssh'] = np.nan
#         i=i+1
        
# print(i)    # 5 profiles throw error of pressure not monotonically increasing

# gives pairwise difference between each data point with each other
# n2 = np.sum((Od[:,np.newaxis] - Od[np.newaxis,:])**2)/(2*n)  ??
# n2 = np.sum([ for ])


# import matplotlib.pyplot as plt
# from matplotlib import colors
# from matplotlib.ticker import PercentFormatter
# import numpy as np

# # Create a random number generator with a fixed seed for reproducibility
# rng = np.random.default_rng(19680801)

# N_points = 100000
# n_bins = 20

# # Generate two normal distributions
# dist1 = rng.standard_normal(N_points)

# # fig, axs = plt.plot(1, 2, tight_layout=True)

# # N is the count in each bin, bins is the lower-limit of the bin
# N, bins, patches = plt.hist(dist1, bins=n_bins)

# # We'll color code by height, but you could use any scalar
# fracs = N / N.max()

# # we need to normalize the data to 0..1 for the full range of the colormap
# norm = colors.Normalize(fracs.min(), fracs.max())

# # Now, we'll loop through our objects and set the color of each accordingly
# for thisfrac, thispatch in zip(fracs, patches):
#     color = plt.cm.viridis(norm(thisfrac))
#     thispatch.set_facecolor(color)

# plt.show()

# import numpy as np
# print(np.linspace(0,400,401))

        # REF: https://faculty.washington.edu/luanne/pages/ocean420/notes/Dynamic.pdf
        # try:
        #     DH = 0
        #     for idx in range(len(D)-1):
        #         dp = P[idx] - P[idx+1]
        #         SA = gsw.SA_from_SP(S[idx], P[idx], latlon[1], latlon[0])
        #         CT = gsw.CT_from_t(SA, T[idx], P[idx])
        #         rho = gsw.density.rho(SA, CT, P[idx])
        #         # integral of specific volume anomaly 
        #         DH = DH + (1/rho)*dp

        #     DH = -DH/9.7
        #     print(DH)
        #     hydr_data[latlon]['ssh'] = DH
        
        
# reference TS   
# S_ref = 34.8    # from Aagard paper for AO
# T_ref = 0
########### using Steele paper
# def get_density(P, lat, lon, S, T):
#     # if S is None: S=S_ref
#     # if T is None: T=T_ref
    
#     SA = gsw.SA_from_SP(S, P, lon, lat)
#     CT = gsw.CT_from_t(SA, T, P)
#     rho = gsw.density.rho(SA, CT, P)
    
#     return rho

# # computing dynamic height from density     
# for latlon,data in hydr_data.items():
#     lat = latlon[0]
#     lon = latlon[1]
    
#     # interpolate upto 400m depth for each 2m step
#     D = [d for d in data['dept'] if d<=400]
#     T = [data['temp'][idx] for idx in range(len(D))]
#     S = [data['sali'][idx] for idx in range(len(D))]
#     P = [data['pres'][idx] for idx in range(len(D))]
    
 
#     if len(D) > 0:
#         xdept = np.linspace(0,400,201)
#         ytemp = interpolate.interp1d(D, T, fill_value='extrapolate')(xdept)
#         ysali = interpolate.interp1d(D, S, fill_value='extrapolate')(xdept)
#         ypres = interpolate.interp1d(D, P, fill_value='extrapolate')(xdept)
        
#         D = xdept[::-1]
#         T = ytemp[::-1]
#         S = ysali[::-1]
#         P = ypres[::-1]
        
#         dh = 0
#         for idx in range(len(D)-1):
#             # if P[idx] <= 500: 
#             rho_ref = get_density(P[idx], lat, lon, S=S_ref, T=T_ref)
#             rho = get_density(P[idx], lat, lon, S=S[idx], T=T[idx])
            
#             dz = D[idx] - D[idx+1]
#             dh += (((rho_ref-rho)/rho_ref) * dz)
#         # print(dh)
            
#         hydr_data[latlon]['ssh'] = dh 
        # print(latlon, hydr_data[latlon]['ssh'])

import numpy as np

# def normalize_lat_lon(lat, lon):
#     # Normalize latitude to the range [-90, 90]
#     lat = (lat + 90) % 180 - 90

#     # Normalize longitude to the range [-180, 180]
#     lon = (lon + 180) % 360 - 180
#     return lat, lon

# # Example latitude and longitude values
# latitude = [95, -100, 45, 88, -92]  # Some latitudes are out of range
# longitude = [190, -190, 360, -370, 720]  # Longitudes are also out of range

# # Normalize the latitudes and longitudes
# latitude, longitude = normalize_lat_lon(np.array(latitude), np.array(longitude))

# print("Normalized Latitude:", latitude)
# print("Normalized Longitude:", longitude)


# print(np.linspace(70, 80, 41))
# print(np.linspace(-155, -125, 61))

import numpy as np
import pandas as pd

def date_difference(dates1, dates2=None):
    # Convert dates1 to pandas datetime objects
    dates1 = pd.to_datetime(dates1)
    
    # Case 1: Difference between two individual dates
    if dates2 is not None:
        dates2 = pd.to_datetime(dates2)
        # If both are single dates, return the difference
        if isinstance(dates1, pd.Timestamp) and isinstance(dates2, pd.Timestamp):
            return abs((dates1 - dates2).days)
        
        # Case 2: Difference between a date array and a single date
        if isinstance(dates1, pd.Series) or isinstance(dates1, pd.DatetimeIndex):
            dates1_np = dates1.to_numpy()
            dates2_np = pd.to_datetime(dates2).to_numpy()  # Convert target date to numpy
            return abs((dates1_np - dates2_np).astype('timedelta64[D]').astype(int))

    # Case 3: Difference between dates in an array (matrix form)
    dates1_np = dates1.to_numpy()  # Convert dates array to numpy array
    date_diff_matrix = (dates1_np[:, None] - dates1_np[None, :]).astype('timedelta64[D]').astype(int)
    
    return abs(date_diff_matrix)

# Example usage

# Case 1: Difference between two dates
date1 = '2023-10-15'
date2 = '2023-10-10'
diff_1 = date_difference(date1, date2)
print(f"Difference between {date1} and {date2}: {diff_1} days")

# Case 2: Difference between a date array and a single date
date_array = ['2023-10-01', '2023-10-05', '2023-10-10', '2023-10-20']
target_date = '2023-10-15'
diff_2 = date_difference(date_array, target_date)
print(f"Difference between {date_array} and {target_date}: {diff_2} days")

# Case 3: Difference between dates in an array (matrix form)
diff_3 = date_difference(date_array)
print(f"Matrix of differences between dates in {date_array}:\n{diff_3}")

