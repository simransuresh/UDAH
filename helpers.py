import numpy as np
from scipy.spatial import cKDTree
from geopy.distance import geodesic
import pandas as pd

# signal variance
def signal(Od, n):
    return np.sum([(dval - np.mean(Od))**2 for dval in Od])/n

# noise variance 
def noise(obs_coords, Od, n):
    # Finding noise between a data point with its nearest neighbour using k-d tree
    tree = cKDTree(np.array(obs_coords))
    n2 = 0
    
    for i, point in enumerate(obs_coords):
        dist, idx = tree.query(point, k=2)  # k=2 because the closest point to a point is itself
        nearest_idx = idx[1]  # idx[0] is the point itself, idx[1] is the nearest neighbor
        n2 = n2 + (Od[i] - Od[nearest_idx])**2
        
    return n2/(2*n)

def D_mat(subset, target=None):
    if type(subset)==tuple: # 2 point distance for selection
        return geodesic(subset, target).km
    
    if target is not None:  # between grid and data points
        return np.array([ 
            geodesic(point, target).km for point in subset
        ])

    return np.array([ [     # between data points
        geodesic(subset[i], subset[j]).km for j in range(len(subset)) ] 
        for i in range(len(subset))
    ])
    
def PV_mat(latd, lond, latg, long, depth_info):
    # Compute Coriolis force for data points and grid points
    fd = 2 * 7.29e-5 * np.sin(np.deg2rad(latd))  # Coriolis force at data points
    fg = 2 * 7.29e-5 * np.sin(np.deg2rad(latg))  # Coriolis force at grid points

    # Ensure latitude and longitude pairs are in the correct shape for depth lookup
    data_coords = np.column_stack((latd, lond))  # Shape (N, 2)
    grid_coords = np.column_stack((latg, long))  # Shape (M, 2)

    # Retrieve bathymetric depth for data points and grid point using depth_info
    Zd = np.array([depth_info[(lat, lon)]['depth'] for lat, lon in data_coords])  # Depth at data points, shape (N,)
    Zg = np.array([depth_info[(lat, lon)]['depth'] for lat, lon in grid_coords])  # Depth at grid point(s), shape (M,)

    # Calculate PV matrix between all data points (n x n)
    return np.abs(fd[:, None] / Zd[:, None] - fg / Zg) / np.sqrt((fd[:, None] / Zd[:, None])**2 + (fg / Zg)**2)
 
# Function to compute covariance based on distance and PV
def covar1(D, PV, s2, L, phi):
    return s2 * np.exp(-(D / L)**2 - (PV / phi)**2)

def covar2(D, PV, t, s2, L, phi, tau):
    return s2 * np.exp(-(D / L)**2 - (PV / phi)**2 - (t / tau)**2)

def tdiff(dates1, dates2=None):
    # Convert dates1 to pandas datetime objects
    dates1 = pd.to_datetime(dates1)
    
    # Case 1: Difference between two individual dates
    if dates2 is not None:
        dates2 = pd.to_datetime(dates2)
        # If both are single dates, return the difference
        if isinstance(dates1, pd.Timestamp) and isinstance(dates2, pd.Timestamp):
            return np.abs((dates1 - dates2).days)
        
        # Case 2: Difference between a date array and a single date
        if isinstance(dates1, pd.Series) or isinstance(dates1, pd.DatetimeIndex):
            dates1_np = dates1.to_numpy()
            dates2_np = pd.to_datetime(dates2).to_numpy()  # Convert target date to numpy
            return np.abs((dates1_np - dates2_np).astype('timedelta64[D]').astype(int))

    # Case 3: Difference between dates in an array (matrix form)
    dates1_np = dates1.to_numpy()  # Convert dates array to numpy array
    date_diff_matrix = (dates1_np[:, None] - dates1_np[None, :]).astype('timedelta64[D]').astype(int)
    
    return np.abs(date_diff_matrix)

    
dp = pd.read_csv('results/depth_500m.csv')
depth_info = {
    (row['Latitude'], row['Longitude']): {
        'depth': row['Depth'],
    }
    for _, row in dp.iterrows()
}

# D_mat(subset=[(82.21415, 39.16856), (82.22155, 39.376728), (82.29131, 39.610153)], target=)
# subset=[(82.21415, 39.16856), (82.22155, 39.376728), (82.29131, 39.610153)]
# latg = 82.25
# long = 39.75
# distances_dd = D_mat(subset)
# # print(distances_dd)
# distances_dg = D_mat(subset, target=(latg, long))
# print(distances_dg)
# print(dist_pt(82.21415, 39.16856, 82.22155, 39.376728))

# given_lats = [lat for lat, lon in subset]
# given_lons = [lon for lat, lon in subset]

# pv_dd = PV_mat(given_lats, given_lons, given_lats, given_lons, depth_info)
# pv_dg = PV_mat(given_lats, given_lons, latg, long, depth_info)

# print(pv_dd, pv_dg)
# print(tdiff('2011-01-01', '2012-01-02'))
# print(tdiff(['2011-01-01', '2011-02-01'], '2012-01-02'))
# print(tdiff(['2011-01-01', '2011-02-01', '2012-01-02']))
