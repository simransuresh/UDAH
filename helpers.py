import numpy as np
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist
from bathymetric_depth import nearest_depth
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

def dist_pt(lat1, lon1, lat2, lon2):
    lat1 = (lat1 + 90) % 180 - 90
    lat2 = (lat2 + 90) % 180 - 90
    lon1 = (lon1 + 180) % 360 - 180
    lon2 = (lon2 + 180) % 360 - 180
    return geodesic((lat1, lon1), (lat2, lon2)).kilometers

# Compute distances between data points (Cdd) and between data points and target grid points (Cdg)
def dist_mat(t1, t2):
    return cdist(np.radians(t1), np.radians(t2), lambda u, v: 2 * 6371.0 * np.arcsin(
    np.sqrt(np.sin((u[0] - v[0]) / 2)**2 + np.cos(u[0]) * np.cos(v[0]) * np.sin((u[1] - v[1]) / 2)**2) ) )
    
# Define the PV function to handle arrays for efficient calculation
def PV(latd, lond, latg, long):
    # Compute Coriolis force for data points and grid points
    fd = 2 * 7.29e-5 * np.sin(np.deg2rad(latd))
    fg = 2 * 7.29e-5 * np.sin(np.deg2rad(latg))

    # Ensure latitude and longitude pairs are in the correct shape for the nearest depth calculation
    data_coords = np.column_stack((latd, lond))  # Shape should be (N, 2)
    grid_coords = np.column_stack((latg, long))  # Shape should be (M, 2)

    # Compute nearest bathymetric depth for given data or grid points
    Zd = np.array([nearest_depth(data_coords[:, 0][idx], data_coords[:, 1][idx]) for idx in range(len(data_coords[:, 0]))])  # Depth at data points
    Zg = np.array([nearest_depth(grid_coords[:, 0][idx], grid_coords[:, 1][idx]) for idx in range(len(grid_coords[:, 0]))])  # Depth at data points

    # Compute potential vorticity matrix using broadcasting
    PV_matrix = np.abs(fd[:, None]/Zd[:, None] - fg/Zg) / np.sqrt((fd[:, None]/Zd[:, None])**2 + (fg/Zg)**2)

    return PV_matrix

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
