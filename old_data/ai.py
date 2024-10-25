import numpy as np
from scipy.spatial.distance import cdist
from scipy.spatial import cKDTree
from bathymetric_depth import nearest_depth

# chosen observation points
subset = [(77.5101, -143.511), (78.6109, -146.961), (77.7388, -141.3697), (77.9923, -134.4294), (77.6371, -143.087), (79.1524, -141.205), (78.4549, -134.889), (74.3868, -145.164), (75.5046, -137.322), (74.9353, -145.461), (76.9131, -143.962), (76.7978, -133.3666), (73.9753, -151.221), (79.1304, -141.008), (78.9887, -140.723), (76.7022, -138.99), (76.8393, -140.263), (78.4866, -146.9573), (78.2498, -141.0715), (78.9082, -135.899), (75.0118, -150.041), (74.9232, -150.064), (75.6289, -151.418), (74.6709, -150.284), (74.671, -150.07), (74.6186, -150.944), (74.7199, -149.402), (75.3178, -153.236), (74.5934, -151.047), (74.7173, -149.228), (74.7306, -148.693), (75.9822, -149.879), (74.8799, -148.133), (74.7101, -148.566), (74.4678, -151.229), (75.9053, -152.731), (74.4181, -151.309), (74.3596, -151.502), (74.3481, -151.556), (74.6516, -147.976), (74.2884, -151.625), (74.6236, -147.778), (74.2657, -151.752), (76.1953, -149.2054), (76.2037, -149.2889), (76.1978, -149.2263), (76.2105, -149.3494), (76.2075, -149.3093), (76.2041, -149.2356), (74.2985, -152.528), (74.2384, -152.032), (76.312, -149.9748), (74.2067, -151.843), (74.2896, -152.675), (74.1991, -152.009), (74.1807, -152.068), (76.3696, -150.1985), (74.2877, -153.05), (76.3422, -149.2951), (76.3411, -149.279)]

# Dynamic height observations at each location
Od = np.array([0.63785078, 0.62819838, 0.6413768, 0.62900474, 0.66928636, 0.36042513,
 0.43502781, 0.83177066, 0.77630484, 0.80460238, 0.73335894, 0.62941059,
 0.6508667, 0.41776979, 0.62302259, 0.67605355, 0.69503378, 0.64214486,
 0.60650402, 0.4141075, 0.77601643, 0.77676383, 0.774071, 0.84156879,
 0.87124863, 0.88107219, 0.83985653, 0.76416769, 0.88403674, 0.86970568,
 0.8768749, 0.75960056, 0.77200857, 0.8647356, 0.89103819, 0.77587633,
 0.94021691, 0.88999367, 0.90748607, 0.75012811, 0.89044561, 0.71676178,
 0.89583259, 0.78228875, 0.7815794, 0.77067151, 0.77500249, 0.77347807,
 0.77072259, 0.87591546, 0.9317338, 0.77990246, 0.88008698, 0.88099138,
 0.89216285, 0.9043266, 0.76826581, 0.87315601, 0.7599587, 0.76241723]
)

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

# Compute signal and noise variances
signal_variance = signal(Od, len(Od))
noise_variance = noise(subset, Od, len(Od)) 

# Define target grid points (latitudes, longitudes)
# target_coords = [(78.0, -145.0), (77.0, -140.0)] 
target_coords = [(75.25, -150.75)] 

# Distance scale for covariance calculation (in km) for stage 1
distance_radius = 600  # Distance radius (Rabe et al., 2011)
pv_scale = 1.0  # Potential vorticity scale

# Function to compute covariance based on distance and PV
def covariance(distance, pv_diff, signal_variance, distance_radius, pv_scale):
    return signal_variance * np.exp(-(distance / distance_radius)**2 - (pv_diff / pv_scale)**2)

# Compute distances between data points (Cdd) and between data points and target grid points (Cdg)
def distance(t1, t2):
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
    Zd = np.array([-nearest_depth(data_coords[:, 0][idx], data_coords[:, 1][idx]) for idx in range(len(data_coords[:, 0]))])  # Depth at data points
    Zg = np.array([-nearest_depth(grid_coords[:, 0][idx], grid_coords[:, 1][idx]) for idx in range(len(grid_coords[:, 0]))])  # Depth at data points

    # Compute potential vorticity matrix using broadcasting
    PV_matrix = np.abs(fd[:, None]/Zd[:, None] - fg/Zg) / np.sqrt((fd[:, None]/Zd[:, None])**2 + (fg/Zg)**2)

    return PV_matrix
   
distances_dd = distance(subset, subset)
distances_dg = distance(subset, target_coords)
# print(distances_dd, distances_dg)

given_lats = [lat for lat, lon in subset]
given_lons = [lon for lat, lon in subset]
target_lats = [lat for lat, lon in target_coords]
target_lons = [lon for lat, lon in target_coords]

pv_values_dd = PV(given_lats, given_lons, given_lats, given_lons)
pv_values_dg = PV(given_lats, given_lons, target_lats, target_lons)
# print(pv_values_dd, pv_values_dg)

# Compute covariance matrices using the covariance function
Cdd = covariance(distances_dd, pv_values_dd, signal_variance, distance_radius, pv_scale)
Cdg = covariance(distances_dg, pv_values_dg, signal_variance, distance_radius, pv_scale)

# Adding noise variance to the diagonal of Cdd for stability
Cdd += noise_variance * np.eye(Cdd.shape[0])
Od_mean = sum(Cdd@Od) / sum(sum(Cdd))
# print(Od_mean)
# print(signal_variance, noise_variance)

# Objective Mapping Calculation Od_mean + Cdg*Cdd*(Od-Od_mean)
# We will estimate the mapped dynamic height on the target grid using the formula: X_target = Cdg.T @ inv(Cdd) @ X_obs
mapped_height = np.linalg.solve(Cdd, Od-Od_mean)  # Solve for Cdd^-1 * dynamic_height
Og1 = np.dot(Cdg.T, mapped_height) + Od_mean  # Compute the interpolated values on the target grid

# Og1_err = np.sqrt( signal_variance - Cdg.T@Cdd@Cdg + (1 - sum(Cdg.T@Cdd))**2/sum(sum(Cdd)) )

print("Og1:", Og1)    # [0.62236907 0.68151931], with -Od_mean, [0.61949844 0.66202294]
# print('Og1 error', Og1_error)

# Cdd = Cdd + noise_variance * np.eye(Cdd.shape[0])  # [0.76995131 0.77012729]

######## stage2
Od = Od-Og1[0]
signal_variance = signal(Od, len(Od))
noise_variance = noise(subset, Od, len(Od)) 

distance_radius = 300  # Distance radius (Rabe et al., 2011)
pv_scale = 0.4  # Potential vorticity scale

# Compute covariance matrices using the covariance function
Cdd = covariance(distances_dd, pv_values_dd, signal_variance, distance_radius, pv_scale)
Cdg = covariance(distances_dg, pv_values_dg, signal_variance, distance_radius, pv_scale)

# Adding noise variance to the diagonal of Cdd for stability
Cdd += noise_variance * np.eye(Cdd.shape[0])
Od_mean = sum(Cdd@Od) / sum(sum(Cdd))
# print(Od_mean)
# print(signal_variance, noise_variance)

# Objective Mapping Calculation Od_mean + Cdg*Cdd*(Od-Od_mean)
# We will estimate the mapped dynamic height on the target grid using the formula: X_target = Cdg.T @ inv(Cdd) @ X_obs
mapped_height = np.linalg.solve(Cdd, Od-Od_mean)  # Solve for Cdd^-1 * dynamic_height
Og2 = np.dot(Cdg.T, mapped_height) + Od_mean  # Compute the interpolated values on the target grid

# Og2_err = np.sqrt( signal_variance - Cdg.T@Cdd@Cdg + (1 - sum(Cdg.T@Cdd))**2/sum(sum(Cdd)) )

print("Og2:", Og2)
# print('Og2 error', Og2_error)

# final 
Og = Og1+Og2
print('Dynamic height', Og)

