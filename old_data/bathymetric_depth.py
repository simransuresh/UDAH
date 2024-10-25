import netCDF4 as nc
import numpy as np
from pyproj import Proj, transform, Transformer
import xarray as xr
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt

 ####### bathymetric depth given in meters from IBCAO dataset
# Open the original NetCDF file
ds = nc.Dataset('data/IBCAO_v4_2_13_400m_ice.nc', 'r')

x = ds.variables['x'][:]
y = ds.variables['y'][:]
z = ds.variables['z'][:]    # in meters

ds.close()

lat_ts = 60.0  # Standard parallel (latitude of true scale)
lon_0 = 0.0  # Central meridian
proj_polar = Proj(proj='stere', lat_ts=lat_ts, lon_0=lon_0, lat_0=90, ellps='WGS84')
proj_wgs84 = Proj(proj='latlong', ellps='WGS84')

transformer = Transformer.from_proj(proj_polar, proj_wgs84)
xx, yy = np.meshgrid(x, y)
lon, lat = transformer.transform(xx, yy)

# form a KD tree which gives nearest depth z to a lat, lon
points = np.column_stack((lat.flatten(), lon.flatten()))
tree = cKDTree(points)

# finding closest depth in km to a grid point (negate depth as it is already negative or else 
# sqrt(neg) will give err)
def nearest_depth(target_lat, target_lon):
    _, idx = tree.query([target_lat, target_lon])
    return z.flatten()[idx]/1000
