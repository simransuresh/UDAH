import netCDF4 as nc
import numpy as np
from pyproj import Proj, transform, Transformer
import xarray as xr
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from scipy.interpolate import griddata
import pandas as pd

# Open the original NetCDF file
ds = nc.Dataset('IBCAO_v4_2_13_400m_ice.nc', 'r')

x = ds.variables['x'][:]
y = ds.variables['y'][:]
z = ds.variables['z'][:]

ds.close()

lat_ts = 60.0  # Standard parallel (latitude of true scale)
lon_0 = 0.0  # Central meridian
proj_polar = Proj(proj='stere', lat_ts=lat_ts, lon_0=lon_0, lat_0=90, ellps='WGS84')
proj_wgs84 = Proj(proj='latlong', ellps='WGS84')

transformer = Transformer.from_proj(proj_polar, proj_wgs84)
xx, yy = np.meshgrid(x, y)
lon, lat = transformer.transform(xx, yy)

