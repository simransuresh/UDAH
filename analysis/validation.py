from SAGA.saga import dot, ug, vg, lats, lons
import pandas as pd
import numpy as np
from scipy.interpolate import griddata
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt


# 2012 mean from saga
print(lats.shape, dot.shape)    # target grid
lons, lats = np.meshgrid(lons, lats) 

# mapped 2012 mean
df = pd.read_csv('grd_dh_2011_01.csv')
# df = df[df['Datetime']=='2015-01-01']
X = df['X_meters'].values
Y = df['Y_meters'].values
grd_lons = df['Longitude'].values
grd_lats = df['Latitude'].values
grd_dhs = df['Surf_DH_err'].values

# Flatten the grid data for interpolation
points = np.column_stack((lats.ravel(), lons.ravel()))
dot = dot.ravel()
ug = ug.ravel()
vg = vg.ravel()

# # Step 3: Perform reverse interpolation
scattered_points = np.column_stack((grd_lats, grd_lons))
interpolated_dot = griddata(points, dot, scattered_points, method='nearest')

# # Step 4: Compute the difference between original and interpolated scattered values
diff_dot = grd_dhs - interpolated_dot
print(max(grd_dhs))

m = Basemap(projection='nplaea', boundinglat=70, lon_0=0, resolution='l', round=True)
fig, ax = plt.subplots(figsize=(8, 8))
m.drawcoastlines(linewidth=1.2)
m.drawparallels([80], labels=[True], linewidth=0.5, color='gray', fontsize=6)
m.drawmeridians([-180, -90, 0, 90], labels=[True, True, True, True], linewidth=0.5, color='gray', fontsize=6)

# print(min(difference), max(difference))
grd_lons, grd_lats = m(grd_lons, grd_lats)
sc = m.scatter(grd_lons, grd_lats, c=grd_dhs, cmap='RdYlBu', s=30)
cbar = plt.colorbar(sc, ax=ax, shrink=0.7)
cbar.set_label("Difference SAGA-UDAH (m)")

# plt.savefig('saga_udah_diff.png', dpi=300)
plt.show()