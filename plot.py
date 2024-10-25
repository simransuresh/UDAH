import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LatitudeFormatter, LongitudeFormatter

data = pd.read_csv('grid_dh_bg_stg2.csv')

# Step 3: Extract latitude, longitude, and dynamic height
lat = data['Latitude'].to_numpy()
lon = data['Longitude'].to_numpy()
dynamic_height = data['Dynamic_Height'].to_numpy()

# Step 4: Create a map with North Polar Stereographic projection
fig, ax = plt.subplots(figsize=(10, 10), subplot_kw={'projection': ccrs.NorthPolarStereo()})

# Step 5: Set map boundaries and limits
ax.set_extent([-160, -120, 65, 85], ccrs.PlateCarree())

# Add coastlines, gridlines, and labels
ax.coastlines()
ax.gridlines(draw_labels=True)

# Step 6: Plot the dynamic height data
scatter = ax.scatter(lon, lat, c=dynamic_height, cmap='YlGnBu', s=10, transform=ccrs.PlateCarree(), vmin=0, vmax=1.2)

# Add a color bar to represent dynamic height values
cbar = plt.colorbar(scatter, ax=ax, orientation='vertical', label='Dynamic Height [m]')

# Title and display the plot
ax.set_title("Dynamic Height of Arctic ocean from hydrographic observations 2011-2018")
plt.show()
# plt.savefig('gridded_dh_bg_20110901.png')
