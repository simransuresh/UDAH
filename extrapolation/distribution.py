import pandas as pd
from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from mpl_toolkits.basemap import Basemap

##### preprocessing
# combined PS122_CTD.tab and PS122_Ocean_City.tab- removed spaces, added comma, replaced * with , and dropped Others column
df = pd.read_csv('profiles.csv')
# df = df.drop(columns=['Others'])
# df.to_csv('profiles.csv', index=False)
# print(df.head)


#### plot mosaic ctd profiles
df['Datetime'] = pd.to_datetime(df['Datetime'])
df['Month'] = df['Datetime'].dt.month

fig = plt.figure(figsize=(6, 6))
m = Basemap(projection='npstere', lon_0=0, lat_0=90, resolution='l', boundinglat=75, round=True)

# Draw coastlines and gridlines
m.drawcoastlines()
m.drawparallels(range(75, 91, 5), labels=[1,0,0,0])  # Parallels from 75N to 90N
m.drawmeridians(range(-180, 180, 60), labels=[0,0,1,1])

x, y = m(df['Longitude'].values, df['Latitude'].values)
sc = m.scatter(x, y, c=df['Month'], cmap='viridis', s=70, edgecolor='black', alpha=0.8)

cb = plt.colorbar(sc, orientation="vertical", label="Month", shrink=0.7)
cb.set_ticks(range(1, 13))  # Ensure all months are shown

plt.show()
# plt.savefig('mosaic_ctd_dist.png', dpi=300)


##### add to this BG profiles from whoi and plot dist


##### consolidate 2019, 2020, BG and TPD, Surf_DH 


##### mapping for 2019, 2020 months